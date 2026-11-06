"""Multi-exponential decay fitting: iterative reconvolution with an IRF, or tail fit.

All times are in nanoseconds. Histogram bins are assumed to be equally spaced.

The reconvolution model is evaluated in the frequency domain: the spectrum of a
periodic exponential is known in closed form, and a fractional IRF shift is a
phase ramp, so one inverse FFT per evaluation gives the full model. Long
histograms are rebinned to at most MAX_FIT_BINS before fitting so that a live
fit stays well under a second.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

MAX_FIT_BINS = 2048


@dataclass
class FitResult:
    n_exp: int
    amps: np.ndarray  # amplitude of each component at the pulse (counts/bin, data grid)
    taus_ns: np.ndarray
    bg: float  # flat background, counts per data bin
    shift_ns: float  # IRF shift (reconvolution only)
    chi2r: float
    model: np.ndarray  # model on the data grid (nan where not fitted)
    mask: np.ndarray  # data bins covered by the fit
    tau_amp_ns: float  # amplitude-weighted mean lifetime
    tau_int_ns: float  # intensity-weighted mean lifetime
    mode: str  # "reconv", "tail" or "none"
    success: bool
    message: str = ""
    params: np.ndarray = field(default_factory=lambda: np.zeros(0))
    resid: np.ndarray = field(default_factory=lambda: np.zeros(0))  # weighted residuals (fit grid)
    resid_t_ns: np.ndarray = field(default_factory=lambda: np.zeros(0))  # their time axis
    rebin: int = 1

    @property
    def fractions(self) -> np.ndarray:
        s = np.nansum(self.amps)
        return self.amps / s if s > 0 else np.zeros_like(self.amps)


def average_lifetimes(amps, taus) -> tuple[float, float]:
    amps = np.asarray(amps, float)
    taus = np.asarray(taus, float)
    sa = amps.sum()
    sat = (amps * taus).sum()
    tau_amp = sat / sa if sa > 0 else float("nan")
    tau_int = (amps * taus**2).sum() / sat if sat > 0 else float("nan")
    return float(tau_amp), float(tau_int)


def periodic_decay(t_ns: np.ndarray, tau_ns: float, period_ns: float) -> np.ndarray:
    """Steady-state response to a pulse train: exp(-t/tau) wrapped over one period."""
    x = np.exp(-t_ns / tau_ns)
    denom = 1.0 - np.exp(-period_ns / tau_ns)
    return x / denom if denom > 1e-12 else x


def rebin(y: np.ndarray, factor: int) -> np.ndarray:
    """Sum consecutive groups of `factor` bins (the tail shorter than factor is dropped)."""
    y = np.asarray(y, float)
    if factor <= 1:
        return y.copy()
    m = (y.size // factor) * factor
    return y[:m].reshape(-1, factor).sum(axis=1)


class ReconvModel:
    """Sum of exponentials convolved (circularly, over the laser period) with the IRF.

    The histogram window is taken to be the first n bins of a grid that spans one
    full laser period, so decay tails that leave the window on the right come back
    in on the left, as they do in the real measurement.
    """

    def __init__(self, irf: np.ndarray, bin_ns: float, period_ns: float | None):
        irf = np.clip(np.asarray(irf, float), 0, None)
        self.n = irf.size
        self.bin_ns = bin_ns
        m = self.n
        if period_ns is not None and period_ns > self.n * bin_ns:
            m = int(round(period_ns / bin_ns))
        self.m = m
        self.period_ns = m * bin_ns
        pad = np.zeros(m)
        pad[: self.n] = irf
        s = pad.sum()
        self.irf = pad / s if s > 0 else pad
        self.t = (np.arange(m) + 0.5) * bin_ns
        self._firf = np.fft.rfft(self.irf)
        self._k = np.arange(self._firf.size)
        self._omega = np.exp(-2j * np.pi * self._k / m)

    def spectrum(self, shift_ns: float, amps, taus) -> np.ndarray:
        phase = np.exp(-2j * np.pi * self._k * (shift_ns / self.bin_ns) / self.m)
        spec = np.zeros_like(self._firf)
        for a, tau in zip(amps, taus):
            r = np.exp(-self.bin_ns / tau)
            spec += a * np.exp(-0.5 * self.bin_ns / tau) / (1.0 - r * self._omega)
        return self._firf * phase * spec

    def __call__(self, bg: float, shift_ns: float, amps, taus) -> np.ndarray:
        out = np.fft.irfft(self.spectrum(shift_ns, amps, taus), n=self.m)
        return bg + out[: self.n]


def _initial_taus(n_exp: int, bin_ns: float, span_ns: float) -> np.ndarray:
    lo = max(0.15, 3 * bin_ns)
    hi = max(lo * 3, span_ns / 4)
    if n_exp == 1:
        return np.array([np.sqrt(lo * hi)])
    return np.geomspace(lo, hi, n_exp)


def fit_decay(
    counts: np.ndarray,
    bin_ns: float,
    n_exp: int = 1,
    irf: np.ndarray | None = None,
    period_ns: float | None = None,
    p0: np.ndarray | None = None,
    max_nfev: int = 200,
    max_bins: int = MAX_FIT_BINS,
) -> FitResult:
    """Fit one histogram with n_exp exponentials.

    With an IRF (same binning as counts) the model is reconvolved with it.
    Without one, a tail fit is done starting shortly after the peak.
    Results (amplitudes, background, model) are expressed on the original grid.
    """
    counts0 = np.asarray(counts, float)
    n0 = counts0.size
    n_exp = int(max(1, n_exp))
    if n0 < 8 or counts0.sum() < 50:
        return _empty(n_exp, counts0, "not enough counts")

    factor = max(1, int(np.ceil(n0 / max_bins)))
    counts = rebin(counts0, factor)
    bw = bin_ns * factor
    n = counts.size
    t = (np.arange(n) + 0.5) * bw
    w = 1.0 / np.sqrt(np.maximum(counts, 1.0))
    peak = int(np.argmax(counts))
    bg0 = float(np.percentile(counts, 5))

    irf_r = None
    if irf is not None and np.asarray(irf).size == n0 and np.nansum(irf) > 0:
        irf_r = rebin(np.nan_to_num(np.asarray(irf, float)), factor)

    if irf_r is not None:
        model = ReconvModel(irf_r, bw, period_ns)
        mask_c = np.ones(n, bool)
        mode = "reconv"
        npar = 2 + 2 * n_exp

        def unpack(p):
            return p[0], p[1], p[2 : 2 + n_exp], p[2 + n_exp :]

        def resid(p):
            bg, sh, a, tau = unpack(p)
            return (model(bg, sh, a, tau) - counts) * w

        if p0 is None or p0.size != npar:
            taus0 = _initial_taus(n_exp, bw, n * bw)
            excess = max(counts.sum() - bg0 * n, 1.0)
            amps0 = excess / (taus0 / bw) / n_exp
            p0 = np.concatenate([[bg0, 0.0], amps0, taus0])
        shift_max = max(5 * bw, 0.5)
        lo = np.concatenate([[0.0, -shift_max], np.zeros(n_exp), np.full(n_exp, 0.3 * bw)])
        hi = np.concatenate([[np.inf, shift_max], np.full(n_exp, np.inf), np.full(n_exp, 10 * model.period_ns)])
        p0 = np.clip(p0, lo + 1e-12, np.where(np.isinf(hi), np.inf, hi - 1e-12))
        try:
            sol = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev)
        except Exception as exc:  # numerical trouble: report and keep the GUI alive
            return _empty(n_exp, counts0, f"fit error: {exc}")
        bg, sh, a, tau = unpack(sol.x)
        m_c = model(bg, sh, a, tau)
        order = np.argsort(tau)
        a, tau = a[order], tau[order]
        amps_peak = a.copy()
    else:
        mode = "tail"
        start = min(peak + max(3, int(round(0.25 / bw))), n - 4)
        mask_c = np.zeros(n, bool)
        mask_c[start:] = True
        tt = t[mask_c] - t[start]
        c = counts[mask_c]
        ww = w[mask_c]
        npar = 1 + 2 * n_exp

        def unpack(p):
            return p[0], p[1 : 1 + n_exp], p[1 + n_exp :]

        def tail_model(p, tt_):
            bg, a, tau = unpack(p)
            return bg + sum(ai * np.exp(-tt_ / ti) for ai, ti in zip(a, tau))

        def resid(p):
            return (tail_model(p, tt) - c) * ww

        if p0 is None or p0.size != npar:
            taus0 = _initial_taus(n_exp, bw, tt[-1])
            amps0 = np.full(n_exp, max(counts[start] - bg0, 1.0) / n_exp)
            p0 = np.concatenate([[bg0], amps0, taus0])
        lo = np.concatenate([[0.0], np.zeros(n_exp), np.full(n_exp, 0.3 * bw)])
        hi = np.concatenate([[np.inf], np.full(n_exp, np.inf), np.full(n_exp, 10 * tt[-1] + 1)])
        p0 = np.clip(p0, lo + 1e-12, np.where(np.isinf(hi), np.inf, hi - 1e-12))
        try:
            sol = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev)
        except Exception as exc:
            return _empty(n_exp, counts0, f"fit error: {exc}")
        bg, a, tau = unpack(sol.x)
        order = np.argsort(tau)
        a, tau = a[order], tau[order]
        sh = 0.0
        m_c = np.full(n, np.nan)
        m_c[mask_c] = tail_model(sol.x, tt)
        amps_peak = a * np.exp((t[start] - t[peak]) / tau)  # referred back to the peak

    dof = max(int(mask_c.sum()) - npar, 1)
    chi2r = float(np.sum(sol.fun**2) / dof)
    tau_amp, tau_int = average_lifetimes(amps_peak, tau)

    # back to the original grid: counts per fine bin
    model_f = np.full(n0, np.nan)
    model_f[: n * factor] = np.repeat(m_c / factor, factor)
    mask_f = np.zeros(n0, bool)
    mask_f[: n * factor] = np.repeat(mask_c, factor)
    return FitResult(
        n_exp=n_exp,
        amps=amps_peak / factor,
        taus_ns=tau,
        bg=float(bg) / factor,
        shift_ns=float(sh),
        chi2r=chi2r,
        model=model_f,
        mask=mask_f,
        tau_amp_ns=tau_amp,
        tau_int_ns=tau_int,
        mode=mode,
        success=bool(sol.success),
        message=str(sol.message),
        params=sol.x,
        resid=sol.fun,
        resid_t_ns=t[mask_c],
        rebin=factor,
    )


def _empty(n_exp: int, counts: np.ndarray, msg: str) -> FitResult:
    nan = float("nan")
    return FitResult(
        n_exp=n_exp,
        amps=np.full(n_exp, nan),
        taus_ns=np.full(n_exp, nan),
        bg=nan,
        shift_ns=nan,
        chi2r=nan,
        model=np.full(counts.size, np.nan),
        mask=np.zeros(counts.size, bool),
        tau_amp_ns=nan,
        tau_int_ns=nan,
        mode="none",
        success=False,
        message=msg,
    )
