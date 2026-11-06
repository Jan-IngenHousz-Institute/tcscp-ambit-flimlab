"""IRF reconvolution fit of a TCSPC decay: the whole decay, its rise included.

Model, on the hardware micro-time grid (N bins over one laser period T)::

    m(t) = floor + [ IRF(t - shift) (*) sum_i a_i exp(-t / tau_i) ]      (circular over T)

- the decay is periodic: the tails of the previous laser pulses are included;
- IRF: the measured one (``blank_dose.load_irf``), floor-subtracted, placed with its peak at the
  same 2.5 ns as the rotated data, normalised to a sum of 1;
- ``shift`` slides the IRF against the fluorescence (the SPAD timing depends on the wavelength,
  and the IRF was measured on scattered laser light);
- the fit compares the model with the data in 47 ps bins from ``START_NS`` before the laser peak
  to ``END_NS`` after it (the echo ~8 ns after the peak is left out), by Poisson maximum
  likelihood, so low-count bins are weighted correctly.

What is reported for each component: its lifetime; its ``peak``, the counts per 47 ps bin of the
convolved component at the peak of the decay, as drawn (the peaks plus the floor add up to the data
peak); its ``share`` of the photons (a tau over the whole decay). The amplitude before the IRF,
``a`` (counts per 47 ps bin at time zero, not normalised), is a measurable quantity only for a
*resolved* component, tau >= ``RESOLVED_FRAC`` x the IRF FWHM: for a shorter one the same photons
fit equally well with a twice as large and tau twice as short, so a comes out far above the data
and means nothing; only a tau (the photons) is determined.

Every decay is fitted with one exponential and with two. The two-exponential model is kept only
when it lowers the deviance (-2 log L) by at least ``DEV_2EXP`` (10 for two more parameters: a
chi2 with 2 degrees of freedom at p < 0.01). Otherwise the decay is a single exponential, reported
in slot 2 (the "slow" one), and slot 1 is *empty*: share 0, no lifetime, no error (``empty``),
left out of ``tau_int``, and out of the covariance, which comes from the one-exponential fit. At
the counts of a 10 ms pause (~1e4 photons) a single exponential is otherwise split into two
nearly equal lifetimes whose errors, propagated through that degeneracy, are 2-20 times the real
scatter; with the choice they match it (synthetic decays). A component of the two-exponential
fit that still ends up with (almost) no photons, share below ``EMPTY_SHARE``, is empty all the
same (its two parameters would make the covariance singular), and so is an unresolved component
below ``EMPTY_SHARE_UNRESOLVED``, which is degenerate with the IRF shift.

Two summary lifetimes are reported, both with their uncertainty from the fit covariance.
``t_avg_all`` (t_amp) = amplitude-weighted lifetime of the components with photons,
sum(a tau) / sum(a): proportional to the fluorescence yield (the same radiative rate for every
component), so phi2 is computed from it; the amplitude of an unresolved component is poorly
determined (see above), so such a component carries most of its error. ``tau_int`` (t_int) =
photon-weighted mean lifetime, sum(a tau^2) / sum(a tau): the model's counterpart of the
model-free mean arrival time, weighted towards the long components, a check of the decay shape.
``t_avg`` (amplitude-weighted over the resolved components only) is in the result for comparison.

The convolution is done in Fourier space: a sampled periodic exponential has a closed-form DFT, so
one inverse FFT per model evaluation.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

REBIN = 20  # 20 x 2.34 ps = 47 ps fit bins, as everywhere else
PEAK_AT_NS = 2.5  # the data and the IRF have their peak here
START_NS = 1.5  # fit window: from this long before the laser peak (floor + rise) ...
END_NS = 7.0  # ... to this long after it; the echo starts ~7.7 ns after the peak
SHIFT_MAX_NS = 1.0  # the shift also absorbs the IRF-to-peak delay (~-0.35 ns at 1.5 ns) and the
# one-bin jitter of the peak at ~1e4 photons: 0.5 pinned 1 fit in 3 there
TAU_BOUNDS = ((0.01, 10.0), (0.05, 30.0))
RESOLVED_FRAC = 0.5  # a component counts as resolved when tau >= this fraction of the IRF FWHM
EMPTY_SHARE = 1e-3  # a component with a smaller share of the photons is empty: not a parameter
EMPTY_SHARE_UNRESOLVED = 1e-2  # the same for an unresolved component (degenerate with the IRF shift)
DEV_2EXP = 10.0  # the second exponential is kept when it lowers the deviance by at least this


def prepare_irf(irf: dict, n_bins: int, bin_s: float, peak_at_ns: float = PEAK_AT_NS) -> dict:
    """The measured IRF on the data's micro-time grid (``n_bins`` of ``bin_s``): one period,
    peak at ``peak_at_ns``, floor (mean from 0.3 ns to 0.6 ns before the peak) subtracted,
    negative bins set to 0, sum 1. Returns the IRF and its real FFT."""
    T = n_bins * bin_s
    b = float(irf["bin_s"])
    c = np.asarray(irf["counts"], float)
    nz = np.flatnonzero(c > 0)  # the TDC-Software export pads outside its range with 0
    n = int(round(T / b))
    seg = c[nz[0] : nz[0] + n]
    k = int(np.argmax(seg))
    seg = np.roll(seg, int(round(peak_at_ns * 1e-9 / b)) - k)
    t_irf = (np.arange(seg.size) + 0.5) * b
    pre = (t_irf > 0.3e-9) & (t_irf < (peak_at_ns - 0.6) * 1e-9)
    y = np.clip(seg - seg[pre].mean(), 0.0, None)
    t = (np.arange(n_bins) + 0.5) * bin_s
    g = np.interp(t, t_irf, y, period=T)  # counts per IRF bin -> density on the data grid
    g = g / g.sum()
    return dict(g=g, G=np.fft.rfft(g), n_bins=n_bins, bin_s=bin_s, name=irf.get("name", ""), fwhm_ps=_fwhm(t, g) * 1e12)


def _fwhm(t, y) -> float:
    i = int(np.argmax(y))
    half = y[i] / 2
    lo, hi = i, i
    while lo > 0 and y[lo] > half:
        lo -= 1
    while hi < y.size - 1 and y[hi] > half:
        hi += 1
    return float(t[hi] - t[lo])


class _Model:
    """Model of one period in raw bins, then summed into 47 ps bins."""

    def __init__(self, irf_p: dict, rebin: int = REBIN):
        self.N, self.dt = irf_p["n_bins"], irf_p["bin_s"] * 1e9  # ns
        self.T = self.N * self.dt
        self.G = irf_p["G"]
        m = np.arange(self.N // 2 + 1)
        self.w = np.exp(-2j * np.pi * m / self.N)  # e^{-2 pi i m / N}
        self.f = np.fft.rfftfreq(self.N, d=self.dt)  # cycles per ns
        self.rebin = rebin
        self.nr = self.N // rebin

    def _exp_hat(self, tau):
        """DFT of the periodic exponential sampled at bin centres, per unit amplitude per raw bin."""
        r = np.exp(-self.dt / tau)
        return np.exp(-self.dt / (2 * tau)) / (1.0 - r * self.w)

    def component(self, a, tau, shift):
        """One convolved component in counts per 47 ps bin (a: counts per 47 ps bin at time zero)."""
        spec = (a / self.rebin) * self._exp_hat(tau) * self.G * np.exp(-2j * np.pi * self.f * shift)
        raw = np.fft.irfft(spec, n=self.N)
        return raw[: self.nr * self.rebin].reshape(self.nr, self.rebin).sum(1)

    def __call__(self, p):
        a1, a2, t1, t2, s, b = p
        return b + self.component(a1, t1, s) + self.component(a2, t2, s)


def _dev_residuals(m, c):
    """Signed Poisson deviance residuals: sum of squares = -2 log-likelihood ratio."""
    m = np.maximum(m, 1e-9)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = 2.0 * (m - c + np.where(c > 0, c * np.log(np.where(c > 0, c, 1.0) / m), 0.0))
    return np.sign(m - c) * np.sqrt(np.maximum(d, 0.0))


def fit(
    h_raw: np.ndarray,
    bin_s: float,
    peak_bin: int,
    irf_p: dict,
    start_ns: float = START_NS,
    end_ns: float = END_NS,
    rebin: int = REBIN,
) -> dict:
    """2-exponential reconvolution fit of one raw hardware histogram (see the module docstring for
    what is reported).

    ``h_raw``: the N raw bins as recorded; ``peak_bin``: the laser peak in 47 ps bins of the
    unrotated data (``blank_dose.decay_summary``'s ``peak_bin``), used to rotate it to 2.5 ns.
    """
    h_raw = np.asarray(h_raw, float)
    N = h_raw.size
    if N != irf_p["n_bins"]:
        raise ValueError(f"histogram has {N} bins, the IRF was prepared for {irf_p['n_bins']}")
    bw = bin_s * rebin * 1e9
    p0 = int(round(PEAK_AT_NS / bw))
    h_rot = np.roll(h_raw, (p0 - peak_bin) * rebin)
    model = _Model(irf_p, rebin)
    h = h_rot[: model.nr * rebin].reshape(model.nr, rebin).sum(1)
    t = (np.arange(model.nr) + 0.5) * bw
    win = (t >= PEAK_AT_NS - start_ns) & (t <= PEAK_AT_NS + end_ns)
    c = h[win]
    out = dict(
        ok=False,
        t=t,
        h=h,
        window=(PEAK_AT_NS - start_ns, PEAK_AT_NS + end_ns),
        fit_mask=win,
        n_counts=float(h.sum()),
        irf_name=irf_p["name"],
        irf_fwhm_ps=irf_p["fwhm_ps"],
    )
    if c.sum() < 300:
        return out
    pre = (t > 0.3) & (t < PEAK_AT_NS - 0.6)
    b0 = float(h[pre].mean())
    hp = max(float(c.max()) - b0, 1.0)

    def resid2(p):
        return _dev_residuals(model(p)[win], c)

    def resid1(q):  # single exponential: the second amplitude is 0
        return _dev_residuals(model(np.array([q[0], 0.0, q[1], 1.0, q[2], q[3]]))[win], c)

    # one exponential first: [a, tau, shift, floor]
    lo1 = np.array([0.0, TAU_BOUNDS[1][0], -SHIFT_MAX_NS, 0.0])
    hi1 = np.array([np.inf, TAU_BOUNDS[1][1], SHIFT_MAX_NS, np.inf])
    best1 = None
    for t_0 in (0.4, 1.5):
        sol = least_squares(resid1, np.array([1.1 * hp, t_0, 0.0, b0]), bounds=(lo1, hi1), x_scale="jac", max_nfev=2000)
        if best1 is None or sol.cost < best1.cost:
            best1 = sol
    # then two: [a1, a2, tau1, tau2, shift, floor]; two starts for the fast lifetime and one from the
    # single exponential with a small fast component added
    lo = np.array([0.0, 0.0, TAU_BOUNDS[0][0], TAU_BOUNDS[1][0], -SHIFT_MAX_NS, 0.0])
    hi = np.array([np.inf, np.inf, TAU_BOUNDS[0][1], TAU_BOUNDS[1][1], SHIFT_MAX_NS, np.inf])
    a_1, t_1, s_1, b_1 = best1.x
    starts = [
        np.array([1.5 * hp, 0.5 * hp, 0.15, 1.6, 0.0, b0]),
        np.array([1.5 * hp, 0.5 * hp, 0.4, 1.6, 0.0, b0]),
        np.array([0.2 * a_1, a_1, 0.15, t_1, s_1, b_1]),
    ]
    best = None
    for p0_ in starts:
        sol = least_squares(resid2, np.clip(p0_, lo, hi), bounds=(lo, hi), x_scale="jac", max_nfev=2000)
        if best is None or sol.cost < best.cost:
            best = sol
    dev1, dev2 = 2.0 * best1.cost, 2.0 * best.cost  # deviance = -2 log L (up to a constant)
    single = bool(dev1 - dev2 < DEV_2EXP)
    sc_all = np.array([1.0, 1.0, 1e-3, 1e-3, 1e-3, 1.0])
    if single:
        # the second exponential does not improve the fit: the decay is one exponential, reported in
        # slot 2 (the "slow" one); slot 1 is empty
        s, b = best1.x[2], best1.x[3]
        taus, amps = np.array([np.nan, best1.x[1]]), np.array([0.0, best1.x[0]])
        p_full = np.array([0.0, amps[1], 1.0, taus[1], s, b])
        empty = np.array([True, False])
        k = np.array([1, 3, 4, 5])  # a2 t2 s b of the sorted basis
        J, ps_k, lo_k, hi_k = best1.jac, best1.x, lo1, hi1
        n_par, cost = 4, best1.cost
    else:
        p = best.x
        a1, a2, t1, t2, s, b = p
        order = [0, 1] if t1 <= t2 else [1, 0]  # component 1 = the shorter lifetime
        taus, amps = np.array([t1, t2])[order], np.array([a1, a2])[order]
        idx = [order[0], order[1], 2 + order[0], 2 + order[1], 4, 5]  # the sorted basis [a1 a2 t1 t2 s b]
        p_full = p[idx]
        n_ph = amps * taus / bw
        share = n_ph / n_ph.sum() if n_ph.sum() > 0 else np.full(2, np.nan)
        wide = taus >= RESOLVED_FRAC * irf_p["fwhm_ps"] * 1e-3
        # a component with (almost) no photons is empty all the same: its amplitude and lifetime are
        # not parameters (their columns of the Jacobian are ~0, J^T J singular, every error nan). An
        # unresolved component with a tiny share too: it is degenerate with the IRF shift, and the
        # covariance through that degeneracy overstates every error
        empty = ~(share >= EMPTY_SHARE) | (~wide & ~(share >= EMPTY_SHARE_UNRESOLVED))
        if empty.all():
            empty[:] = False
        active = np.array([~empty[0], ~empty[1], ~empty[0], ~empty[1], True, True])
        k = np.flatnonzero(active)
        J, ps_k, lo_k, hi_k = best.jac[:, idx][:, k], p_full[k], lo[idx][k], hi[idx][k]
        n_par, cost = 6, best.cost
    dof = max(int(win.sum()) - n_par, 1)
    chi2r = float(2.0 * cost / dof)
    # covariance from the Jacobian over the active parameters, each scaled to its own size first:
    # amplitudes of ~1e5 next to lifetimes of ~0.05 ns make J^T J too ill-conditioned to invert directly
    sc = np.maximum(np.abs(ps_k), sc_all[k])
    Js = J * sc[None, :]
    cov = np.zeros((6, 6))  # 0 for the parameters of an empty component
    try:
        cov[np.ix_(k, k)] = np.linalg.inv(Js.T @ Js) * np.outer(sc, sc) * max(chi2r, 1.0)
    except np.linalg.LinAlgError:
        cov[:] = np.nan
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    err[[i for i in range(6) if i not in k]] = np.nan
    taus_err, amps_err = err[2:4], err[0:2]
    full = model(p_full)
    comps = [np.zeros(model.nr) if empty[i] else model.component(amps[i], taus[i], s) for i in range(2)]
    pk = int(np.argmax(full))
    peak_amps = np.array([cc[pk] for cc in comps])  # counts per 47 ps bin at the decay peak, as drawn
    tt = np.where(empty, 1.0, taus)  # a finite stand-in for an empty component
    n_ph = np.where(empty, 0.0, amps * tt / bw)  # photons of each component over the whole decay
    n_ph_err = np.array(
        [
            0.0
            if empty[i]
            else np.sqrt(
                max(
                    (tt[i] / bw) ** 2 * cov[i, i]
                    + (amps[i] / bw) ** 2 * cov[2 + i, 2 + i]
                    + 2 * amps[i] * tt[i] / bw**2 * cov[i, 2 + i],
                    0.0,
                )
            )
            for i in range(2)
        ]
    )
    share = n_ph / n_ph.sum() if n_ph.sum() > 0 else np.full(2, np.nan)
    resolved_ns = RESOLVED_FRAC * irf_p["fwhm_ps"] * 1e-3
    resolved = (tt >= resolved_ns) & ~empty

    def weighted(mask, power):
        """sum(a tau^power) / sum(a tau^(power-1)) over the components in ``mask``, with its error."""
        if not mask.any():
            return float("nan"), float("nan")
        w = np.where(mask, amps, 0.0)
        num, den = np.sum(w * tt**power), np.sum(w * tt ** (power - 1))
        if den <= 0:
            return float("nan"), float("nan")
        val = num / den
        g = np.zeros(6)
        g[0:2] = mask * (tt**power - val * tt ** (power - 1)) / den
        g[2:4] = mask * (power * w * tt ** (power - 1) - val * (power - 1) * w * tt ** (power - 2)) / den
        return float(val), float(np.sqrt(max(g @ cov @ g, 0.0)))

    t_avg, t_avg_err = weighted(resolved, 1)  # amplitude-weighted, resolved components only
    t_avg_all, t_avg_all_err = weighted(~empty, 1)  # the same over both (meaningless if one is unresolved)
    tau_int, tau_int_err = weighted(~empty, 2)  # photon-weighted mean lifetime, components with photons
    at_bound = np.isclose(ps_k, lo_k, atol=1e-6) | np.isclose(ps_k, hi_k, atol=1e-6)
    out.update(
        ok=bool((best1 if single else best).success),
        taus=taus.tolist(),
        taus_err=taus_err.tolist(),
        amps=amps.tolist(),
        amps_err=amps_err.tolist(),
        peak_amps=peak_amps.tolist(),
        n_ph=n_ph.tolist(),
        n_ph_err=n_ph_err.tolist(),
        share=share.tolist(),
        resolved=resolved.tolist(),
        empty=empty.tolist(),
        n_comp=int((~empty).sum()),
        single=single,
        deviance=(float(dev1), float(dev2)),
        resolved_ns=float(resolved_ns),
        t_avg=t_avg,
        t_avg_err=t_avg_err,
        t_avg_all=t_avg_all,
        t_avg_all_err=t_avg_all_err,
        tau_int=tau_int,
        tau_int_err=tau_int_err,
        shift_ns=float(s),
        shift_err_ns=float(err[4]),
        floor=float(b),
        floor_err=float(err[5]),
        chi2r=chi2r,
        n_fit=int(win.sum()),
        model=full,
        comps=comps,
        peak_bin=pk,
        cov=cov,
        at_bound=bool(at_bound.any()),
    )
    return out


def plot(decays: list, title: str = "", save_to=None, show: bool = True, ylim=None, yscale: str = "linear"):
    """One column per decay, ``decays`` = [(name, fit result, colour), ...]: raw counts per 47 ps bin,
    the model over the fit window (solid) and the rest of the period (thin), and each component
    convolved with the IRF, drawn from zero under the data, with the floor as a dashed line. At t0,
    the laser peak (dotted line), the two component heights plus the floor add up to the fit, so on
    the linear axis they can be summed by eye. A text box with lifetime, height at t0, photon share
    and, for a resolved component, the amplitude before the IRF, then t_int; Pearson residuals
    below, grey outside the fit window. ``yscale``: "linear" (from 0 to
    ``ylim[1]``, or to the data peak if None) or "log" (``ylim`` = its range; heights do not add by
    eye there, but the tails are visible)."""
    import matplotlib.pyplot as plt

    decays = [x for x in decays if x[1] is not None and x[1].get("ok")]
    if not decays:
        return None
    n = len(decays)
    fig, axes = plt.subplots(
        2, n, figsize=(6.2 * n, 8.2), squeeze=False, sharex=True, gridspec_kw=dict(height_ratios=[3, 1])
    )
    for j, (name, r, col) in enumerate(decays):
        ax, axr = axes[0, j], axes[1, j]
        t, h, m, win = r["t"], r["h"], r["model"], r["fit_mask"]
        w0, w1 = r["window"]
        ax.axvspan(0, w0, color="0.92", lw=0)
        ax.axvspan(w1, t[-1] + t[0], color="0.92", lw=0)
        axr.axvspan(0, w0, color="0.92", lw=0)
        axr.axvspan(w1, t[-1] + t[0], color="0.92", lw=0)
        # each component convolved with the IRF, from zero; at t0 (the laser peak, the top of the fit)
        # the two heights plus the floor add up to the fit
        lin = yscale != "log"
        pk = int(r["peak_bin"])
        t0 = t[pk]
        ax.axhline(r["floor"], color="0.5", lw=0.8, ls="--", label=f"floor {r['floor']:,.1f} /bin")
        for i, (cc, cl) in enumerate(zip(r["comps"], ("tab:blue", "tab:purple"))):
            if r["empty"][i]:
                continue
            if lin:
                ax.fill_between(t, 0, cc, color=cl, alpha=0.15, lw=0)
            ax.plot(
                t,
                cc,
                color=cl,
                lw=1.2,
                label=f"component {i + 1}: tau {r['taus'][i]:.3f} ns, {r['peak_amps'][i]:,.0f} /bin at t0, "
                f"{100 * r['share'][i]:.0f} % of photons",
            )
            ax.plot(t0, cc[pk], "o", color=cl, mec="k", mew=0.6, ms=5, zorder=5)
        ax.step(t, h if lin else np.maximum(h, 0.5), where="mid", color=col, lw=0.8, label="data (47 ps bins)")
        ax.plot(t, m, color="k", lw=0.6, alpha=0.5)
        ax.plot(
            t[win], m[win], color="k", lw=1.3, label=f"fit: t_amp {r['t_avg_all']:.3f} ns, t_int {r['tau_int']:.3f} ns"
        )
        ax.plot(t0, m[pk], "o", color="k", ms=5, zorder=5)
        ax.axvline(t0, color="k", ls=":", lw=0.8)
        txt = [f"{h[win].sum():,.0f} counts in the window, chi2r {r['chi2r']:.2f}"]
        for i in range(2):
            if r["empty"][i]:
                txt.append(f"{i + 1}: none (no photons: the decay is a single exponential)")
                continue
            tag = (
                f"a {r['amps'][i]:,.0f}+-{r['amps_err'][i]:,.0f} /bin before the IRF"
                if r["resolved"][i]
                else f"unresolved (< {r['resolved_ns']:.2f} ns): no amplitude"
            )
            txt.append(
                f"{i + 1}: tau {r['taus'][i]:.3f}+-{r['taus_err'][i]:.3f} ns   {r['peak_amps'][i]:,.0f} /bin "
                f"at t0   {100 * r['share'][i]:.0f} % of photons"
            )
            txt.append(f"   {tag}")
        txt += [
            f"at t0 = {t0:.2f} ns: {r['peak_amps'][0]:,.0f} + {r['peak_amps'][1]:,.0f} + floor {r['floor']:,.0f}"
            f" = fit {m[pk]:,.0f}  (data {h[pk]:,.0f})",
            f"t_amp (amplitude-weighted lifetime) = {r['t_avg_all']:.3f} +- {r['t_avg_all_err']:.3f} ns",
            f"t_int (photon-weighted mean lifetime) = {r['tau_int']:.3f} +- {r['tau_int_err']:.3f} ns",
            f"IRF shift {r['shift_ns'] * 1e3:+.1f} +- {r['shift_err_ns'] * 1e3:.1f} ps"
            + ("   (a parameter at its bound)" if r["at_bound"] else ""),
        ]
        # the numbers go above the panel, so that they never cover the data (fixed count axis)
        ax.text(
            0.5,
            1.01,
            "\n".join(txt),
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=7.5,
            family="monospace",
            multialignment="left",
            bbox=dict(fc="white", ec="0.8", alpha=0.9),
        )
        if yscale == "log":
            ax.set_yscale("log")
            ax.set_ylim(*(ylim or (max(0.5, 0.3 * r["floor"]), 2.0 * h.max())))
        else:
            ax.set_ylim(0, ylim[1] if ylim else 1.12 * max(h[win].max(), m[pk]))
        ax.set_title(f"{title}{': ' if title else ''}{name}", fontsize=10, pad=10 + 10.5 * len(txt))
        ax.set_ylabel("counts per 47 ps bin")
        ax.legend(fontsize=7.5, loc="upper right" if lin else "lower left")
        pr = (h - m) / np.sqrt(np.maximum(m, 1.0))
        axr.plot(t[win], pr[win], color="k", lw=0.7)
        axr.plot(t, np.where(win, np.nan, pr), color="0.6", lw=0.6)
        axr.axhline(0, color="0.5", lw=0.6)
        axr.set_ylim(-6, 6)
        axr.set_xlim(0, t[-1] + t[0])
        axr.set_xlabel("micro-time [ns] (laser peak at 2.5 ns; grey: outside the fit)")
        axr.set_ylabel("residual [sigma]")
        for x in (ax, axr):
            x.grid(alpha=0.3)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=120)
    if show:
        plt.show()
    return fig
