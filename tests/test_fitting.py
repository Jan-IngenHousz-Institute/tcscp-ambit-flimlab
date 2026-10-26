"""Fit synthetic decays and check the recovered lifetimes.

Run with:  python -m pytest tests  (or simply  python tests/test_fitting.py)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from flimgui.fitting import ReconvModel, average_lifetimes, fit_decay  # noqa: E402
from flimgui.irf import IRF  # noqa: E402

RNG = np.random.default_rng(1)


def make_irf(n: int, bin_ns: float, pos_ns: float, sigma_ns: float, counts: float = 2e5) -> np.ndarray:
    t = (np.arange(n) + 0.5) * bin_ns
    g = np.exp(-0.5 * ((t - pos_ns) / sigma_ns) ** 2)
    g = g / g.sum() * counts
    return RNG.poisson(g).astype(float)


def make_decay(irf, bin_ns, amps, taus, period_ns, bg=2.0, total=2e6):
    model = ReconvModel(irf, bin_ns, period_ns)
    pure = model(0.0, 0.0, np.asarray(amps, float), np.asarray(taus, float))
    pure *= total / pure.sum()
    return RNG.poisson(pure + bg).astype(float)


def test_average_lifetimes():
    ta, ti = average_lifetimes([1.0, 1.0], [1.0, 3.0])
    assert abs(ta - 2.0) < 1e-12
    assert abs(ti - (1 + 9) / 4) < 1e-12


def test_reconv_single_exp():
    n, bw = 2000, 0.025                      # 50 ns window, 25 ps bins
    irf = make_irf(n, bw, 8.0, 0.12)
    data = make_decay(irf, bw, [1.0], [3.6], period_ns=50.0)
    r = fit_decay(data, bw, 1, irf=irf, period_ns=50.0)
    assert r.success
    assert abs(r.taus_ns[0] - 3.6) / 3.6 < 0.02, r.taus_ns
    assert r.chi2r < 1.5, r.chi2r


def test_reconv_double_exp():
    n, bw = 2000, 0.025
    irf = make_irf(n, bw, 8.0, 0.12)
    amps = np.array([0.3 / 0.8, 0.7 / 3.6])    # fractions 30 % / 70 % in amplitude... scaled below
    amps = np.array([3.0, 7.0])
    taus = np.array([0.8, 3.6])
    data = make_decay(irf, bw, amps, taus, period_ns=50.0, total=5e6)
    r = fit_decay(data, bw, 2, irf=irf, period_ns=50.0)
    assert r.success
    assert abs(r.taus_ns[0] - 0.8) / 0.8 < 0.10, r.taus_ns
    assert abs(r.taus_ns[1] - 3.6) / 3.6 < 0.05, r.taus_ns
    ta_true, _ = average_lifetimes(amps, taus)
    assert abs(r.tau_amp_ns - ta_true) / ta_true < 0.05, (r.tau_amp_ns, ta_true)
    assert abs(r.fractions[1] - 0.7) < 0.08, r.fractions


def test_reconv_window_shorter_than_period():
    """25 ns window at 20 MHz: the tail wraps around to the start of the window."""
    n, bw = 1000, 0.025
    period = 50.0
    full_irf = make_irf(2000, bw, 5.0, 0.12)
    full = make_decay(full_irf, bw, [1.0], [6.0], period_ns=period, total=4e6)
    data, irf = full[:n], full_irf[:n]
    r = fit_decay(data, bw, 1, irf=irf, period_ns=period)
    assert r.success
    assert abs(r.taus_ns[0] - 6.0) / 6.0 < 0.03, r.taus_ns


def test_tail_fit():
    n, bw = 2000, 0.025
    irf = make_irf(n, bw, 8.0, 0.12)
    data = make_decay(irf, bw, [1.0], [3.6], period_ns=50.0)
    r = fit_decay(data, bw, 1, irf=None)
    assert r.success and r.mode == "tail"
    assert abs(r.taus_ns[0] - 3.6) / 3.6 < 0.03, r.taus_ns


def test_speed_and_rebin_on_long_histogram():
    """5333 bins of 2.34 ps like the Felix at 12.5 ns: a live fit must stay fast."""
    import time
    n, bw = 5333, 0.00234375
    irf = make_irf(n, bw, 2.0, 0.08)
    data = make_decay(irf, bw, [2.0, 1.0], [0.4, 1.0], period_ns=12.5, total=3e6)
    t0 = time.perf_counter()
    r = fit_decay(data, bw, 2, irf=irf, period_ns=12.5)
    dt = time.perf_counter() - t0
    print(f"   5333-bin reconv fit: {dt*1000:.0f} ms, rebin x{r.rebin}, taus {r.taus_ns.round(3)}")
    assert r.success and r.rebin == 3
    assert abs(r.taus_ns[0] - 0.4) / 0.4 < 0.10 and abs(r.taus_ns[1] - 1.0) < 0.10, r.taus_ns
    assert r.model.size == n and r.mask.size == n and r.resid.size == r.resid_t_ns.size
    assert dt < 3.0, dt
    # amplitudes/background are per original bin: model integral matches the data
    assert abs(np.nansum(r.model) - data.sum()) / data.sum() < 0.02


def test_not_enough_counts():
    r = fit_decay(np.zeros(100), 0.05, 1)
    assert not r.success and r.mode == "none"


def test_irf_resample_keeps_shape():
    irf = IRF("x", make_irf(1000, 0.05, 8.0, 0.12).astype(np.int64), left_ns=0.0, bin_ns=0.05)
    same, exact = irf.on_axis(0.0, 0.05, 1000)
    assert exact
    coarse, exact2 = irf.on_axis(0.0, 0.1, 500)
    assert not exact2
    assert abs(np.argmax(coarse) * 0.1 - np.argmax(same) * 0.05) < 0.15


if __name__ == "__main__":
    import inspect
    failures = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and inspect.isfunction(fn):
            try:
                fn()
                print("PASS", name)
            except AssertionError as exc:
                failures += 1
                print("FAIL", name, exc)
    sys.exit(1 if failures else 0)
