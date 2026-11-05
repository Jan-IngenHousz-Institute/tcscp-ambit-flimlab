"""Fv/Fm from the AMBIT and from the TCSPC lifetime, measured at the same time.

The AMBIT runs its usual Fv/Fm protocol (dark, saturating pulse, dark) while the Felix
records short lifetime windows. The SPAD must stay linear during the pulse (ND filter in
front of it, about 3 Mcps or less); the laser decay then rides on a flat floor of
pulse-excited fluorescence, which is subtracted here before fitting.

Three estimates of Fv/Fm are compared:
  AMBIT      (Fm - Fo) / Fm from the s_630 / r_630 ratio
  intensity  1 - Io / Im from the laser-correlated Felix counts per second
  lifetime   1 - tau_o / tau_m, with the amplitude-weighted lifetime of a 2-exponential
             tail fit (the fluorescence yield is proportional to it when the excitation
             is constant) and, as a cruder check, with the 1-exponential lifetime

Use:
    python fvfm_lifetime.py                      # run the protocol, analyse, save
    python fvfm_lifetime.py --run runs/<dir>     # analyse a saved run
    python fvfm_lifetime.py --port COM13 --window 0.1 --fm-last 0.2 --fo-last 1.0
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambit_gate as ag  # noqa: E402

FVFM_PROTOCOL = [
    dict(type=2, far_red=0, pulses=30, freq=10, actinic=0, subsampling=0),    # 3 s dark
    dict(type=2, far_red=0, pulses=80, freq=100, actinic=255, subsampling=0),  # 0.8 s saturating pulse
    dict(type=2, far_red=0, pulses=30, freq=100, actinic=0, subsampling=0),   # 0.3 s dark
]
REBIN = 20            # 20 x 2.34 ps = 47 ps bins for the fits
FIT_START_NS = 0.25   # tail fit starts this long after the peak
TAU_BOUNDS = (0.05, 10.0)
MIN_RESOLVED_NS = 0.15   # a tail fit starting 0.25 ns after the peak cannot resolve faster components


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def load_run(run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    z = np.load(run_dir / "felix_windows.npz")
    meta = json.loads((run_dir / "run.json").read_text())
    amb = list(csv.DictReader(open(run_dir / "ambit.csv")))
    _, edges = ag.time_axis(meta["protocol"])
    flash = [(a, b) for a, b, s in edges if s["actinic"] > ag.ACTINIC_ON]
    fm = meta.get("flash_measured_s") or list(flash[0])
    return dict(
        dir=run_dir, bins=z["bins"].astype(float), t0=z["t_start_s"], phase=z["phase"].astype(str),
        bin_s=float(z["bin_s"]), window_s=float(z["window_s"]), laser_hz=float(meta["laser_hz"]),
        amb_t=np.array([float(r["t_s"]) for r in amb]),
        amb_ratio=np.array([float(r["ratio"]) for r in amb]),
        amb_actinic=np.array([int(r["actinic"]) for r in amb]),
        flash_onset=float(fm[0]), flash_end=float(fm[1]), protocol=meta["protocol"], meta=meta,
    )


# ---------------------------------------------------------------------------
# one histogram: floor, laser counts, tail fits
# ---------------------------------------------------------------------------
def _rebin(h: np.ndarray, k: int) -> np.ndarray:
    n = (h.size // k) * k
    return h[:n].reshape(-1, k).sum(1)


def _rebin_matrix(m: np.ndarray, k: int) -> np.ndarray:
    n = (m.shape[1] // k) * k
    return m[:, :n].reshape(m.shape[0], -1, k).sum(2)


def _tail_fit(t: np.ndarray, h: np.ndarray, floor: float, start: int, peak: int, n_exp: int):
    tt = t[start:] - t[start]
    c = h[start:]
    w = 1.0 / np.sqrt(np.maximum(c, 1.0))          # Poisson weights of the raw counts
    taus0 = {1: [1.0], 2: [0.3, 1.5]}[n_exp]
    a0 = [max(c[0] - floor, 1.0) / n_exp] * n_exp
    p0 = np.array(a0 + taus0, float)
    lo = np.array([0.0] * n_exp + [TAU_BOUNDS[0]] * n_exp)
    hi = np.array([np.inf] * n_exp + [TAU_BOUNDS[1]] * n_exp)

    def model(p):
        return floor + sum(p[i] * np.exp(-tt / p[n_exp + i]) for i in range(n_exp))

    sol = least_squares(lambda p: (model(p) - c) * w, p0, bounds=(lo, hi), x_scale="jac", max_nfev=400)
    a, tau = sol.x[:n_exp], sol.x[n_exp:]
    order = np.argsort(tau)
    a, tau = a[order], tau[order]
    a_peak = a * np.exp((t[start] - t[peak]) / tau)   # amplitudes referred back to the peak
    chi2r = float(np.sum(sol.fun ** 2) / max(c.size - 2 * n_exp, 1))
    tau_amp = float(np.sum(a_peak * tau) / np.sum(a_peak)) if a_peak.sum() > 0 else float("nan")
    tau_int = float(np.sum(a_peak * tau ** 2) / np.sum(a_peak * tau)) if (a_peak * tau).sum() > 0 else float("nan")
    return dict(taus=tau, amps=a_peak, fractions=a_peak / a_peak.sum() if a_peak.sum() > 0 else a_peak,
                tau_amp=tau_amp, tau_int=tau_int, chi2r=chi2r, ok=bool(sol.success), model=model(sol.x), t_fit=t[start:])


def fit_window(bins: np.ndarray, bin_s: float, rebin: int = REBIN) -> dict:
    """Floor from the pre-peak region, laser counts above it, 1- and 2-exponential tail fits."""
    h = _rebin(np.asarray(bins, float), rebin)
    bw = bin_s * rebin * 1e9
    t = (np.arange(h.size) + 0.5) * bw
    peak = int(np.argmax(h))
    pre = (t > 0.3) & (t < t[peak] - 0.6)
    if pre.sum() < 10 or h.sum() < 200:
        return dict(ok=False, n_counts=float(h.sum()))
    floor = float(h[pre].mean())
    floor_err = float(h[pre].std(ddof=1) / np.sqrt(pre.sum()))
    rise = max(peak - 3, 0)
    laser = float((h[rise:] - floor).sum())
    laser_err = float(np.sqrt(h[rise:].sum() + (floor_err * (h.size - rise)) ** 2))
    sel = (t >= t[peak]) & (t < t[peak] + 5.0)
    wts = np.clip(h[sel] - floor, 0, None)
    t_mean = float(((t[sel] - t[peak]) * wts).sum() / wts.sum()) if wts.sum() > 0 else float("nan")
    start = min(peak + max(3, int(round(FIT_START_NS / bw))), h.size - 6)
    f1 = _tail_fit(t, h, floor, start, peak, 1)
    f2 = _tail_fit(t, h, floor, start, peak, 2)
    # a component faster than the tail fit can resolve would get an arbitrary amplitude and
    # wreck the amplitude-weighted mean, so fall back to the 1-exp lifetime in that case
    degenerate = bool(f2["taus"][0] < MIN_RESOLVED_NS) or not np.isfinite(f2["tau_amp"])
    tau_amp2 = float(f1["taus"][0]) if degenerate else f2["tau_amp"]
    return dict(ok=True, n_counts=float(h.sum()), floor=floor, floor_err=floor_err, floor_per_s=None,
                laser_counts=laser, laser_err=laser_err, peak_ns=float(t[peak]), t_mean=t_mean,
                tau1=float(f1["taus"][0]), chi2_1=f1["chi2r"],
                taus2=f2["taus"], fractions2=f2["fractions"], tau_amp2=tau_amp2, tau_int2=f2["tau_int"],
                chi2_2=f2["chi2r"], fit2_ok=f2["ok"], fit2_degenerate=degenerate, t=t, h=h, fit1=f1, fit2=f2)


# ---------------------------------------------------------------------------
# Fv/Fm
# ---------------------------------------------------------------------------
def fvfm(run: dict, fm_last_s: float = 0.2, fo_last_s: float = 1.0, log=print) -> dict:
    W = run["window_s"]
    bins, t0, phase = run["bins"], run["t0"], run["phase"]
    onset, end = run["flash_onset"], run["flash_end"]

    # AMBIT ------------------------------------------------------------------
    at, ar, act = run["amb_t"], run["amb_ratio"], run["amb_actinic"]
    fo_sel = (act == 0) & (at >= onset - fo_last_s) & (at < onset)
    pulse = (act > 0) & (at >= onset) & (at <= end + 0.02)
    fm_sel = pulse & (at >= end - fm_last_s)
    Fo_A = float(ar[fo_sel].mean()); Fm_A_max = float(ar[pulse].max()); Fm_A = float(ar[fm_sel].mean())
    ambit = dict(Fo=Fo_A, Fm_last=Fm_A, Fm_max=Fm_A_max, FvFm=(Fm_A - Fo_A) / Fm_A,
                 FvFm_max=(Fm_A_max - Fo_A) / Fm_A_max, n_fo=int(fo_sel.sum()), n_fm=int(fm_sel.sum()))

    # Felix windows -----------------------------------------------------------
    pre_k = [k for k in range(len(bins)) if phase[k] == "pre" and t0[k] + W > onset - fo_last_s - 1e-6]
    flash_k = [k for k in range(len(bins)) if phase[k] == "flash"]
    n_fm = max(1, int(round(fm_last_s / W)))
    fm_k = flash_k[-n_fm:]
    if not pre_k or not flash_k:
        raise RuntimeError("no pre or flash windows (phases: %s)" % sorted(set(phase)))
    per_window = []
    for k in range(len(bins)):
        if phase[k] == "outside":
            continue
        f = fit_window(bins[k], run["bin_s"])
        row = dict(k=int(k), t_mid=float(t0[k] + W / 2), phase=str(phase[k]), total=float(bins[k].sum()))
        if f["ok"]:
            row.update(laser_per_s=f["laser_counts"] / W, floor_per_s=f["floor"] * f["h"].size / W,
                       tau1=f["tau1"], tau_amp2=f["tau_amp2"], tau_int2=f["tau_int2"],
                       taus2=[float(x) for x in f["taus2"]], chi2_2=f["chi2_2"], t_mean=f["t_mean"])
        per_window.append(row)

    def summed(ks):
        f = fit_window(bins[ks].sum(0), run["bin_s"])
        f["laser_per_s"] = f["laser_counts"] / (len(ks) * W)
        f["laser_per_s_err"] = f["laser_err"] / (len(ks) * W)
        f["floor_per_s"] = f["floor"] * f["h"].size / W
        f["windows"] = [int(k) for k in ks]
        f["seconds"] = len(ks) * W
        return f

    Fo = summed(pre_k); Fm = summed(fm_k)
    intensity = dict(Io=Fo["laser_per_s"], Im=Fm["laser_per_s"], FvFm=1 - Fo["laser_per_s"] / Fm["laser_per_s"],
                     Io_err=Fo["laser_per_s_err"], Im_err=Fm["laser_per_s_err"])
    lifetime = dict(tau1_o=Fo["tau1"], tau1_m=Fm["tau1"], FvFm_tau1=1 - Fo["tau1"] / Fm["tau1"],
                    tau_amp_o=Fo["tau_amp2"], tau_amp_m=Fm["tau_amp2"], FvFm_tau_amp=1 - Fo["tau_amp2"] / Fm["tau_amp2"],
                    tau_int_o=Fo["tau_int2"], tau_int_m=Fm["tau_int2"],
                    taus2_o=[float(x) for x in Fo["taus2"]], fractions2_o=[float(x) for x in Fo["fractions2"]],
                    taus2_m=[float(x) for x in Fm["taus2"]], fractions2_m=[float(x) for x in Fm["fractions2"]],
                    chi2_o=Fo["chi2_2"], chi2_m=Fm["chi2_2"], t_mean_o=Fo["t_mean"], t_mean_m=Fm["t_mean"])
    res = dict(run=str(run["dir"]), window_s=W, fm_last_s=fm_last_s, fo_last_s=fo_last_s,
               flash_onset=onset, flash_end=end, fo_windows=Fo["windows"], fm_windows=Fm["windows"],
               fo_seconds=Fo["seconds"], fm_seconds=Fm["seconds"],
               spad_rate_dark=Fo["n_counts"] / Fo["seconds"], spad_rate_flash=Fm["n_counts"] / Fm["seconds"],
               ambit=ambit, intensity=intensity, lifetime=lifetime, per_window=per_window,
               _Fo=Fo, _Fm=Fm)
    log(report(res))
    return res


def report(r: dict) -> str:
    a, i, L = r["ambit"], r["intensity"], r["lifetime"]
    lines = [
        f"Fo from {r['fo_seconds']:.1f} s of dark windows, Fm from the last {r['fm_seconds']:.1f} s of the pulse "
        f"(SPAD {r['spad_rate_dark']/1e3:.0f} kcps dark, {r['spad_rate_flash']/1e6:.2f} Mcps in the pulse)",
        f"{'':>22} {'Fo':>10} {'Fm':>10} {'Fv/Fm':>7}",
        f"{'AMBIT ratio':>22} {a['Fo']:10.4f} {a['Fm_last']:10.4f} {a['FvFm']:7.3f}   (Fm = max during pulse: {a['Fm_max']:.4f} -> {a['FvFm_max']:.3f})",
        f"{'Felix laser counts/s':>22} {i['Io']:10.0f} {i['Im']:10.0f} {i['FvFm']:7.3f}   (+-{i['Io_err']:.0f}, +-{i['Im_err']:.0f}; no dead-time correction)",
        f"{'1-exp lifetime [ns]':>22} {L['tau1_o']:10.3f} {L['tau1_m']:10.3f} {L['FvFm_tau1']:7.3f}",
        f"{'2-exp amp-weighted [ns]':>22} {L['tau_amp_o']:10.3f} {L['tau_amp_m']:10.3f} {L['FvFm_tau_amp']:7.3f}   "
        f"(Fo: {np.round(L['taus2_o'],2)} ns {np.round(L['fractions2_o'],2)}, chi2r {L['chi2_o']:.2f}; "
        f"Fm: {np.round(L['taus2_m'],2)} ns {np.round(L['fractions2_m'],2)}, chi2r {L['chi2_m']:.2f})",
        f"{'2-exp int-weighted [ns]':>22} {L['tau_int_o']:10.3f} {L['tau_int_m']:10.3f}",
        f"{'mean time after peak':>22} {L['t_mean_o']:10.3f} {L['t_mean_m']:10.3f}",
    ]
    return "\n".join(lines)


def save(res: dict, run_dir: Path) -> Path:
    out = {k: v for k, v in res.items() if not k.startswith("_")}
    p = Path(run_dir) / "fvfm.json"
    p.write_text(json.dumps(out, indent=1, default=float))
    return p


def plot(res: dict, run: dict, save_to: Path | None = None, show: bool = False):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a, i, L = res["ambit"], res["intensity"], res["lifetime"]
    Fo, Fm = res["_Fo"], res["_Fm"]
    fig, axs = plt.subplots(2, 2, figsize=(14, 9))
    ax = axs[0, 0]
    ax.plot(run["amb_t"], run["amb_ratio"], ".-", ms=3, color="tab:blue")
    ax.axvspan(res["flash_onset"], res["flash_end"], color="orange", alpha=0.2)
    ax.axhline(a["Fo"], color="k", lw=0.7, ls="--"); ax.axhline(a["Fm_last"], color="k", lw=0.7, ls="--")
    ax.set_xlabel("time since AMBIT sample 0 [s]"); ax.set_ylabel("AMBIT s_630 / r_630")
    ax.set_title(f"AMBIT: Fo {a['Fo']:.4f}, Fm {a['Fm_last']:.4f}, Fv/Fm {a['FvFm']:.3f}")
    ax = axs[0, 1]
    pw = [r for r in res["per_window"] if "tau1" in r]
    col = {"pre": "tab:blue", "flash": "tab:red", "post": "tab:green", "edge": "grey"}
    for r in pw:
        ax.plot(r["t_mid"], r["tau1"], "o", color=col.get(r["phase"], "k"))
    ax.axvspan(res["flash_onset"], res["flash_end"], color="orange", alpha=0.2)
    ax.set_ylim(0, max(2.5, 1.2 * max(r["tau1"] for r in pw))); ax.set_xlabel("time since AMBIT sample 0 [s]")
    ax.set_ylabel("1-exp lifetime per window [ns]"); ax.set_title("lifetime per window (blue dark, red pulse, green after)")
    ax2 = ax.twinx(); ax2.plot([r["t_mid"] for r in pw], [r["laser_per_s"] for r in pw], "x", color="tab:purple", ms=4)
    ax2.set_ylabel("laser-correlated counts / s", color="tab:purple")
    ax = axs[1, 0]
    for f, name, c in ((Fo, "Fo (dark)", "tab:blue"), (Fm, "Fm (last %.1f s of the pulse)" % res["fm_seconds"], "tab:red")):
        d = np.clip(f["h"] - f["floor"], 1e-3, None)
        ax.plot(f["t"], d / d.max(), lw=0.9, color=c, label=f"{name}: 1-exp {f['tau1']:.2f} ns, amp-w {f['tau_amp2']:.2f} ns")
        ax.plot(f["fit2"]["t_fit"], (f["fit2"]["model"] - f["floor"]) / d.max(), lw=1.0, color="k", ls="--")
    ax.set_yscale("log"); ax.set_ylim(2e-3, 1.5); ax.set_xlim(max(0, Fo["peak_ns"] - 1.5), 12.5)
    ax.set_xlabel("time after SYNC [ns]"); ax.set_ylabel("normalised, floor subtracted"); ax.legend(fontsize=8)
    ax.set_title("decays with 2-exp tail fits (dashed)")
    ax = axs[1, 1]
    names = ["AMBIT ratio", "Felix intensity", "lifetime 1-exp", "lifetime amp-w"]
    vals = [a["FvFm"], i["FvFm"], L["FvFm_tau1"], L["FvFm_tau_amp"]]
    bars = ax.bar(names, vals, color=["tab:blue", "tab:purple", "tab:orange", "tab:red"])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.01, f"{v:.3f}", ha="center")
    ax.set_ylim(0, 1); ax.set_ylabel("Fv/Fm"); ax.set_title("Fv/Fm by method")
    fig.suptitle(Path(res["run"]).name)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=90)
    if show:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
def run_protocol(port: str | None, window_s: float, out_dir: Path, protocol=None, log=print) -> Path:
    protocol = protocol or FVFM_PROTOCOL
    ambit = ag.find_ambit(port, log=log)
    felix = ag.FelixLink(log=log).open()
    try:
        ag.run_gated(ambit, felix, protocol, window_s=window_s, out_dir=out_dir, log=log)
    finally:
        felix.close()
        ambit.close()
    return out_dir


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", help="analyse this saved run instead of measuring")
    ap.add_argument("--port", help="AMBIT serial port (default: scan)")
    ap.add_argument("--window", type=float, default=0.1, help="Felix window length in s (default 0.1)")
    ap.add_argument("--fm-last", type=float, default=0.2, help="Fm = last this many seconds of the pulse")
    ap.add_argument("--fo-last", type=float, default=1.0, help="Fo = last this many seconds of dark before the pulse")
    ap.add_argument("--out", help="output folder for a new run (default runs/<timestamp>_fvfm)")
    ap.add_argument("--show", action="store_true", help="open the figure window")
    args = ap.parse_args(argv)
    if args.run:
        run_dir = Path(args.run)
    else:
        run_dir = Path(args.out) if args.out else HERE / "runs" / (time.strftime("%Y%m%d_%H%M%S") + "_fvfm")
        run_protocol(args.port, args.window, run_dir)
    run = load_run(run_dir)
    res = fvfm(run, fm_last_s=args.fm_last, fo_last_s=args.fo_last)
    print("saved", save(res, run_dir))
    plot(res, run, save_to=run_dir / "fvfm.png", show=args.show)
    print("figure", run_dir / "fvfm.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
