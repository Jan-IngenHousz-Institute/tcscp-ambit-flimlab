"""Phi2 (PSII quantum yield in the light) from the AMBIT and from the TCSPC lifetime.

Each repeat is the AMBIT Phi2 protocol: 5 s at actinic X (10 Hz sampling), a 0.1 s
saturating pulse at 255 (100 Hz), 0.5 s at X (100 Hz). The repeats run back to back with
the LED kept on between commands (arrun persist flag), so the leaf stays in the light.
By default the actinic byte steps 10, 20, ... 100 over ten repeats, which spreads Phi2
over a range and lets a slope be told from an offset; pass --actinic 50 for ten repeats
at one level.

The Felix records 20 ms windows through each repeat. The pulse windows are found from
the window totals (the SPAD rate jumps during the pulse), so no timing alignment is
needed. Fs is the last second before the pulse, Fm' the windows fully inside the pulse.
Every histogram is fitted after subtracting the flat floor of LED-excited fluorescence.

Phi2 estimates compared:
  AMBIT        (Fm' - Fs) / Fm' from s_630 / r_630
  tau global   1 - tau_s / tau_m, amplitude-weighted lifetimes from a 2-exponential fit
               with the two lifetimes shared between Fs and Fm' (only amplitudes differ)
  tau amp      same but with independent 2-exponential fits
  tau 1-exp    1 - tau_s / tau_m from single-exponential fits
  intensity    1 - Is / Im from laser-correlated counts per second (no dead-time correction)

Needs the ND filter in front of the SPAD (pulse below about 3 Mcps). Use:
    python phi2_lifetime.py                         # ramp 10..100, ten repeats
    python phi2_lifetime.py --actinic 50 --repeats 10
    python phi2_lifetime.py --run runs/<dir>        # re-analyse a saved series
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
import ambit_gate as ag            # noqa: E402
import fvfm_lifetime as fl         # noqa: E402

WINDOW_S = 0.02
PRE_S = 0.3
POST_S = 0.4
FS_LAST_S = 1.0
REBIN = fl.REBIN


def phi2_protocol(x: int) -> list[dict]:
    return [
        dict(type=2, far_red=0, pulses=50, freq=10, actinic=int(x), subsampling=0),    # 5 s at X
        dict(type=2, far_red=0, pulses=10, freq=100, actinic=255, subsampling=0),      # 0.1 s saturating pulse
        dict(type=2, far_red=0, pulses=50, freq=100, actinic=int(x), subsampling=0),   # 0.5 s at X
    ]


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------
def measure(port: str | None, actinic: list[int], out_dir: Path, window_s=WINDOW_S, log=print) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ambit = ag.find_ambit(port, log=log)
    felix = ag.FelixLink(log=log).open()
    felix.set_thresholds(ag.THRESHOLD_V)
    rates = felix.count_rates(0.5)
    log("Felix input rates: " + ", ".join(f"ch{c}={r:,.0f}/s" for c, r in rates.items()))
    felix.configure_histogram(ag.SYNC_CH, ag.SPAD_CH, 0.0, 1.0 / ag.LASER_REP_RATE_HZ)
    led_on = False
    try:
        for i, x in enumerate(actinic):
            protocol = phi2_protocol(x)
            last = i == len(actinic) - 1
            cmd = ag.build_arrun(protocol, persist=0 if last else 1)
            n_win = int(np.ceil((PRE_S + ag.protocol_duration(protocol) + POST_S) / window_s))
            actual = felix.start_windows(window_s, n_win)
            time.sleep(PRE_S)
            log(f"repeat {i + 1}/{len(actinic)}: actinic {x}, {n_win} windows of {actual*1e3:.1f} ms")
            led_on = not last
            data, t_send, t_done = ag.run_arrun(ambit, cmd, protocol, log=log)
            felix.wait_done(timeout_s=POST_S + 5.0)
            wins = list(felix.windows)
            bins = np.vstack([w.bins for w in wins]).astype(np.float64)
            flags = np.array([w.lost or w.overflow for w in wins])
            s = data.get("s_630", np.zeros(0)).astype(float)
            r = data.get("r_630", np.zeros(0)).astype(float)
            np.savez_compressed(out_dir / f"repeat_{i:02d}.npz", bins47=fl._rebin_matrix(bins, REBIN).astype(np.uint32),
                                totals=bins.sum(1), flags=flags, t_soa=felix.t_soa, t_send=t_send, t_done=t_done,
                                window_s=felix.window_s, bin_s=felix.bin_s, s_630=s, r_630=r, actinic=x,
                                protocol=json.dumps(protocol), cmd=cmd)
            # full-resolution sums of the Fs and Fm' windows, chosen from the totals
            sel = classify(bins.sum(1), flags, felix.window_s)
            np.savez_compressed(out_dir / f"repeat_{i:02d}_decays.npz", fs=bins[sel["fs"]].sum(0),
                                fm=bins[sel["fm"]].sum(0), n_fs=len(sel["fs"]), n_fm=len(sel["fm"]),
                                bin_s=felix.bin_s, window_s=felix.window_s)
            log(f"   SPAD {bins.sum() / (len(wins) * felix.window_s) / 1e3:.0f} kcps mean, pulse windows {sel['fm']}, "
                f"{len(sel['fs'])} Fs windows, {int(flags.sum())} flagged")
    finally:
        if led_on:
            try:   # a one-sample dark command with persist 0 switches the LED off
                ag.run_arrun(ambit, ag.build_arrun([dict(type=2, far_red=0, pulses=1, freq=10, actinic=0, subsampling=0)]),
                             [dict(type=2, far_red=0, pulses=1, freq=10, actinic=0, subsampling=0)], log=log)
            except Exception as exc:
                log(f"could not switch the LED off: {exc}")
        felix.close()
        ambit.close()
    (out_dir / "series.json").write_text(json.dumps(dict(actinic=actinic, window_s=window_s, pre_s=PRE_S, post_s=POST_S,
                                                           rates_before=rates, felix=felix.description, ambit=ambit.name,
                                                           ambit_fw=ambit.fw), indent=1))
    return out_dir


# ---------------------------------------------------------------------------
# window selection from the totals
# ---------------------------------------------------------------------------
def classify(totals: np.ndarray, flags: np.ndarray, window_s: float, fs_last_s: float = FS_LAST_S) -> dict:
    """Pulse windows: within 10 % of the maximum total. Fs: the second before the pulse."""
    tot = np.asarray(totals, float)
    n = tot.size
    i_max = int(np.argmax(tot))
    before = tot[max(0, i_max - int(1.5 / window_s)): max(1, i_max - 8)]
    base = float(np.median(before)) if before.size else float(np.median(tot))
    level = tot[i_max]
    full = [i for i in range(max(0, i_max - 10), min(n, i_max + 11))
            if tot[i] >= base + 0.9 * (level - base) and not flags[i]]
    first = min(full)
    n_fs = int(round(fs_last_s / window_s))
    fs = [i for i in range(max(0, first - 2 - n_fs), first - 2)
          if not flags[i] and tot[i] < base + 0.3 * (level - base)]
    return dict(fm=full, fs=fs, base=base, level=level, i_max=i_max)


# ---------------------------------------------------------------------------
# fits
# ---------------------------------------------------------------------------
def _prep(h_full: np.ndarray, bin_s: float):
    h = fl._rebin(np.asarray(h_full, float), REBIN)
    bw = bin_s * REBIN * 1e9
    t = (np.arange(h.size) + 0.5) * bw
    return h, t, bw


def global_fit(h_s_full: np.ndarray, h_m_full: np.ndarray, bin_s: float) -> dict:
    """2-exponential tail fit of Fs and Fm' with shared lifetimes and free amplitudes."""
    h_s, t, bw = _prep(h_s_full, bin_s)
    h_m, _, _ = _prep(h_m_full, bin_s)
    peak = int(np.argmax(h_s))
    pre = (t > 0.3) & (t < t[peak] - 0.6)
    floors = [float(h_s[pre].mean()), float(h_m[pre].mean())]
    start = min(peak + max(3, int(round(fl.FIT_START_NS / bw))), h_s.size - 6)
    tt = t[start:] - t[start]
    cs = [h_s[start:], h_m[start:]]
    ws = [1.0 / np.sqrt(np.maximum(c, 1.0)) for c in cs]

    def model(p, j):
        a1, a2 = p[2 * j], p[2 * j + 1]
        return floors[j] + a1 * np.exp(-tt / p[4]) + a2 * np.exp(-tt / p[5])

    def resid(p):
        return np.concatenate([(model(p, j) - cs[j]) * ws[j] for j in range(2)])

    p0 = np.array([max(cs[0][0] - floors[0], 1) / 2] * 2 + [max(cs[1][0] - floors[1], 1) / 2] * 2 + [0.3, 1.5])
    lo = np.array([0, 0, 0, 0, 0.1, 0.1])      # a tail fit from 0.25 ns cannot resolve faster components
    hi = np.array([np.inf] * 4 + [fl.TAU_BOUNDS[1]] * 2)
    sol = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac", max_nfev=600)
    taus = sol.x[4:6]
    order = np.argsort(taus)
    taus = taus[order]
    out = dict(taus=taus, chi2r=float(np.sum(sol.fun ** 2) / max(2 * tt.size - 6, 1)), ok=bool(sol.success), floors=floors)
    for j, name in enumerate(("s", "m")):
        a = sol.x[2 * j: 2 * j + 2][order] * np.exp((t[start] - t[peak]) / taus)   # referred to the peak
        tau_amp = float(np.sum(a * taus) / np.sum(a)) if a.sum() > 0 else float("nan")
        out[f"amps_{name}"] = a
        out[f"fractions_{name}"] = a / a.sum() if a.sum() > 0 else a
        out[f"tau_amp_{name}"] = tau_amp
        out[f"tau_int_{name}"] = float(np.sum(a * taus ** 2) / np.sum(a * taus)) if (a * taus).sum() > 0 else float("nan")
    return out


# ---------------------------------------------------------------------------
# analysis of a saved series
# ---------------------------------------------------------------------------
def analyse(run_dir: Path, fs_last_s: float = FS_LAST_S, log=print) -> dict:
    run_dir = Path(run_dir)
    series = json.loads((run_dir / "series.json").read_text())
    rows = []
    for f in sorted(run_dir.glob("repeat_??.npz")):
        z = np.load(f, allow_pickle=False)
        d = np.load(f.with_name(f.stem + "_decays.npz"))
        x = int(z["actinic"]); W = float(z["window_s"]); bin_s = float(z["bin_s"])
        protocol = json.loads(str(z["protocol"]))
        s, r = z["s_630"], z["r_630"]
        ratio = np.divide(s, r, out=np.full(s.size, np.nan), where=r > 0)
        n0, n1 = protocol[0]["pulses"], protocol[1]["pulses"]
        fs_n = int(round(fs_last_s * protocol[0]["freq"]))
        Fs_A = float(np.nanmean(ratio[n0 - fs_n:n0]))
        Fm_A = float(np.nanmax(ratio[n0:n0 + n1]))
        sel = classify(z["totals"], z["flags"], W, fs_last_s)
        hs, hm = d["fs"].astype(float), d["fm"].astype(float)
        n_fs, n_fm = int(d["n_fs"]), int(d["n_fm"])
        g = global_fit(hs, hm, bin_s)
        fs_fit, fm_fit = fl.fit_window(hs, bin_s), fl.fit_window(hm, bin_s)
        row = dict(
            repeat=len(rows), actinic=x, Fs_ambit=Fs_A, Fm_ambit=Fm_A, phi2_ambit=(Fm_A - Fs_A) / Fm_A,
            spad_fs_kcps=float(hs.sum() / (n_fs * W) / 1e3), spad_fm_kcps=float(hm.sum() / (n_fm * W) / 1e3),
            n_fs=n_fs, n_fm=n_fm, fm_windows=sel["fm"],
            taus_global=[float(v) for v in g["taus"]], chi2_global=g["chi2r"],
            tau_amp_s=g["tau_amp_s"], tau_amp_m=g["tau_amp_m"], fractions_s=[float(v) for v in g["fractions_s"]],
            fractions_m=[float(v) for v in g["fractions_m"]], phi2_tau_global=1 - g["tau_amp_s"] / g["tau_amp_m"],
            tau_amp2_s=fs_fit["tau_amp2"], tau_amp2_m=fm_fit["tau_amp2"], phi2_tau_amp=1 - fs_fit["tau_amp2"] / fm_fit["tau_amp2"],
            tau1_s=fs_fit["tau1"], tau1_m=fm_fit["tau1"], phi2_tau1=1 - fs_fit["tau1"] / fm_fit["tau1"],
            laser_s=fs_fit["laser_counts"] / (n_fs * W), laser_m=fm_fit["laser_counts"] / (n_fm * W),
            laser_counts_fs=fs_fit["laser_counts"], laser_counts_fm=fm_fit["laser_counts"],
        )
        row["phi2_intensity"] = 1 - row["laser_s"] / row["laser_m"]
        row["_ratio"] = ratio; row["_decays"] = (hs, hm, g, fs_fit, fm_fit, bin_s)
        rows.append(row)
    A = np.array([r["phi2_ambit"] for r in rows])
    res = dict(run=str(run_dir), series=series, rows=rows, primary=PRIMARY)
    # pooled: all Fs decays summed, all Fm' decays summed (meaningful when X is the same)
    hs_all = sum(r["_decays"][0] for r in rows); hm_all = sum(r["_decays"][1] for r in rows)
    bin_s = rows[0]["_decays"][5]
    fs_p, fm_p = fl.fit_window(hs_all, bin_s), fl.fit_window(hm_all, bin_s)
    g_p = global_fit(hs_all, hm_all, bin_s)
    Fs_all = float(np.mean([r["Fs_ambit"] for r in rows])); Fm_all = float(np.mean([r["Fm_ambit"] for r in rows]))
    res["pooled"] = dict(
        phi2_ambit=(Fm_all - Fs_all) / Fm_all, phi2_ambit_mean=float(np.mean(A)),
        phi2_tau_amp=1 - fs_p["tau_amp2"] / fm_p["tau_amp2"], tau_amp2_s=fs_p["tau_amp2"], tau_amp2_m=fm_p["tau_amp2"],
        taus2_s=[float(v) for v in fs_p["taus2"]], taus2_m=[float(v) for v in fm_p["taus2"]],
        fractions2_s=[float(v) for v in fs_p["fractions2"]], fractions2_m=[float(v) for v in fm_p["fractions2"]],
        chi2_s=fs_p["chi2_2"], chi2_m=fm_p["chi2_2"],
        phi2_tau1=1 - fs_p["tau1"] / fm_p["tau1"], tau1_s=fs_p["tau1"], tau1_m=fm_p["tau1"],
        phi2_tau_global=1 - g_p["tau_amp_s"] / g_p["tau_amp_m"], taus_global=[float(v) for v in g_p["taus"]],
        laser_counts_fs=fs_p["laser_counts"], laser_counts_fm=fm_p["laser_counts"],
        phi2_intensity=1 - (fs_p["laser_counts"] / sum(r["n_fs"] for r in rows)) / (fm_p["laser_counts"] / sum(r["n_fm"] for r in rows)),
    )
    for key in ("phi2_tau_amp", "phi2_tau_global", "phi2_tau1"):
        L = np.array([r[key] for r in rows])
        fit = np.polyfit(A, L, 1) if len(rows) >= 3 and np.ptp(A) > 0.02 else [np.nan, np.nan]
        res[key] = dict(offset_mean=float(np.mean(L - A)),
                        offset_sd=float(np.std(L - A, ddof=1)) if len(rows) > 1 else float("nan"),
                        slope=float(fit[0]), intercept=float(fit[1]))
    log(report(res))
    return res


PRIMARY = "phi2_tau_amp"      # what the user asked for: independent 2-exp, amplitude-weighted
LABEL = {"phi2_tau_amp": "lifetime, 2-exp amp-weighted", "phi2_tau_global": "lifetime, shared-tau 2-exp",
         "phi2_tau1": "lifetime, 1-exp"}


def report(res: dict) -> str:
    L = [f"{'rep':>3} {'X':>3} {'SPAD Fs/Fm kcps':>15} {'Fs_A':>7} {'Fm_A':>7} {'Phi2 AMBIT':>10} {'tau_s':>6} {'tau_m':>6} "
         f"{'Phi2 tau glob':>13} {'Phi2 tau amp':>12} {'Phi2 1-exp':>10} {'Phi2 int':>8} {'taus glob':>14} {'cts Fm':>6}"]
    for r in res["rows"]:
        L.append(f"{r['repeat']:3d} {r['actinic']:3d} {r['spad_fs_kcps']:7.0f}/{r['spad_fm_kcps']:<7.0f} {r['Fs_ambit']:7.4f} {r['Fm_ambit']:7.4f} "
                 f"{r['phi2_ambit']:10.3f} {r['tau_amp_s']:6.3f} {r['tau_amp_m']:6.3f} {r['phi2_tau_global']:13.3f} "
                 f"{r['phi2_tau_amp']:12.3f} {r['phi2_tau1']:10.3f} {r['phi2_intensity']:8.3f} "
                 f"{str(np.round(r['taus_global'], 2)):>14} {r['laser_counts_fm']:6.0f}")
    for key in ("phi2_tau_amp", "phi2_tau_global", "phi2_tau1"):
        s = res[key]
        L.append(f"{LABEL[key]:>30} minus AMBIT: mean offset {s['offset_mean']:+.3f}, sd {s['offset_sd']:.3f}; "
                 f"linear fit = {s['slope']:.3f} x Phi2_AMBIT {s['intercept']:+.3f}")
    p = res["pooled"]
    L.append(f"pooled over all repeats ({p['laser_counts_fs']:.0f} laser counts in Fs, {p['laser_counts_fm']:.0f} in Fm'): "
             f"AMBIT {p['phi2_ambit']:.3f} (mean of repeats {p['phi2_ambit_mean']:.3f}) | "
             f"2-exp amp-weighted {p['phi2_tau_amp']:.3f} (tau_s {p['tau_amp2_s']:.3f}: {np.round(p['taus2_s'], 2)} {np.round(p['fractions2_s'], 2)}, "
             f"tau_m {p['tau_amp2_m']:.3f}: {np.round(p['taus2_m'], 2)} {np.round(p['fractions2_m'], 2)}) | "
             f"1-exp {p['phi2_tau1']:.3f} ({p['tau1_s']:.3f} -> {p['tau1_m']:.3f} ns) | shared-tau {p['phi2_tau_global']:.3f} | "
             f"intensity {p['phi2_intensity']:.3f}")
    return "\n".join(L)


def save(res: dict, run_dir: Path) -> Path:
    rows = [{k: v for k, v in r.items() if not k.startswith("_")} for r in res["rows"]]
    out = {k: v for k, v in res.items() if k != "rows"}
    out["rows"] = rows
    p = Path(run_dir) / "phi2.json"
    p.write_text(json.dumps(out, indent=1, default=float))
    with open(Path(run_dir) / "phi2.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)
    return p


def fvfm_point(fvfm_dir: Path | None) -> dict | None:
    """Dark-adapted Fv/Fm from a saved fvfm run, with the same global fit, for the scatter plot."""
    if fvfm_dir is None or not (Path(fvfm_dir) / "fvfm.json").exists():
        return None
    j = json.loads((Path(fvfm_dir) / "fvfm.json").read_text())
    run = fl.load_run(fvfm_dir)
    hs = run["bins"][j["fo_windows"]].sum(0); hm = run["bins"][j["fm_windows"]].sum(0)
    g = global_fit(hs, hm, run["bin_s"])
    fo, fm = fl.fit_window(hs, run["bin_s"]), fl.fit_window(hm, run["bin_s"])
    return dict(phi2_ambit=j["ambit"]["FvFm"], phi2_tau_global=1 - g["tau_amp_s"] / g["tau_amp_m"],
                phi2_tau_amp=1 - fo["tau_amp2"] / fm["tau_amp2"], phi2_tau1=1 - fo["tau1"] / fm["tau1"],
                tau_amp_s=fo["tau_amp2"], tau_amp_m=fm["tau_amp2"], taus=g["taus"].tolist(), run=str(fvfm_dir))


def plot(res: dict, save_to: Path | None = None, show: bool = False, fvfm: dict | None = None):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = res["rows"]
    x = np.array([r["actinic"] for r in rows]); idx = np.arange(len(rows))
    A = np.array([r["phi2_ambit"] for r in rows]); G = np.array([r["phi2_tau_global"] for r in rows])
    T2 = np.array([r["phi2_tau_amp"] for r in rows]); T1 = np.array([r["phi2_tau1"] for r in rows])
    fig, axs = plt.subplots(2, 2, figsize=(14, 9))
    ax = axs[0, 0]
    xl = x if np.ptp(x) > 0 else idx
    ax.plot(xl, A, "o-", color="tab:blue", label="AMBIT (Fm'-Fs)/Fm'")
    ax.plot(xl, G, "s-", color="tab:red", label="lifetime, global 2-exp amp-weighted")
    ax.plot(xl, T2, "^--", color="tab:orange", ms=5, label="lifetime, independent 2-exp amp-weighted")
    ax.plot(xl, T1, "v:", color="grey", ms=5, label="lifetime, 1-exp")
    ax.set_xlabel("actinic byte X" if np.ptp(x) > 0 else "repeat"); ax.set_ylabel("Phi2"); ax.set_ylim(0, 1)
    ax.legend(fontsize=8); ax.grid(lw=0.3); ax.set_title("Phi2 per repeat")
    ax = axs[0, 1]
    P = np.array([r[PRIMARY] for r in rows]); st = res[PRIMARY]
    ax.plot(A, P, "s", color="tab:red", label=LABEL[PRIMARY])
    ax.plot(A, G, "o", mfc="none", color="tab:orange", label=LABEL["phi2_tau_global"])
    for a, p, xx in zip(A, P, x):
        ax.annotate(str(xx), (a, p), fontsize=7, xytext=(3, 3), textcoords="offset points")
    if fvfm:
        ax.plot(fvfm["phi2_ambit"], fvfm[PRIMARY], "D", color="k", label="dark Fv/Fm, earlier run (same estimator)")
    lim = [0, 1]; ax.plot(lim, lim, "k--", lw=0.7, label="1:1")
    if np.isfinite(st["slope"]):
        xx = np.linspace(min(A.min(), 0.2), max(A.max(), 0.9), 10)
        ax.plot(xx, st["slope"] * xx + st["intercept"], "-", color="tab:red", lw=0.8,
                label=f"fit: {st['slope']:.2f} x {st['intercept']:+.3f}")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel("Phi2 AMBIT"); ax.set_ylabel("Phi2 from lifetime"); ax.legend(fontsize=8)
    ax.grid(lw=0.3); ax.set_title(f"{LABEL[PRIMARY]}: offset {st['offset_mean']:+.3f} +- {st['offset_sd']:.3f}")
    ax = axs[1, 0]
    ax.plot(xl, [r["tau_amp2_s"] for r in rows], "s-", color="tab:green", label="tau_s (Fs, 2-exp amp-weighted)")
    ax.plot(xl, [r["tau_amp2_m"] for r in rows], "s-", color="tab:red", label="tau_m (Fm', 2-exp amp-weighted)")
    ax.plot(xl, [r["tau_amp_s"] for r in rows], "o--", mfc="none", color="tab:green", label="tau_s (shared-tau fit)")
    ax.plot(xl, [r["tau_amp_m"] for r in rows], "o--", mfc="none", color="tab:red", label="tau_m (shared-tau fit)")
    ax.set_xlabel(axs[0, 0].get_xlabel()); ax.set_ylabel("lifetime [ns]"); ax.legend(fontsize=8); ax.grid(lw=0.3)
    ax.set_title("lifetimes")
    ax = axs[1, 1]
    t_off = 0.0
    for r in rows:
        ratio = r["_ratio"]; t_amb, _ = ag.time_axis(phi2_protocol(r["actinic"]))
        n = min(len(ratio), len(t_amb))
        ax.plot(t_off + t_amb[:n], ratio[:n], "-", lw=0.8, color="tab:blue")
        ax.text(t_off + 2.0, np.nanmax(ratio) * 1.02, f"X={r['actinic']}", fontsize=7)
        t_off += t_amb[n - 1] + 1.0
    ax.set_xlabel("time [s] (repeats concatenated, 1 s gap)"); ax.set_ylabel("AMBIT s_630 / r_630"); ax.grid(lw=0.3)
    ax.set_title("AMBIT trace")
    fig.suptitle(Path(res["run"]).name); fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=90)
    if show:
        plt.show()
    return fig


def _latest_fvfm() -> Path | None:
    runs = sorted((HERE / "runs").glob("*_fvfm"))
    return runs[-1] if runs else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", help="analyse this saved series instead of measuring")
    ap.add_argument("--port", help="AMBIT serial port (default: scan)")
    ap.add_argument("--actinic", default="10,20,30,40,50,60,70,80,90,100",
                    help="actinic byte(s): one value, or a comma list used in order (default ramp 10..100)")
    ap.add_argument("--repeats", type=int, help="number of repeats (default: length of the list, or 10)")
    ap.add_argument("--window", type=float, default=WINDOW_S)
    ap.add_argument("--fs-last", type=float, default=FS_LAST_S, help="Fs = last this many seconds before the pulse")
    ap.add_argument("--fvfm", help="fvfm run folder for the dark point in the scatter (default: latest *_fvfm)")
    ap.add_argument("--out")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args(argv)
    if args.run:
        run_dir = Path(args.run)
    else:
        vals = [int(v) for v in args.actinic.split(",") if v.strip()]
        n = args.repeats or (len(vals) if len(vals) > 1 else 10)
        actinic = (vals * n)[:n] if len(vals) == 1 else vals[:n]
        run_dir = Path(args.out) if args.out else HERE / "runs" / (time.strftime("%Y%m%d_%H%M%S") + "_phi2")
        measure(args.port, actinic, run_dir, window_s=args.window)
    res = analyse(run_dir, fs_last_s=args.fs_last)
    print("saved", save(res, run_dir))
    fv = fvfm_point(Path(args.fvfm) if args.fvfm else _latest_fvfm())
    if fv:
        print(f"dark Fv/Fm point from {fv['run']}: AMBIT {fv['phi2_ambit']:.3f}, lifetime 2-exp amp-weighted "
              f"{fv['phi2_tau_amp']:.3f} (tau_o {fv['tau_amp_s']:.3f}, tau_m {fv['tau_amp_m']:.3f}), "
              f"shared-tau {fv['phi2_tau_global']:.3f}, 1-exp {fv['phi2_tau1']:.3f}")
    plot(res, save_to=run_dir / "phi2.png", show=args.show, fvfm=fv)
    print("figure", run_dir / "phi2.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
