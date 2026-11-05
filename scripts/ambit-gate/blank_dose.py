"""Blank and laser-dose tests of the AMBIT + Felix set-up (notebook ``3_blank_laser_dose.ipynb``).

The Felix runs in histogram mode here (micro-times referenced to the laser SYNC by the hardware):
no chopping, no gating and no micro-time reconstruction, so the decays are the instrument's own.

Blank
    What both instruments read without a leaf. Felix: laser-correlated photons and flat floor
    with the leaf removed, laser on and laser blocked. AMBIT: the same chopped protocols as the
    leaf runs; its s/r per window is the offset to subtract from the leaf runs.

Laser dose
    A leaf dark-adapted with the laser blocked. The AMBIT samples its yield while the laser is
    unblocked by hand; the Felix records the transient from the first laser photon on. If the
    laser closes PSII centres in its spot, the AMBIT yield rises after laser-on and the Felix
    lifetime of the first 100 ms is shorter than the steady one.
"""
from __future__ import annotations

import csv
import json
import sys
import threading
import time
from pathlib import Path

import numpy as np

import ambit_chop as ac
import ambit_gate as ag
import fvfm_lifetime as fl

HERE = Path(__file__).resolve().parent
IRF_DIR = HERE.parent.parent / "irf-library"
DEFAULT_IRF_NAME = "261104-IRF-blackplastic-50p-10minwarmup_Histogram 0.csv"   # shipped with the phi2 GUI
REBIN = fl.REBIN              # 20 x 2.34 ps = 47 ps bins for floors, counts and fits
PEAK_AT_NS = 2.5              # every decay is rotated so its laser peak sits here
RISE_NS = 0.5                 # laser-correlated counts from this long before the peak: a scatter
                              # (IRF-shaped) peak sits ~0.2 ns before a fluorescence peak
DEAD_TIME_S = ac.DEAD_TIME_S
SPAD_MAX_CPS = 3.5e6          # SPAD-050 rated maximum count rate (3-3.5 Mcps, max 50 uW incident)
MIN_FAST_NS = fl.MIN_RESOLVED_NS   # a tail fit from 0.25 ns after the peak cannot resolve faster components
# time since laser-on (s) of the kinetics bins: fine at the start, coarse later. The bins before 0
# hold the unblocking ramp (laser-on = the window at half of the steady level, see _laser_on)
KIN_EDGES_S = np.array([-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2,
                        0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 60.0])


# ---------------------------------------------------------------------------
# One decay: floor, laser-correlated counts, mean arrival, tail fits
# ---------------------------------------------------------------------------
def _rotate(h: np.ndarray, bin_s: float, peak_bin: int | None = None):
    """Rebinned decay rotated so the laser peak sits at PEAK_AT_NS; returns (h, bin_ns, peak_bin).

    ``peak_bin`` (rebinned index) forces the rotation, e.g. a blank on the leaf's peak.
    """
    hr = fl._rebin(np.asarray(h, float), REBIN)
    bw = bin_s * REBIN * 1e9
    if peak_bin is None:
        sm = np.convolve(np.r_[hr[-2:], hr, hr[:2]], np.ones(5) / 5, mode="valid")
        peak_bin = int(np.argmax(sm))
    return np.roll(hr, int(round(PEAK_AT_NS / bw)) - peak_bin), bw, int(peak_bin)


def _fwhm_ps(h: np.ndarray, bin_s: float, shift: int) -> float:
    y = np.convolve(np.roll(np.asarray(h, float), shift), np.ones(9) / 9, mode="same")
    pre = int((PEAK_AT_NS - 0.6) * 1e-9 / bin_s)
    y = y - np.median(y[int(0.3e-9 / bin_s):pre])
    pk = int(np.argmax(y))
    half = y[pk] / 2
    lo, hi = pk, pk
    while lo > 0 and y[lo] > half:
        lo -= 1
    while hi < y.size - 1 and y[hi] > half:
        hi += 1
    left = lo + (half - y[lo]) / max(y[lo + 1] - y[lo], 1e-12)
    right = hi - 1 + (y[hi - 1] - half) / max(y[hi - 1] - y[hi], 1e-12)
    return float((right - left) * bin_s * 1e12)


def decay_summary(h: np.ndarray, bin_s: float, seconds: float, peak_bin: int | None = None,
                  dead_time_s: float = DEAD_TIME_S) -> dict:
    """Rates (dead-time corrected, non-paralyzable), peak, FWHM, mean arrival and tail fits.

    Floor = mean of the rotated decay from 0.3 ns to 0.6 ns before the peak, as in
    ``ambit_chop``; laser-correlated = everything above it from RISE_NS before the peak on.
    ``significant``: laser-correlated counts above 5 sigma; only then are the mean arrival,
    FWHM and fits given.
    """
    if seconds <= 0:
        raise ValueError("decay_summary: no windows selected")
    h = np.asarray(h, float)
    hr, bw, pk = _rotate(h, bin_s, peak_bin)
    t = (np.arange(hr.size) + 0.5) * bw
    p0 = int(round(PEAK_AT_NS / bw))
    pre = (t > 0.3) & (t < PEAK_AT_NS - 0.6)
    floor = float(hr[pre].mean())
    floor_err = float(hr[pre].std(ddof=1) / np.sqrt(pre.sum()))
    rise = p0 - int(round(RISE_NS / bw))
    lc = float((hr[rise:] - floor).sum())
    lc_err = float(np.sqrt(max(hr[rise:].sum(), 1.0) + (floor_err * (hr.size - rise)) ** 2))
    m = float(h.sum()) / seconds if seconds > 0 else float("nan")
    corr = 1.0 / (1.0 - m * dead_time_s) if m * dead_time_s < 0.9 else float("nan")
    sel = (t >= t[p0]) & (t < t[p0] + 6.0)
    w = hr[sel] - floor            # not clipped at 0: clipping the noise in the tail biases t_mean up
    out = dict(seconds=float(seconds), counts=float(h.sum()), total_cps=m, dead_time_factor=corr,
               laser_cps=lc / seconds * corr, laser_cps_err=lc_err / seconds * corr,
               floor_cps=floor * hr.size / seconds * corr, floor_per_bin=floor,
               laser_fraction=lc / max(hr.sum(), 1.0), significant=bool(lc > 5 * lc_err), peak_bin=pk,
               peak_ns=float((pk + 0.5) * bw), t_mean_ns=float("nan"), t_mean_err_ns=float("nan"),
               fwhm_ps=float("nan"), tau1_ns=float("nan"), tau_int2_ns=float("nan"), tau_avg2_ns=float("nan"),
               amps2=[], taus2_ns=[],
               fractions2=[], chi2_1=float("nan"), chi2_2=float("nan"), h=hr, bin_ns=bw, fit=None)
    if out["significant"]:
        if w.sum() > 0:
            dt = t[sel] - t[p0]
            tm = float((dt * w).sum() / w.sum())
            # Poisson noise of every bin plus the uncertainty of the subtracted floor
            var = ((dt - tm) ** 2 * hr[sel]).sum() + (floor_err * (dt - tm).sum()) ** 2
            out["t_mean_ns"], out["t_mean_err_ns"] = tm, float(np.sqrt(var) / w.sum())
        out["fwhm_ps"] = _fwhm_ps(h, bin_s, (p0 - pk) * REBIN)
        if lc > 300:
            f = fl.fit_window(hr, bin_s * REBIN, rebin=1)
            if f.get("ok"):
                out.update(tau1_ns=float(f["tau1"]), tau_int2_ns=float(f["tau_int2"]),
                           taus2_ns=[float(x) for x in f["taus2"]],
                           fractions2=[float(x) for x in f["fractions2"]],
                           chi2_1=float(f["chi2_1"]), chi2_2=float(f["chi2_2"]),
                           # the fitted curves, for plot_fits: same time axis as ``h``
                           fit=dict(t=f["fit1"]["t_fit"], model1=f["fit1"]["model"], model2=f["fit2"]["model"],
                                    floor=float(f["floor"]), degenerate=bool(f["fit2_degenerate"])))
                # 2-exp amplitudes as fitted: counts per (rebinned) bin at the laser peak, and the
                # amplitude-weighted lifetime (tau1 a1 + tau2 a2) / (a1 + a2); always the 2-exp
                # value (no fallback to the 1-exp lifetime)
                a2 = np.array([float(x) for x in f["fit2"]["amps"]])
                out["amps2"] = a2.tolist()
                out["tau_avg2_ns"] = float(np.sum(np.array(out["taus2_ns"]) * a2) / a2.sum()) if a2.sum() > 0 \
                    else float("nan")
    return out


def _public(s: dict) -> dict:
    """decay_summary without the arrays, for JSON and tables."""
    return {k: v for k, v in s.items() if k not in ("h", "fit")}


def _shape(s: dict):
    """(t_ns, floor-subtracted decay normalised to its peak) of a decay_summary."""
    t = (np.arange(s["h"].size) + 0.5) * s["bin_ns"]
    y = s["h"] - s["h"][(t > 0.3) & (t < PEAK_AT_NS - 0.6)].mean()
    return t, np.clip(y / y.max(), 1e-4, None)


def irf_dirs(extra=()) -> list[Path]:
    """Folders searched for IRF files, first wins: ``extra`` (e.g. the user's own folder), the
    folder shipped inside the frozen phi2 GUI, the repository's irf-library."""
    dirs = [Path(p) for p in extra]
    if getattr(sys, "frozen", False):
        dirs.append(Path(getattr(sys, "_MEIPASS", "")) / "irf-library")
    dirs.append(IRF_DIR)
    return [d for d in dirs if d.is_dir()]


def list_irfs(extra=()) -> list[Path]:
    """The available IRF files (*.csv), one per file name, the user's own folders first."""
    out, seen = [], set()
    for d in irf_dirs(extra):
        for p in sorted(d.glob("*.csv")):
            if p.name not in seen:
                seen.add(p.name)
                out.append(p)
    return out


def load_irf(path=None, extra=()) -> dict | None:
    """The measured IRF (TDC-Software CSV export: '# ...' header, 'seconds,counts' rows).

    ``path``: a file, or a file name looked up in ``irf_dirs(extra)``. None: ``DEFAULT_IRF_NAME``
    when available, else the first file found. Returns None when there is no IRF at all.
    """
    if path is None or not Path(path).is_file():
        found = list_irfs(extra)
        if not found:
            return None
        want = Path(path).name if path else DEFAULT_IRF_NAME
        path = next((p for p in found if p.name == want), found[0])
    t, c = np.loadtxt(path, delimiter=",", comments="#", unpack=True)
    return dict(bin_s=float(t[1] - t[0]), counts=c.astype(float), name=Path(path).name, path=str(path))


def _irf_curve(irf: dict, period_s: float = 1.0 / ag.LASER_REP_RATE_HZ):
    """(t_ns, floor-subtracted IRF normalised to its peak), one laser period rotated like the decays."""
    n = int(round(period_s / irf["bin_s"]))
    nz = np.flatnonzero(irf["counts"] > 0)         # the export pads outside the TDC range with 0
    c = irf["counts"][nz[0]:nz[0] + n]
    k = int(np.argmax(c))
    c = np.roll(c, int(round(PEAK_AT_NS * 1e-9 / irf["bin_s"])) - k)
    t = (np.arange(c.size) + 0.5) * irf["bin_s"] * 1e9
    c = c - np.median(c[(t > 0.3) & (t < PEAK_AT_NS - 0.6)])
    return t, np.clip(c / c.max(), 1e-4, None)


# ---------------------------------------------------------------------------
# Felix histogram recordings (no AMBIT)
# ---------------------------------------------------------------------------
def record_felix(felix: ag.FelixLink, out_dir, condition: str, seconds: float = 10.0,
                 window_s: float = 0.1, note: str = "", plot: bool = False, threshold_v: float = ag.THRESHOLD_V,
                 laser_hz: float = ag.LASER_REP_RATE_HZ, log=print) -> Path:
    """``seconds`` of contiguous hardware histograms; saves felix.npz + run.json in ``out_dir``.

    ``plot``: also draw the summed histogram with its FWHM (``plot_recording``), shown and
    saved as histogram_fwhm.png in ``out_dir``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    felix.set_thresholds(threshold_v)
    rates = felix.count_rates(0.5)
    felix.configure_histogram(ag.SYNC_CH, ag.SPAD_CH, 0.0, 1.0 / laser_hz)
    n = max(1, int(round(seconds / window_s)))
    w = felix.start_windows(window_s, n)
    ok = felix.wait_done(timeout_s=n * w + 10.0)
    wins = list(felix.windows)
    if not ok or len(wins) < n:
        log(f"WARNING: the Felix delivered {len(wins)} of {n} windows")
    if not wins:
        raise RuntimeError("no Felix windows")
    bins = np.vstack([x.bins for x in wins]).astype(np.uint32)
    np.savez_compressed(out_dir / "felix.npz", bins=bins, window_s=w, bin_s=felix.bin_s, left_s=felix.left_s)
    meta = dict(kind="felix", condition=condition, note=note, seconds=len(wins) * w, window_s=w,
                bin_s=felix.bin_s, laser_hz=laser_hz, rates={str(k): v for k, v in rates.items()},
                lost=int(sum(x.lost for x in wins)), felix=felix.description,
                saved=time.strftime("%Y-%m-%d %H:%M:%S"))
    (out_dir / "run.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    s = decay_summary(bins.sum(0), felix.bin_s, len(wins) * w)
    peak = (f"laser-correlated {s['laser_cps'] / 1e3:.2f} +- {s['laser_cps_err'] / 1e3:.2f} kcps"
            if s["significant"] else "no laser peak of its own")
    log(f"{condition}: {s['total_cps'] / 1e3:.1f} kcps total, {peak}, floor {s['floor_cps'] / 1e3:.1f} kcps"
        f" -> {out_dir}")
    if plot:
        plot_recording(out_dir, save=True, show=True)
    return out_dir


def _cross(x: np.ndarray, y: np.ndarray, i: int, level: float, step: int) -> float:
    """Walk from index ``i`` (the peak) by ``step`` until ``y`` drops below ``level``; linear
    interpolation between the two bins around the crossing."""
    j = i
    while 0 <= j + step < y.size and y[j + step] >= level:
        j += step
    k = j + step
    if not 0 <= k < y.size:
        return float("nan")
    return float(x[k] + (level - y[k]) * (x[j] - x[k]) / (y[j] - y[k]))


def peak_fwhm(h: np.ndarray, bin_s: float, left_s: float = 0.0, rebin_ps: float = 9.4) -> dict:
    """Peak, flat background and FWHM of one histogram, as in ``scripts/compute_fwhm.ipynb``.

    The crossings are read on a light rebin (~``rebin_ps``, mean so it stays in counts per raw
    bin) so the bin-to-bin TDC ripple does not move them. Unlike compute_fwhm the background is
    the mean of 4 ns to 0.6 ns *before* the peak (wrapped over the period): a leaf decay is not
    flat 3.5 ns after its peak.
    """
    c = np.asarray(h, float)
    bw = bin_s * 1e9
    t = left_s * 1e9 + (np.arange(c.size) + 0.5) * bw
    k = max(1, int(round(rebin_ps / (bw * 1e3)))) if rebin_ps > 0 else 1
    n = (c.size // k) * k
    tr, cr = t[:n].reshape(-1, k).mean(1), c[:n].reshape(-1, k).mean(1)
    ip = int(np.argmax(cr))
    tp, peak = float(tr[ip]), float(cr[ip])
    before = (tp - t) % (c.size * bw)
    flat = (before > 0.6) & (before < 4.0)
    bg = float(c[flat].mean())
    half = bg + 0.5 * (peak - bg)
    t_l, t_r = _cross(tr, cr, ip, half, -1), _cross(tr, cr, ip, half, +1)
    return dict(t=t, c=c, bw_ns=bw, tr=tr, cr=cr, k=k, peak_ns=tp, peak=peak, bg=bg, half=half,
                t_left=t_l, t_right=t_r, fwhm_ps=(t_r - t_l) * 1e3)


def plot_histogram(h: np.ndarray, bin_s: float, title: str, seconds: float | None = None, left_s: float = 0.0,
                   view_ns=(-0.7, 1.1), rebin_ps: float = 9.4, save_to=None, show: bool = True):
    """Raw histogram around the peak with the FWHM drawn on it, plus the full period on a log
    scale as an inset (the figure of ``scripts/compute_fwhm.ipynb``). Without a laser peak
    (laser blocked, background only) the whole period is drawn and no FWHM is given."""
    import matplotlib.pyplot as plt

    r = peak_fwhm(h, bin_s, left_s, rebin_ps)
    t, c, bw = r["t"], r["c"], r["bw_ns"]
    has_peak = decay_summary(h, bin_s, seconds or 1.0)["significant"]
    r["has_peak"] = has_peak
    if not has_peak:
        r["fwhm_ps"] = float("nan")
    tp, peak, half = r["peak_ns"], r["peak"], r["half"]
    lo, hi = (tp + view_ns[0], tp + view_ns[1]) if has_peak else (t[0] - bw / 2, t[-1] + bw / 2)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    m = (t > lo) & (t < hi)
    ax.step(t[m], c[m], where="mid", lw=0.6, color="0.45", label=f"raw data, {bw * 1e3:.2f} ps bins")
    if r["k"] > 1:
        mr = (r["tr"] > lo) & (r["tr"] < hi)
        ax.plot(r["tr"][mr], r["cr"][mr], lw=1.4, color="tab:blue",
                label=f"{r['k'] * bw * 1e3:.1f} ps rebin (used for the crossings)")
    ax.axhline(r["bg"], color="0.3", lw=0.8, ls=":", alpha=0.6, label=f"background {r['bg']:.1f} counts/bin")
    if has_peak:
        ax.axhline(peak, color="tab:blue", lw=0.8, ls=":", alpha=0.6)
        ax.axhline(half, color="tab:red", lw=0.8, ls="--", alpha=0.6)
        ax.vlines([r["t_left"], r["t_right"]], 0, half, color="tab:red", lw=0.8, ls="--", alpha=0.6)
        ax.annotate("", xy=(r["t_right"], half), xytext=(r["t_left"], half),
                    arrowprops=dict(arrowstyle="<->", color="tab:red", lw=2))
        ax.text((r["t_left"] + r["t_right"]) / 2, half, f"FWHM = {r['fwhm_ps']:.0f} ps",
                color="tab:red", ha="center", va="bottom", fontsize=14, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="tab:red", alpha=0.9))
        ax.text(tp + 0.02, peak, f"peak {peak:,.0f} counts/bin at {tp:.2f} ns", fontsize=8, color="tab:blue",
                va="bottom")
    else:
        ax.text(0.5, 0.9, "no laser peak: background only, no FWHM", transform=ax.transAxes, ha="center",
                color="tab:red", fontsize=13, fontweight="bold")
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, max(peak, c[m].max()) * 1.12)
    ax.set_xlabel("time after SYNC [ns]")
    ax.set_ylabel(f"counts per {bw * 1e3:.2f} ps bin" + (f" in {seconds:.1f} s" if seconds else ""))
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    if has_peak:
        ins = ax.inset_axes([0.66, 0.38, 0.32, 0.28])
        ins.step(t, np.clip(c, 1, None), where="mid", lw=0.4, color="0.45")
        ins.set_yscale("log")
        ins.set_xlim(t[0] - bw / 2, t[-1] + bw / 2)
        ins.set_ylim(max(1, 0.5 * np.percentile(c[c > 0], 5)) if (c > 0).any() else 1, peak * 3)
        ins.set_title("full period, log", fontsize=8)
        ins.tick_params(labelsize=7)
        ins.grid(alpha=0.3)
        ins.axvspan(lo, hi, color="tab:red", alpha=0.08)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=110)
    if show:
        plt.show()
    return fig, r


def plot_recording(run_dir, view_ns=(-0.7, 1.1), rebin_ps: float = 9.4, save: bool = True, show: bool = True):
    """``plot_histogram`` of a ``record_felix`` folder; saved there as histogram_fwhm.png."""
    run_dir = Path(run_dir)
    d = load_felix(run_dir)
    left_s = float(np.load(run_dir / "felix.npz")["left_s"])
    title = f"{d['meta'].get('condition', '')}  ({run_dir.parent.name}/{run_dir.name})"
    fig, r = plot_histogram(d["bins"].sum(0), d["bin_s"], title, seconds=d["seconds"], left_s=left_s,
                            view_ns=view_ns, rebin_ps=rebin_ps,
                            save_to=run_dir / "histogram_fwhm.png" if save else None, show=show)
    return r


def load_felix(run_dir) -> dict:
    run_dir = Path(run_dir)
    z = np.load(run_dir / "felix.npz")
    meta = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    return dict(dir=run_dir, bins=z["bins"], window_s=float(z["window_s"]), bin_s=float(z["bin_s"]),
                seconds=z["bins"].shape[0] * float(z["window_s"]), meta=meta)


def blank_felix(recs: dict, leaf: str = "leaf", blank: str = "no leaf") -> dict:
    """Table of the Felix recordings, each relative to the ``leaf`` one, and leaf minus blank.

    ``recs``: condition -> run folder. Rates come from decays rotated on the leaf's laser peak,
    so a blank without a peak of its own is counted at the same micro-time. Peak position, FWHM,
    mean arrival and tau1 come from each recording's own peak (``shapes``): a scatter peak sits
    ~0.2 ns before a fluorescence peak.
    """
    data = {name: load_felix(d) for name, d in recs.items()}
    ref_peak, warnings = None, []
    if leaf in data:
        L = data[leaf]
        s = decay_summary(L["bins"].sum(0), L["bin_s"], L["seconds"])
        ref_peak = s["peak_bin"] if s["significant"] else None
        if ref_peak is None:
            warnings.append(f"the '{leaf}' recording has no laser peak: the laser was blocked or does not reach "
                            f"the leaf. Repeat it with the laser on; ratios to it and '{leaf} - {blank}' are meaningless")
    rows, stats, shapes = [], {}, {}
    shape_keys = ("peak_ns", "fwhm_ps", "t_mean_ns", "t_mean_err_ns", "tau1_ns", "tau_int2_ns")
    for name, d in data.items():
        s = decay_summary(d["bins"].sum(0), d["bin_s"], d["seconds"], peak_bin=ref_peak)
        own = decay_summary(d["bins"].sum(0), d["bin_s"], d["seconds"])
        stats[name], shapes[name] = s, own
        row = dict(condition=name, **{k: v for k, v in _public(s).items() if not isinstance(v, (list, np.ndarray))})
        row.update({k: (own[k] if own["significant"] else float("nan")) for k in shape_keys})
        rows.append(row)
    if leaf in stats:
        for r in rows:
            r["laser_vs_leaf"] = r["laser_cps"] / stats[leaf]["laser_cps"]
            r["floor_vs_leaf"] = r["floor_cps"] / stats[leaf]["floor_cps"]
    corrected = None
    if leaf in data and blank in data:
        # only the blank's laser-correlated part is subtracted (its flat floor stays in the leaf's
        # own floor), so the leaf's tail keeps its raw Poisson counts and the fit weights stay right
        L, B = data[leaf], data[blank]
        hb = B["bins"].sum(0).astype(float) - stats[blank]["floor_per_bin"] / REBIN
        h = L["bins"].sum(0).astype(float) - hb * (L["seconds"] / B["seconds"])
        corrected = decay_summary(h, L["bin_s"], L["seconds"], peak_bin=ref_peak)
    if ref_peak is None:
        corrected = None
    return dict(rows=rows, stats=stats, shapes=shapes, corrected=corrected, leaf=leaf, blank=blank,
                ref_peak=ref_peak, warnings=warnings)


def blank_felix_report(b: dict) -> str:
    out = [f"{'condition':<26}{'total kcps':>11}{'laser kcps':>16}{'floor kcps':>11}{'laser/leaf':>11}"
           f"{'peak ns':>9}{'FWHM ps':>9}{'t_mean ns':>16}{'tau1 ns':>9}"]
    for r in b["rows"]:
        out.append(f"{r['condition']:<26}{r['total_cps'] / 1e3:>11.2f}"
                   f"{r['laser_cps'] / 1e3:>9.2f} +-{r['laser_cps_err'] / 1e3:>5.2f}{r['floor_cps'] / 1e3:>11.2f}"
                   f"{r.get('laser_vs_leaf', float('nan')):>11.3f}{r['peak_ns']:>9.2f}{r['fwhm_ps']:>9.0f}"
                   f"{r['t_mean_ns']:>8.3f} +-{r['t_mean_err_ns']:>6.3f}{r['tau1_ns']:>9.3f}"
                   + ("" if r["significant"] else "   (no laser peak)"))
    out += [f"WARNING: {w}" for w in b.get("warnings", [])]
    c = b["corrected"]
    if c is not None:
        out.append(f"{b['leaf'] + ' - ' + b['blank']:<26}{'':>11}{c['laser_cps'] / 1e3:>9.2f} +-"
                   f"{c['laser_cps_err'] / 1e3:>5.2f}{'':>11}{'':>11}{'':>9}{'':>9}{c['t_mean_ns']:>8.3f} +-"
                   f"{c['t_mean_err_ns']:>6.3f}{c['tau1_ns']:>9.3f}")
    return "\n".join(out)


def plot_blank_felix(b: dict, irf: dict | None = None, save_to=None, show: bool = True):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    for name, s in b["stats"].items():
        t = (np.arange(s["h"].size) + 0.5) * s["bin_ns"]
        ax.semilogy(t, np.maximum(s["h"] / s["seconds"], 1e-2), lw=0.9,
                    label=f"{name}: laser {s['laser_cps'] / 1e3:.2f} kcps, floor {s['floor_cps'] / 1e3:.1f} kcps")
    ax.set_xlim(0, 12.5)
    ax.set_xlabel("micro-time, laser peak at 2.5 ns [ns]")
    ax.set_ylabel(f"counts per second per {b['stats'][next(iter(b['stats']))]['bin_ns'] * 1e3:.0f} ps bin")
    ax.set_title("Felix decays, absolute (same rotation for all)")
    ax.legend(fontsize=7)
    ax = axes[1]
    curves = []
    if b["leaf"] in b["stats"]:
        curves.append((b["leaf"], b["stats"][b["leaf"]]))
    if b["corrected"] is not None:
        curves.append((f"{b['leaf']} - {b['blank']}", b["corrected"]))
    if b["blank"] in b["shapes"] and b["shapes"][b["blank"]]["significant"]:
        curves.append((b["blank"] + " (own peak)", b["shapes"][b["blank"]]))
    for name, s in curves:
        ax.semilogy(*_shape(s), lw=0.9, label=f"{name}: t_mean {s['t_mean_ns']:.3f} ns, tau1 {s['tau1_ns']:.3f} ns")
    if irf is not None:
        ax.semilogy(*_irf_curve(irf), "k:", lw=0.8, label=f"IRF ({irf['name'][:28]})")
    ax.set_xlim(0, 12.5)
    ax.set_ylim(1e-3, 1.5)
    ax.set_xlabel("micro-time [ns]")
    ax.set_ylabel("floor-subtracted, normalised")
    ax.set_title("shapes")
    ax.legend(fontsize=7)
    for a in axes:
        a.grid(alpha=0.3)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=120)
    if show:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# AMBIT windows of a chopped run (ambit_chop.run_chopped folders)
# ---------------------------------------------------------------------------
def _amb_rows(run_dir) -> dict:
    rows = list(csv.DictReader(open(Path(run_dir) / "ambit.csv")))
    col = lambda k, f=float: np.array([f(r[k]) for r in rows])
    return dict(line=col("line", int), k=col("k", int), s=col("s_630"), r=col("r_630"), ratio=col("ratio"),
                actinic=col("actinic", int), freq=col("freq", int))


def _stat(x: np.ndarray) -> dict:
    x = x[np.isfinite(x)]
    return dict(n=int(x.size), mean=float(x.mean()) if x.size else float("nan"),
                sem=float(x.std(ddof=1) / np.sqrt(x.size)) if x.size > 1 else float("nan"))


def ambit_windows(run_dir, low_last_s: float = 1.0, fm_skip_s: float = 0.25) -> dict:
    """AMBIT s, r and s/r in the windows ``ambit_chop.analyse`` uses.

    Pulse = the first line with actinic >= 200. ``low`` = the line before it: Fo when dark,
    Fs when lit; its last ``low_last_s``. Fm = the pulse from ``fm_skip_s`` on. ``after`` = the
    line after the pulse. Also the mean of every line.
    """
    meta = json.loads((Path(run_dir) / "run.json").read_text(encoding="utf-8"))
    prot = meta["protocol"]
    a = _amb_rows(run_dir)
    tk = a["k"] / np.maximum(a["freq"], 1)
    p = next((i for i, s in enumerate(prot) if s["actinic"] >= 200), int(np.argmax([s["actinic"] for s in prot])))
    low = p - 1
    win = {}
    if low >= 0:
        dur = prot[low]["pulses"] / prot[low]["freq"]
        name = "Fo" if prot[low]["actinic"] <= ag.ACTINIC_ON else "Fs"
        m = (a["line"] == low) & (tk >= dur - low_last_s)
        win[name] = dict(line=low, **{q: _stat(a[q][m]) for q in ("s", "r", "ratio")})
    m = (a["line"] == p) & (tk >= fm_skip_s)
    win["Fm"] = dict(line=p, **{q: _stat(a[q][m]) for q in ("s", "r", "ratio")})
    if p + 1 < len(prot):
        m = a["line"] == p + 1
        win["after"] = dict(line=p + 1, **{q: _stat(a[q][m]) for q in ("s", "r", "ratio")})
    lines = [dict(line=i, actinic=s["actinic"], freq=s["freq"],
                  **{q: _stat(a[q][a["line"] == i]) for q in ("s", "r", "ratio")}) for i, s in enumerate(prot)]
    return dict(run_dir=str(run_dir), protocol=prot, windows=win, lines=lines,
                low=next((n for n in ("Fo", "Fs") if n in win), None))


def _same_protocol(a: list, b: list) -> bool:
    key = lambda p: [(s["pulses"], s["freq"], s["actinic"]) for s in p]
    return key(a) == key(b)


def ambit_blank_correct(leaf_run, blank_run, felix_blank_laser_cps: float | None = None) -> dict:
    """Fv/Fm or Phi2 of a leaf run before and after subtracting the blank run.

    AMBIT: s/r of the blank (same protocol) is subtracted window by window.
    Felix (optional): ``felix_blank_laser_cps``, the laser-correlated rate with no leaf, is
    subtracted from the leaf run's gap intensities (chop_results.json).
    """
    L, B = ambit_windows(leaf_run), ambit_windows(blank_run)
    low = L["low"]
    out = dict(leaf=str(leaf_run), blank=str(blank_run), low=low,
               same_protocol=_same_protocol(L["protocol"], B["protocol"]))
    lo, fm = L["windows"][low]["ratio"]["mean"], L["windows"]["Fm"]["ratio"]["mean"]
    blo, bfm = B["windows"].get(low, {}).get("ratio", {}).get("mean", float("nan")), B["windows"]["Fm"]["ratio"]["mean"]
    out["ambit"] = dict(leaf_low=lo, leaf_Fm=fm, blank_low=blo, blank_Fm=bfm,
                        blank_share_low=blo / lo, blank_share_Fm=bfm / fm,
                        raw=1 - lo / fm, corrected=1 - (lo - blo) / (fm - bfm))
    res_file = Path(leaf_run) / "chop_results.json"
    if res_file.exists():
        r = json.loads(res_file.read_text(encoding="utf-8"))
        w = r.get("windows", {})
        if low in w and "Fm" in w:
            il, im = w[low]["laser_cps"], w["Fm"]["laser_cps"]
            fe = dict(leaf_low_cps=il, leaf_Fm_cps=im, raw=1 - il / im)
            est = r.get("estimates", {}).get("fvfm" if low == "Fo" else "phi2", {})
            fe["tau_shared"] = est.get("tau_shared", float("nan"))
            fe["tau_1exp"] = est.get("tau_1exp", float("nan"))
            if felix_blank_laser_cps is not None:
                fe["blank_cps"] = felix_blank_laser_cps
                fe["corrected"] = 1 - (il - felix_blank_laser_cps) / (im - felix_blank_laser_cps)
            out["felix"] = fe
    return out


def ambit_lines_report(run_dir) -> str:
    w = ambit_windows(run_dir)
    out = [f"{Path(run_dir).name}"]
    for ln in w["lines"]:
        out.append(f"  line {ln['line']}  actinic {ln['actinic']:>3}  {ln['freq']:>3} Hz  n {ln['s']['n']:>3}  "
                   f"s {ln['s']['mean']:8.1f}  r {ln['r']['mean']:8.1f}  s/r {ln['ratio']['mean']:.4f} "
                   f"+- {ln['ratio']['sem']:.4f}")
    for name, x in w["windows"].items():
        out.append(f"  window {name:<5} (line {x['line']}): s {x['s']['mean']:8.1f}  s/r {x['ratio']['mean']:.4f} "
                   f"+- {x['ratio']['sem']:.4f}  (n {x['ratio']['n']})")
    return "\n".join(out)


def spad_rate(run_dir, bin_s: float = 1e-3) -> tuple[np.ndarray, np.ndarray]:
    """SPAD counts per second in ``bin_s`` bins from a run's ttm.bin: (t_s, rate)."""
    head, rec = ac.read_ttm(Path(run_dir) / "ttm.bin")
    lsb = head["lsb_s"]
    acc = np.zeros(0)
    t0 = None
    for ch, t in ac._chunks(rec):
        sp = t[ch == ag.SPAD_CH]
        if not sp.size:
            continue
        if t0 is None:
            t0 = int(sp[0])
        c = np.bincount(((sp - t0) * (lsb / bin_s)).astype(np.int64))
        if c.size > acc.size:
            acc = np.r_[acc, np.zeros(c.size - acc.size)]
        acc[:c.size] += c
    return (np.arange(acc.size) + 0.5) * bin_s, acc / bin_s


# ---------------------------------------------------------------------------
# Laser dose: AMBIT + Felix histograms around a manual laser unblock
# ---------------------------------------------------------------------------
def dose_protocol(dark_s=5.0, laser_s=20.0, pulse_s=0.8, post_s=2.0, f_ambit=20) -> list[dict]:
    """Line 0: measuring pulses only, laser blocked for ``dark_s`` then unblocked by hand;
    line 1: saturating pulse; line 2: dark."""
    return [ac.line(dark_s + laser_s, f_ambit, 0), ac.line(pulse_s, 100, 255), ac.line(post_s, 10, 0)]


def _beep():
    try:
        import winsound
        winsound.Beep(1500, 400)
    except Exception:
        pass


def wait_minutes(minutes: float) -> None:
    """Dark-adaptation timer, prints the remaining time every minute (interrupt to stop)."""
    end = time.monotonic() + minutes * 60
    while (left := end - time.monotonic()) > 0:
        print(f"  {left / 60:4.1f} min left", flush=True)
        time.sleep(min(60.0, left))
    print("dark adaptation done", flush=True)
    _beep()


def run_dose(ambit: ag.AmbitLink, felix: ag.FelixLink, out_dir, label: str, dark_s: float = 5.0,
             laser_s: float = 20.0, pulse_s: float = 0.8, post_s: float = 2.0, f_ambit: int = 20,
             window_s: float = 0.02, cue: bool = True, note: str = "", pre_s: float = 1.0,
             tail_s: float = 1.0, threshold_v: float = ag.THRESHOLD_V,
             laser_hz: float = ag.LASER_REP_RATE_HZ, log=print) -> Path:
    """One laser-dose run. The laser must be blocked when this starts.

    The Felix records contiguous histograms of ``window_s`` for the whole run. The AMBIT runs
    ``dose_protocol``; ``dark_s`` after its first sample the cell prints UNBLOCK and beeps.
    The real laser-on moment is measured afterwards from the Felix. ``cue=False`` for the
    control run with the laser blocked throughout.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    protocol = dose_protocol(dark_s, laser_s, pulse_s, post_s, f_ambit)
    cmd = ag.build_arrun(protocol)
    duration = ag.protocol_duration(protocol)
    try:
        ac.disarm_chop(ambit)        # no gaps: a chop armed and never consumed would apply here
    except Exception:
        pass
    felix.set_thresholds(threshold_v)
    rates = felix.count_rates(0.5)
    log("Felix input rates: " + ", ".join(f"ch{c}={r:,.0f}/s" for c, r in rates.items()))
    felix.configure_histogram(ag.SYNC_CH, ag.SPAD_CH, 0.0, 1.0 / laser_hz)
    n_win = int(np.ceil((pre_s + ag.NOMINAL_LATENCY_S + duration + tail_s + 0.5) / window_s))
    felix.start_windows(window_s, n_win)
    time.sleep(pre_s)
    result: dict = {}

    def worker():
        try:
            result["data"], result["t_send"], result["t_done"] = ag.run_arrun(ambit, cmd, protocol, log=log)
        except Exception as exc:            # re-raised in the main thread
            result["error"] = exc

    th = threading.Thread(target=worker, daemon=True)
    log("AMBIT: " + cmd.strip())
    t_cmd = time.monotonic()
    th.start()
    t_cue = t_cmd + ag.NOMINAL_LATENCY_S + dark_s
    if cue:
        for k in (3, 2, 1):
            while time.monotonic() < t_cue - k:
                time.sleep(0.01)
            print(f"   unblock in {k} ...", flush=True)
        while time.monotonic() < t_cue:
            time.sleep(0.005)
        print(">>> UNBLOCK THE LASER NOW <<<", flush=True)
        _beep()
    th.join(timeout=duration + 30)
    if "error" in result:
        felix.wait_done(timeout_s=1.0)
        raise result["error"]
    if "data" not in result:
        felix.wait_done(timeout_s=1.0)
        raise TimeoutError("no AMBIT data after the run")
    if not felix.wait_done(timeout_s=tail_s + 10.0):
        log(f"WARNING: the Felix delivered {len(felix.windows)} of {n_win} windows")
    res = ag.RunResult(protocol=protocol, cmd=cmd, ambit=result["data"], windows=list(felix.windows),
                       window_s=felix.window_s, t_soa=felix.t_soa, t_send=result["t_send"],
                       t_done=result["t_done"], bin_s=felix.bin_s, left_s=felix.left_s, laser_hz=laser_hz,
                       rates_before=rates, ambit_name=ambit.name, ambit_fw=ambit.fw,
                       felix_desc=felix.description)
    ag.align(res, log=log)
    save_dose(out_dir, res, label=label, cue_s=t_cue - res.run_start if cue else float("nan"),
              dark_s=dark_s, laser_s=laser_s, note=note)
    log(f"saved to {out_dir}")
    return out_dir


def save_dose(out_dir: Path, res: ag.RunResult, **meta) -> None:
    """felix.npz (windows, start times since AMBIT sample 0), ambit.csv, run.json."""
    out_dir = Path(out_dir)
    bins = res.bins_matrix().astype(np.uint32)
    np.savez_compressed(out_dir / "felix.npz", bins=bins, t_start_s=res.win_start - res.run_start,
                        window_s=res.window_s, bin_s=res.bin_s, left_s=res.left_s)
    t_amb = res.ambit_time() - res.run_start
    s = np.asarray(res.ambit.get("s_630", np.zeros(0)), float)
    r = np.asarray(res.ambit.get("r_630", np.zeros(0)), float)
    n = min(t_amb.size, s.size, r.size)
    with open(out_dir / "ambit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["line", "k", "t_s", "s_630", "r_630", "ratio", "actinic", "freq"])
        i = 0
        for li, seg in enumerate(res.protocol):
            for k in range(seg["pulses"]):
                if i >= n:
                    break
                wr.writerow([li, k, f"{t_amb[i]:.5f}", int(s[i]), int(r[i]),
                             f"{s[i] / r[i]:.6f}" if r[i] > 0 else "nan", seg["actinic"], seg["freq"]])
                i += 1
    info = dict(kind="dose", protocol=res.protocol, cmd=res.cmd.strip(), window_s=res.window_s,
                bin_s=res.bin_s, laser_hz=res.laser_hz, align_mode=res.align_mode,
                latency_s=res.run_start - res.t_send, rates_before={str(k): v for k, v in res.rates_before.items()},
                lost=int(sum(w.lost for w in res.windows)), ambit=dict(name=res.ambit_name, fw=res.ambit_fw),
                felix=res.felix_desc, saved=time.strftime("%Y-%m-%d %H:%M:%S"), **meta)
    (out_dir / "run.json").write_text(json.dumps(info, indent=1, default=float), encoding="utf-8")


def _laser_on(lc: np.ndarray, ts: np.ndarray, W: float, level: float, usable: np.ndarray, noise1: float):
    """Laser-on time from the per-window laser-correlated counts ``lc``, or None; +-1 window.

    Step search on 5-window running sums (100 ms at 20 ms windows), above 20 % of the steady
    ``level`` and 6 sigma of the laser-blocked noise (``noise1`` per window): a leaf that
    starts much dimmer than its steady state, or a weak laser, is still found at its first
    photons. The edge is then placed in the first window above half of the level just after it.
    """
    k5 = np.convolve(np.where(usable, lc, 0.0), np.ones(5), mode="full")[4:]   # sum of windows i..i+4
    hits = np.flatnonzero(usable & (k5 > max(0.2 * 5 * level, 6 * noise1 * np.sqrt(5))))
    if not hits.size:
        return None, False
    i0 = int(hits[0])
    already = bool(i0 == np.flatnonzero(usable)[0])
    nxt = lc[i0 + 5:i0 + 10]
    local = float(np.median(nxt)) if nxt.size else level           # level right after the step
    i = next((j for j in range(i0, min(i0 + 5, lc.size)) if lc[j] > 0.5 * local), i0)
    frac = float(np.clip(lc[i] / max(local, 1e-9), 0.0, 1.0))
    if i > 0 and frac > 0.9 and lc[i - 1] > max(0.1 * local, 3 * noise1):   # edge in the previous window
        i -= 1
        frac = float(np.clip(lc[i] / max(local, 1e-9), 0.0, 1.0))
    return float(ts[i] + (1.0 - frac) * W), already


def analyse_dose(run_dir, steady_s: float = 5.0, initial_s: float | None = None, fm_skip_s: float = 0.25,
                 initial_photons: float = 5000.0) -> dict:
    """Laser-on time, Felix kinetics after it, initial vs steady decays, AMBIT yields.

    ``steady_s``: the "steady" decay and the AMBIT F_ss are the last ``steady_s`` before the pulse.
    ``initial_s``: length of the Felix "first" decay after laser-on. None: long enough for
    ``initial_photons`` laser photons at the steady rate, 0.1 .. 1 s (0.1 s at ~90 kcps, ~0.5 s
    at a 10x weaker laser), so its mean arrival is good to ~0.05 ns.
    The Felix windows inside the saturating pulse are summarised apart (``pulse``), without
    dead-time correction: the SPAD is usually far above its rating there.
    """
    run_dir = Path(run_dir)
    meta = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    z = np.load(run_dir / "felix.npz")
    bins = z["bins"].astype(float)
    ts = z["t_start_s"].astype(float)
    W, bin_s = float(z["window_s"]), float(z["bin_s"])
    prot = meta["protocol"]
    t_pulse = prot[0]["pulses"] / prot[0]["freq"]
    tc = ts + W / 2
    before_pulse = ts + W <= t_pulse - 0.02
    warnings = []
    if meta.get("align_mode") != "flash-edges":
        warnings.append(f"AMBIT time axis from the nominal latency ({meta.get('align_mode')}), "
                        "not from the flash edges: AMBIT vs Felix timing uncertain by ~10 ms")
    # steady reference decay, its peak fixes the rotation of every window
    ss = before_pulse & (ts >= t_pulse - 0.02 - steady_s)
    if not ss.any():
        raise ValueError(f"{run_dir.name}: no Felix windows in the last {steady_s} s before the pulse")
    steady = decay_summary(bins[ss].sum(0), bin_s, ss.sum() * W)
    t_on, detected, already = float(meta.get("cue_s", float("nan"))), False, False
    lc = np.full(ts.size, np.nan)
    if steady["significant"]:
        pk = steady["peak_bin"]
        hr = fl._rebin_matrix(bins, REBIN)
        bw = bin_s * REBIN * 1e9
        p0 = int(round(PEAK_AT_NS / bw))
        hr = np.roll(hr, p0 - pk, axis=1)
        t = (np.arange(hr.shape[1]) + 0.5) * bw
        pre = (t > 0.3) & (t < PEAK_AT_NS - 0.6)
        rise = p0 - int(round(RISE_NS / bw))
        lc = (hr[:, rise:] - hr[:, pre].mean(1)[:, None]).sum(1)
        level = float(np.median(lc[ss]))
        # noise per window with the laser surely blocked: from the Felix start to 1 s before the cue
        cue = meta.get("cue_s")
        cue = float(cue) if cue is not None and np.isfinite(cue) else float(meta.get("dark_s", 5.0))
        quiet = before_pulse & (ts + W < cue - 1.0)
        q = lc[quiet]
        noise1 = float(1.4826 * np.median(np.abs(q - np.median(q)))) if q.size >= 10 else np.sqrt(max(level, 1.0))
        found, already = _laser_on(lc, ts, W, level, before_pulse, max(noise1, 1.0))
        if found is not None:
            t_on, detected = found, True
        if already:
            warnings.append("the laser was already on when the recording started: no dark reference")
    else:
        pk = None
        warnings.append("no laser peak before the pulse (laser blocked throughout?): Felix not analysed, "
                        "AMBIT times relative to the cue")
    if not np.isfinite(t_on):
        t_on = float(meta.get("dark_s", 5.0))
        warnings.append("no laser-on time: using dark_s as time zero")
    rel = tc - t_on
    after = before_pulse & (ts >= t_on)
    kin = []
    if pk is not None:
        for a, b in zip(KIN_EDGES_S[:-1], KIN_EDGES_S[1:]):
            m = before_pulse & (ts >= 0) & (rel >= a) & (rel < b)
            if not m.any():
                continue
            s = decay_summary(bins[m].sum(0), bin_s, m.sum() * W, peak_bin=pk)
            kin.append(dict(t_s=float(rel[m].mean()), n=int(m.sum()), laser_cps=s["laser_cps"],
                            laser_cps_err=s["laser_cps_err"], floor_cps=s["floor_cps"],
                            t_mean_ns=s["t_mean_ns"], t_mean_err_ns=s["t_mean_err_ns"], tau1_ns=s["tau1_ns"]))
    if initial_s is None:
        initial_s = float(np.clip(initial_photons / steady["laser_cps"], 0.1, 1.0)) if pk is not None else 0.1
    init_m = after & (ts + W <= t_on + initial_s + 1e-9)
    initial = decay_summary(bins[init_m].sum(0), bin_s, init_m.sum() * W, peak_bin=pk) \
        if (init_m.any() and pk is not None) else None
    bg_m = before_pulse & (ts + W <= t_on) & (ts >= 0)
    background = decay_summary(bins[bg_m].sum(0), bin_s, bg_m.sum() * W, peak_bin=pk) if bg_m.any() else None
    # saturating pulse: windows fully inside it, raw rates (the dead-time model does not hold there)
    pulse_len = prot[1]["pulses"] / prot[1]["freq"] if len(prot) > 1 else 0.0
    in_p = (ts >= t_pulse + 0.01) & (ts + W <= t_pulse + pulse_len - 0.01)
    pulse = decay_summary(bins[in_p].sum(0), bin_s, in_p.sum() * W, peak_bin=pk, dead_time_s=0.0) \
        if (in_p.any() and pk is not None) else None
    if pulse is not None and pulse["total_cps"] > SPAD_MAX_CPS:
        warnings.append(f"the SPAD counted {pulse['total_cps'] / 1e6:.1f} Mcps in the saturating pulse, above its "
                        f"rated maximum ({SPAD_MAX_CPS / 1e6:.1f} Mcps, max 50 uW incident): more ND in front of "
                        "the SPAD, or no pulse")
    felix_total = bins.sum(1) / W
    # AMBIT
    rows = list(csv.DictReader(open(run_dir / "ambit.csv")))
    at = np.array([float(r["t_s"]) for r in rows])
    ar = np.array([float(r["ratio"]) for r in rows])
    al = np.array([int(r["line"]) for r in rows])
    st = lambda m: _stat(ar[m])
    dark = st((al == 0) & (at >= max(0.2, t_on - 3.0)) & (at < t_on - 0.05))
    first = st((al == 0) & (at >= t_on + 0.05) & (at < t_on + 0.35))
    sss = st((al == 0) & (at >= t_pulse - steady_s) & (at < t_pulse - 0.01))
    fm = st((al == 1) & (at >= t_pulse + fm_skip_s))
    amb = dict(F_dark=dark, F_first=first, F_ss=sss, Fm=fm,
               fvfm_dark=1 - dark["mean"] / fm["mean"], phi2_laser=1 - sss["mean"] / fm["mean"],
               step=first["mean"] / dark["mean"] - 1, rise=sss["mean"] / dark["mean"] - 1)
    return dict(run_dir=str(run_dir), label=meta.get("label", run_dir.name), note=meta.get("note", ""),
                t_on=t_on, detected=detected, warnings=warnings, window_s=W, t_pulse=t_pulse - t_on,
                pulse_s=pulse_len, steady_s=steady_s,
                kin={k: np.array([r[k] for r in kin]) for k in (kin[0] if kin else {})},
                initial=initial, steady=steady if pk is not None else None, background=background, pulse=pulse,
                ambit=amb, amb_t=at - t_on, amb_ratio=ar, amb_line=al, felix_t=rel, felix_lc_cps=lc / W,
                felix_total_cps=felix_total)


def dose_report(d: dict) -> str:
    a = d["ambit"]
    f = lambda s: f"{s['mean']:.4f} +- {s['sem']:.4f} (n {s['n']})"
    out = [f"{d['label']}  ({Path(d['run_dir']).name})",
           f"  laser on at {d['t_on']:.3f} s after AMBIT sample 0" + ("" if d["detected"] else " (NOT detected)")]
    out += [f"  WARNING: {w}" for w in d["warnings"]]
    ss = d.get("steady_s", 3.0)
    out += [f"  AMBIT s/r  before laser {f(a['F_dark'])}",
            f"             first 0.05-0.35 s with laser {f(a['F_first'])}   step {100 * a['step']:+.1f} %",
            f"             last {ss:g} s with laser {f(a['F_ss'])}   rise {100 * a['rise']:+.1f} %",
            f"             Fm {f(a['Fm'])}",
            f"  AMBIT Fv/Fm from the dark level {a['fvfm_dark']:.3f};  1 - F(laser)/Fm {a['phi2_laser']:.3f}"]
    b = d["background"]
    if b is not None:
        out.append(f"  Felix before laser-on: {b['total_cps'] / 1e3:.1f} kcps total, laser-correlated "
                   f"{b['laser_cps'] / 1e3:.2f} +- {b['laser_cps_err'] / 1e3:.2f} kcps")
    for name in ("initial", "steady"):
        s = d[name]
        if s is not None:
            out.append(f"  Felix {name:<8} ({s['seconds']:.2f} s): laser {s['laser_cps'] / 1e3:6.1f} kcps, "
                       f"t_mean {s['t_mean_ns']:.3f} +- {s['t_mean_err_ns']:.3f} ns, tau1 {s['tau1_ns']:.3f} ns")
    p = d.get("pulse")
    if p is not None:
        lp = (f"t_mean {p['t_mean_ns']:.3f} +- {p['t_mean_err_ns']:.3f} ns" if p["significant"]
              else "no laser peak")
        out.append(f"  Felix in the pulse ({p['seconds']:.2f} s): {p['total_cps'] / 1e6:.2f} Mcps total, laser-correlated "
                   f"{p['laser_cps'] / 1e3:.0f} +- {p['laser_cps_err'] / 1e3:.0f} kcps (raw), {lp}")
    return "\n".join(out)


def plot_fits(decays: list, irf: dict | None = None, title: str = "", save_to=None, show: bool = True,
              models=(1, 2), ylim=None):
    """Raw counts and tail fits, one column per decay: ``decays`` = [(name, decay_summary, colour), ...].

    Top: counts per 47 ps bin as recorded (rotated so the laser peak sits at 2.5 ns), the flat
    floor, the tail fits over their range and the IRF scaled to the peak.
    Bottom: residuals in Poisson sigma, (data - fit) / sqrt(data), as the fit weights them.
    ``models``: which fits to draw, 1- and / or 2-exponential. The 2-exp reports its lifetimes,
    their amplitudes as fitted (counts per 47 ps bin at the laser peak) and
    t_avg = (tau1 a1 + tau2 a2) / (a1 + a2).
    ``ylim``: the same (log) count range in every panel, e.g. (1, 1e5); None: per panel.
    """
    import matplotlib.pyplot as plt

    decays = [x for x in decays if x[1] is not None and x[1]["significant"]]
    if not decays:
        return None
    n = len(decays)
    fig, axes = plt.subplots(2, n, figsize=(6.2 * n, 6.6), squeeze=False, sharex=True,
                             gridspec_kw=dict(height_ratios=[3, 1]))
    irf_c = _irf_curve(irf) if irf is not None else None
    for j, (name, s, col) in enumerate(decays):
        ax, axr = axes[0, j], axes[1, j]
        t = (np.arange(s["h"].size) + 0.5) * s["bin_ns"]
        h = s["h"]
        ax.step(t, np.maximum(h, 0.5), where="mid", color=col, lw=0.8, label="data (47 ps bins)")
        ax.axhline(s["floor_per_bin"], color="0.5", lw=0.8, ls="--", label=f"floor {s['floor_per_bin']:.1f} /bin")
        f = s.get("fit")
        txt = [f"{s['seconds']:.2f} s, {s['laser_cps'] * s['seconds'] / s['dead_time_factor']:,.0f} laser counts",
               f"t_mean {s['t_mean_ns']:.3f} +- {s['t_mean_err_ns']:.3f} ns"]
        if f is not None:
            sig = np.sqrt(np.maximum(h[-f["t"].size:], 1.0))
            c = h[-f["t"].size:]
            ax.axvline(f["t"][0], color="0.3", lw=0.6, ls=":")
            if 1 in models:
                ax.plot(f["t"], f["model1"], color="k", lw=1.0, label=f"1-exp: tau {s['tau1_ns']:.3f} ns")
                txt.append(f"1-exp: tau {s['tau1_ns']:.3f} ns, chi2r {s['chi2_1']:.2f}")
                axr.plot(f["t"], (c - f["model1"]) / sig, color="k", lw=0.7, label="1-exp")
            if 2 in models:
                (t1, t2), (a1, a2) = s["taus2_ns"], s["amps2"]
                ax.plot(f["t"], f["model2"], color="tab:cyan", lw=1.2, ls="--",
                        label=f"2-exp: t_avg {s['tau_avg2_ns']:.3f} ns")
                txt += [f"2-exp: tau1 {t1:.3f} ns, a1 {a1:,.0f};  tau2 {t2:.3f} ns, a2 {a2:,.0f}",
                        "       (a = counts per 47 ps bin at the peak)",
                        f"t_avg = (tau1*a1 + tau2*a2)/(a1 + a2) = {s['tau_avg2_ns']:.3f} ns, chi2r {s['chi2_2']:.2f}"]
                if t1 < MIN_FAST_NS:
                    txt.append(f"tau1 < {MIN_FAST_NS} ns: below what a tail fit resolves")
                axr.plot(f["t"], (c - f["model2"]) / sig, color="tab:cyan", lw=0.7, label="2-exp")
        else:
            txt.append("too few laser counts for the tail fits")
        if irf_c is not None:
            pk = float(h.max() - s["floor_per_bin"])
            ax.plot(irf_c[0], irf_c[1] * pk + s["floor_per_bin"], "k:", lw=0.8, label="IRF (scaled)")
        ax.set_yscale("log")
        ax.set_ylim(*(ylim or (max(0.5, 0.3 * s["floor_per_bin"]), 2.0 * h.max())))
        ax.set_title(f"{title}{': ' if title else ''}{name}")
        ax.set_ylabel("counts per 47 ps bin")
        ax.text(0.98, 0.97, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                family="monospace", bbox=dict(fc="white", ec="0.8", alpha=0.9))
        ax.legend(fontsize=7.5, loc="lower left")
        axr.axhline(0, color="0.5", lw=0.6)
        axr.set_ylim(-6, 6)
        axr.set_xlim(0, 12.5)
        axr.set_xlabel("micro-time [ns] (laser peak at 2.5 ns)")
        axr.set_ylabel("residual [sigma]")
        if f is not None:
            axr.legend(fontsize=7.5, loc="upper right", ncol=2)
        for x in (ax, axr):
            x.grid(alpha=0.3)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=120)
    if show:
        plt.show()
    return fig


def plot_dose(d: dict, irf: dict | None = None, save_to=None, show: bool = True, fits: bool = True):
    """The dose figure; ``fits``: also the raw decays with their tail fits (``plot_fits``), saved
    next to ``save_to`` as <name>_fits.png."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    a = d["ambit"]
    sl = f"last {d.get('steady_s', 3.0):g} s"
    ax = axes[0, 0]
    ax.plot(d["amb_t"], d["amb_ratio"], ".", ms=3, color="tab:blue")
    for key, c in (("F_dark", "0.4"), ("F_ss", "tab:orange"), ("Fm", "tab:red")):
        ax.axhline(a[key]["mean"], color=c, lw=0.8, ls="--", label=f"{key} {a[key]['mean']:.4f}")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("time since laser on [s]")
    ax.set_ylabel("AMBIT s/r")
    ax.set_title(f"AMBIT: {d['label']}")
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    ft, lcps = d["felix_t"], d["felix_lc_cps"]
    tp0, tp1 = d["t_pulse"], d["t_pulse"] + d.get("pulse_s", 0.0)
    inside = (ft > tp0) & (ft < tp1)
    m = np.isfinite(lcps) & ~inside
    ax.plot(ft[m], lcps[m] / 1e3, ".", ms=2, color="0.6", label=f"{d['window_s'] * 1e3:.0f} ms windows (raw)")
    p = d.get("pulse")
    if inside.any() and np.isfinite(lcps[inside]).any():
        ax.plot(ft[inside], lcps[inside] / 1e3, ".", ms=2, color="tab:orange", label="in the pulse (raw)")
    if tp1 > tp0:
        ax.axvspan(tp0, tp1, color="tab:orange", alpha=0.15, lw=0)
        if p is not None:
            msg = (f"pulse: SPAD {p['total_cps'] / 1e6:.1f} Mcps total\nlaser-correlated "
                   f"{p['laser_cps'] / 1e3:.0f} +- {p['laser_cps_err'] / 1e3:.0f} kcps"
                   + ("" if p["significant"] else "\n= no laser peak, SPAD blind"))
            ax.text(tp0 - 0.3, 0.55, msg, transform=ax.get_xaxis_transform(), ha="right", va="center",
                    fontsize=7.5, color="tab:orange", bbox=dict(fc="white", ec="tab:orange", alpha=0.9))
    k = d["kin"]
    if k:
        ax.errorbar(k["t_s"], k["laser_cps"] / 1e3, k["laser_cps_err"] / 1e3, fmt="o", ms=3, color="tab:red",
                    label="kinetics bins (dead-time corr.)")
        ax3 = ax.twinx()
        ax3.errorbar(k["t_s"], k["t_mean_ns"], k["t_mean_err_ns"], fmt="^", ms=4, color="tab:purple")
        ax3.set_ylabel("mean arrival after peak [ns]", color="tab:purple")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("time since laser on [s]")
    ax.set_ylabel("laser-correlated kcps", color="tab:red")
    ax.set_title("Felix")
    ax.legend(fontsize=8, loc="lower right")
    ax = axes[1, 0]
    if k:
        ss = d["steady"]
        pos = k["t_s"] > 0
        ax.semilogx(k["t_s"][pos], k["laser_cps"][pos] / ss["laser_cps"], "o-", ms=3, color="tab:red",
                    label="Felix laser-correlated rate / steady")
        ax.semilogx(k["t_s"][pos], k["t_mean_ns"][pos] / ss["t_mean_ns"], "^-", ms=3, color="tab:purple",
                    label="Felix mean arrival / steady")
    on = d["amb_line"] == 0
    rel = d["amb_t"]
    yb = []
    for lo, hi in zip(KIN_EDGES_S[:-1], KIN_EDGES_S[1:]):
        mm = on & (rel >= max(lo, 0.0)) & (rel < hi) & (rel < d["t_pulse"])
        if hi > 0 and mm.any():
            yb.append((rel[mm].mean(), np.nanmean(d["amb_ratio"][mm])))
    if yb:
        yb = np.array(yb)
        ax.semilogx(yb[:, 0], yb[:, 1] / a["F_ss"]["mean"], "s-", ms=3, color="tab:blue", label="AMBIT s/r / steady")
    ax.axhline(a["F_dark"]["mean"] / a["F_ss"]["mean"], color="tab:blue", ls="--", lw=0.8,
               label="AMBIT before laser-on / steady")
    ax.axhline(1, color="k", lw=0.5)
    ax.set_xlabel("time since laser on [s]")
    ax.set_ylabel(f"relative to the {sl}")
    ax.set_title("yield after laser-on, both instruments")
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    first = f"first {d['initial']['seconds']:.2f} s" if d["initial"] is not None else "first"
    for name, s, c in ((first, d["initial"], "tab:green"), (sl, d["steady"], "tab:red")):
        if s is None or not s["significant"]:
            continue
        ax.semilogy(*_shape(s), color=c, lw=0.9,
                    label=f"{name}: t_mean {s['t_mean_ns']:.3f} ns, tau1 {s['tau1_ns']:.3f} ns")
    if irf is not None:
        ax.semilogy(*_irf_curve(irf), "k:", lw=0.8, label="IRF")
    ax.set_xlim(0, 12.5)
    ax.set_ylim(1e-3, 1.5)
    ax.set_xlabel("micro-time [ns]")
    ax.set_ylabel("floor-subtracted, normalised")
    ax.set_title("Felix decays after laser-on")
    ax.legend(fontsize=8)
    for x in axes.flat:
        x.grid(alpha=0.3)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=120)
    if show:
        plt.show()
    if fits:
        fsave = None
        if save_to:
            fsave = Path(save_to).with_name(Path(save_to).stem + "_fits.png")
        plot_fits([(first, d["initial"], "tab:green"), (sl, d["steady"], "tab:red"),
                   ("in the pulse", d.get("pulse"), "tab:orange")], irf, title=d["label"],
                  save_to=fsave, show=show)
    return fig


def compare_dose(results: list[dict], save_to=None, show: bool = True):
    """Overlay of several dose runs: AMBIT yield and Felix rate / lifetime vs time since laser-on."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for d in results:
        a = d["ambit"]
        on = d["amb_line"] == 0
        axes[0].plot(d["amb_t"][on], d["amb_ratio"][on] / a["F_dark"]["mean"], ".", ms=3, label=d["label"])
        k = d["kin"]
        if k:
            pos = k["t_s"] > 0                  # log time axis: the unblocking ramp before 0 is left out
            axes[1].semilogx(k["t_s"][pos], k["laser_cps"][pos] / 1e3, "o-", ms=3, label=d["label"])
            axes[2].errorbar(k["t_s"][pos], k["t_mean_ns"][pos], k["t_mean_err_ns"][pos], fmt="^-", ms=3,
                             label=d["label"])
            axes[2].set_xscale("log")
    axes[0].axvline(0, color="k", lw=0.8)
    axes[0].axhline(1, color="k", lw=0.5)
    axes[0].set_xlabel("time since laser on [s]")
    axes[0].set_ylabel("AMBIT s/r / level before laser-on")
    axes[0].set_title("AMBIT")
    axes[1].set_xlabel("time since laser on [s]")
    axes[1].set_ylabel("laser-correlated kcps")
    axes[1].set_title("Felix intensity")
    axes[2].set_xlabel("time since laser on [s]")
    axes[2].set_ylabel("mean arrival after peak [ns]")
    axes[2].set_title("Felix lifetime (model free)")
    for x in axes:
        x.grid(alpha=0.3)
        x.legend(fontsize=8)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=120)
    if show:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# SPAD during the saturating pulse of a dose run
# ---------------------------------------------------------------------------
def pulse_spad(run_dir, before_s: float = 1.0, after_s: float = 1.5) -> dict:
    """SPAD rates and micro-time histograms around the saturating pulse of a dose run.

    The Felix data are 20 ms hardware histograms, so the time resolution is one window.
    ``ref`` is a stretch of the same length just before the pulse (laser on, no LED).
    """
    run_dir = Path(run_dir)
    meta = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    z = np.load(run_dir / "felix.npz")
    bins = z["bins"].astype(float)
    ts = z["t_start_s"].astype(float)
    W, bin_s = float(z["window_s"]), float(z["bin_s"])
    prot = meta["protocol"]
    if len(prot) < 2:
        raise ValueError(f"{run_dir.name}: no pulse line in the protocol")
    tp0 = prot[0]["pulses"] / prot[0]["freq"]
    tp1 = tp0 + prot[1]["pulses"] / prot[1]["freq"]
    in_p = (ts >= tp0 + 0.01) & (ts + W <= tp1 - 0.01)
    if not in_p.any():
        raise ValueError(f"{run_dir.name}: no Felix window inside the pulse")
    n_p = int(in_p.sum())
    ref_m = (ts + W <= tp0 - 0.02) & (ts >= tp0 - 0.02 - n_p * W)
    own = decay_summary(bins[ref_m].sum(0), bin_s, ref_m.sum() * W)
    pk = own["peak_bin"] if own["significant"] else None
    ref = decay_summary(bins[ref_m].sum(0), bin_s, ref_m.sum() * W, peak_bin=pk, dead_time_s=0.0)
    pulse = decay_summary(bins[in_p].sum(0), bin_s, n_p * W, peak_bin=pk, dead_time_s=0.0)
    show_m = (ts >= tp0 - before_s) & (ts < tp1 + after_s)
    tot = bins[show_m].sum(1) / W
    lc = np.full(tot.size, np.nan)
    if pk is not None:
        hr = fl._rebin_matrix(bins[show_m], REBIN)
        bw = bin_s * REBIN * 1e9
        p0 = int(round(PEAK_AT_NS / bw))
        hr = np.roll(hr, p0 - pk, axis=1)
        t = (np.arange(hr.shape[1]) + 0.5) * bw
        pre = (t > 0.3) & (t < PEAK_AT_NS - 0.6)
        rise = p0 - int(round(RISE_NS / bw))
        lc = (hr[:, rise:] - hr[:, pre].mean(1)[:, None]).sum(1) / W
    rows = list(csv.DictReader(open(run_dir / "ambit.csv")))
    amb = {k: np.array([float(r[k]) for r in rows if r["line"] == "1"]) for k in ("t_s", "s_630", "r_630")}
    return dict(run_dir=str(run_dir), label=meta.get("label", run_dir.name), t_pulse=(tp0, tp1), window_s=W,
                t=ts[show_m] + W / 2, total_cps=tot, laser_cps=lc, in_pulse=in_p[show_m], pulse=pulse, ref=ref,
                ambit=amb)


def plot_pulse_spad(run_dir, before_s: float = 1.0, after_s: float = 1.5, save_to=None, show: bool = True):
    """What the SPAD counts during the saturating pulse of a dose run (``pulse_spad``).

    Top left: total and laser-correlated rate per 20 ms window, with the SPAD's rated maximum.
    Top right: inside the pulse, the SPAD rate next to the AMBIT fluorescence (s_630) and LED
    reference (r_630), each relative to the last 0.3 s of the pulse: a linear SPAD would rise
    with s_630.
    Bottom: micro-time histograms of the pulse and of as long a stretch just before it,
    raw (left) and floor-subtracted (right), in counts per second per 47 ps bin.
    """
    import matplotlib.pyplot as plt

    p = pulse_spad(run_dir, before_s, after_s)
    tp0, tp1 = p["t_pulse"]
    P, R = p["pulse"], p["ref"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    ax = axes[0, 0]
    ax.semilogy(p["t"] - tp0, p["total_cps"], ".-", ms=3, lw=0.6, color="k", label="SPAD total")
    lc = p["laser_cps"]
    pos = lc > 0
    ax.semilogy((p["t"] - tp0)[pos], lc[pos], ".", ms=3, color="tab:red", label="laser-correlated (where > 0)")
    ax.axhline(SPAD_MAX_CPS, color="tab:orange", ls="--", lw=0.9, label=f"SPAD rated max {SPAD_MAX_CPS / 1e6:.1f} Mcps")
    ax.axvspan(0, tp1 - tp0, color="tab:orange", alpha=0.15, lw=0)
    ax.set_xlabel("time since the pulse started [s]")
    ax.set_ylabel(f"counts/s per {p['window_s'] * 1e3:.0f} ms window")
    ax.set_title(f"SPAD around the saturating pulse: {p['label']}")
    ax.legend(fontsize=8, loc="center right")
    ax = axes[0, 1]
    a = p["ambit"]
    end = a["t_s"] >= tp1 - 0.3
    ins = p["in_pulse"]
    tf = p["t"][ins]
    tot_in = p["total_cps"][ins]
    ax.plot(tf - tp0, tot_in / tot_in[tf >= tp1 - 0.3].mean(), "o-", ms=3, color="k",
            label=f"SPAD total ({tot_in.mean() / 1e6:.1f} Mcps)")
    ax.plot(a["t_s"] - tp0, a["s_630"] / a["s_630"][end].mean(), "s-", ms=2.5, lw=0.8, color="tab:green",
            label="AMBIT fluorescence s_630")
    ax.plot(a["t_s"] - tp0, a["r_630"] / a["r_630"][end].mean(), "^-", ms=2.5, lw=0.8, color="tab:blue",
            label="AMBIT LED reference r_630")
    ax.axhline(1, color="0.5", lw=0.5)
    ax.set_xlabel("time since the pulse started [s]")
    ax.set_ylabel("relative to the last 0.3 s of the pulse")
    ax.set_title("inside the pulse: a linear SPAD would follow the fluorescence")
    ax.legend(fontsize=8, loc="lower right")
    t = (np.arange(P["h"].size) + 0.5) * P["bin_ns"]
    for name, s, c in ((f"pulse ({P['seconds']:.2f} s)", P, "tab:orange"),
                       (f"just before ({R['seconds']:.2f} s, laser on)", R, "tab:red")):
        axes[1, 0].semilogy(t, s["h"] / s["seconds"], color=c, lw=0.8,
                            label=f"{name}: {s['total_cps'] / 1e6:.2f} Mcps total")
        y = (s["h"] - s["floor_per_bin"]) / s["seconds"]
        axes[1, 1].plot(t, y / 1e3, color=c, lw=0.8,
                        label=f"{name}: laser-correlated {s['laser_cps'] / 1e3:.0f} +- {s['laser_cps_err'] / 1e3:.0f} kcps")
    sig = np.sqrt(np.maximum(P["h"], 1.0)) / P["seconds"] / 1e3
    axes[1, 1].fill_between(t, -sig, sig, color="tab:orange", alpha=0.2, lw=0, label="pulse: +-1 sigma per bin")
    axes[1, 0].set_ylabel("counts/s per 47 ps bin")
    axes[1, 0].set_title("micro-time histograms, raw")
    axes[1, 1].axhline(0, color="0.5", lw=0.5)
    axes[1, 1].set_ylabel("kcounts/s per 47 ps bin, floor subtracted")
    axes[1, 1].set_title("laser-correlated part")
    for ax in axes[1]:
        ax.set_xlim(0, 12.5)
        ax.set_xlabel("micro-time [ns] (laser peak at 2.5 ns)")
        ax.legend(fontsize=8, loc="upper right")
    for x in axes.flat:
        x.grid(alpha=0.3)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=120)
    if show:
        plt.show()
    return fig
