"""Actinic pauses: what the SPAD counts when the AMBIT light stops (notebook ``4_ambit_gate_2.ipynb``).

The AMBIT runs ordinary arrun lines at 50 Hz. Per actinic level: 1 s of light, a short dark
pause, 1 s of light, a longer pause, 1 s of light, the longest pause. No chop firmware: a pause
is a dark line of a few samples. The Felix records contiguous hardware histograms (micro-times
from the laser SYNC) of 5 ms for the whole run.

Questions
    1. After LED-off, how fast does the SPAD come back to the laser-only floor (level 0)?
    2. Is the laser decay in a pause usable, and does it differ from the level-0 decay?
    3. What does the SPAD count while the LED is on?

Timing
    Every arrun line carries ~9 ms of fixed time (set-up before its first sample, read-out, LED
    switching), so a pause line of n samples at f Hz is dark for about (n - 1) / f + 9 ms: 29, 49
    and 89 ms for 2, 3 and 5 samples at 50 Hz. ``ag.time_axis`` ignores those overheads (fine for
    3 lines, ~0.1 s off for 13), so every lit line is placed on the Felix clock from its own
    measured LED edges, and every pause is the measured dark time between two lit lines.
    The AMBIT takes at most 16 lines (``ag.MAX_SEGMENTS``): an actinic level 0 has no LED to
    pause, so it is one dark line as long as a lit block.

All rates here are raw counts per second (no dead-time correction): in the lit windows the SPAD
is far outside the range of the dead-time model, and in the dark ones the correction is ~4 %.

Light series (``run_light_series``): the same measurement after a pre-illumination at several
levels. A one-sample arrun with the persist flag leaves the LED on while the AMBIT is idle (no
measuring flashes), the PC waits, switches the LED off and runs the train. ``analyse_light_series``
gives per level the AMBIT Fo, Fm and phi2, the lifetimes of the reference line (Fo) and of every
pause (Fm), and phi2 from the lifetimes; ``plot_light_series`` draws them.
"""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np

import ambit_chop as ac
import ambit_gate as ag
import blank_dose as bd
import fvfm_lifetime as fl
import reconv as rc

LEVELS = (0, 100, 250)        # actinic bytes, in run order
BLOCK_S = 1.0                 # light between two pauses
PAUSE_SAMPLES = (2, 3, 5)     # AMBIT samples per pause line: ~29, 49, 89 ms dark at 50 Hz
FREQ = 50                     # AMBIT measuring rate of every line
WINDOW_S = 0.005              # Felix histogram window
LED_ON_DELAY_S = 0.001        # a lit line switches its LED on ~1 ms after the line starts
LEVEL_COLOURS = {0: "0.4", 100: "tab:orange", 250: "tab:red"}


def _colour(level: int) -> str:
    return LEVEL_COLOURS.get(level, "tab:purple")


# ---------------------------------------------------------------------------
# Protocol and its predicted timing
# ---------------------------------------------------------------------------
def pause_protocol(levels=LEVELS, block_s: float = BLOCK_S, pause_samples=PAUSE_SAMPLES,
                   freq: int = FREQ) -> list[dict]:
    """arrun lines: per lit level, (block_s of light, a pause of n samples) for every n in
    ``pause_samples``; a level <= 3 (LED off) is one dark line of the same number of samples."""
    n_block = int(round(block_s * freq))
    prot = []
    for x in levels:
        if x <= ag.ACTINIC_ON:
            prot.append(dict(ac.line(1.0, freq, 0), pulses=len(pause_samples) * n_block + int(sum(pause_samples))))
        else:
            for n in pause_samples:
                prot.append(ac.line(block_s, freq, int(x)))
                prot.append(dict(ac.line(1.0, freq, 0), pulses=int(n)))
    if len(prot) > ag.MAX_SEGMENTS:
        raise ValueError(f"{len(prot)} lines, the AMBIT takes at most {ag.MAX_SEGMENTS}: fewer lit levels or pauses")
    return prot


def train_protocol(level: int, dark_points: int, block_s: float, pause_samples, freq: int = FREQ) -> list[dict]:
    """One actinic level as a train of lit blocks and short pauses, inside dark:

        dark (dark_points) | lit block_s, pause n1 | lit block_s, pause n2 | ... | dark (dark_points)

    One (lit, pause) pair per entry of ``pause_samples`` (AMBIT points in that pause; 1 point is
    the shortest pause, ~9.6 ms dark at 100 Hz). The last pause runs straight into the dark tail.
    The first dark line is the reference of ``analyse_pauses``. 2 + 2 x len(pause_samples) lines."""
    if int(level) <= ag.ACTINIC_ON:
        raise ValueError(f"level {level}: the LED switches on only above {ag.ACTINIC_ON}")
    if dark_points < 1 or any(int(n) < 1 for n in pause_samples):
        raise ValueError("dark_points and every pause need at least 1 AMBIT point")
    dark = lambda n: dict(ac.line(1.0, freq, 0), pulses=int(n))
    prot = [dark(dark_points)]
    for n in pause_samples:
        prot.append(ac.line(block_s, freq, int(level)))
        prot.append(dark(n))
    prot.append(dark(dark_points))
    if len(prot) > ag.MAX_SEGMENTS:
        raise ValueError(f"{len(prot)} lines, the AMBIT takes at most {ag.MAX_SEGMENTS}: at most "
                         f"{(ag.MAX_SEGMENTS - 2) // 2} pauses")
    n_pts = sum(s["pulses"] for s in prot)
    if n_pts > ag.MAX_POINTS:
        raise ValueError(f"{n_pts} AMBIT points, at most {ag.MAX_POINTS} per run: fewer dark points or shorter blocks")
    return prot


def line_roles(protocol: list[dict]) -> list[dict]:
    """Per line: role 'lit', 'pause' (a dark line right after a lit one) or 'dark'; its level;
    for a pause, its index within the level (0 = first pause of that level)."""
    out, count = [], {}
    for i, s in enumerate(protocol):
        if s["actinic"] > ag.ACTINIC_ON:
            out.append(dict(role="lit", level=int(s["actinic"]), pause=None))
        elif i > 0 and protocol[i - 1]["actinic"] > ag.ACTINIC_ON:
            lv = int(protocol[i - 1]["actinic"])
            j = count.get(lv, 0)
            count[lv] = j + 1
            out.append(dict(role="pause", level=lv, pause=j))
        else:
            out.append(dict(role="dark", level=0, pause=None))
    return out


def predict(protocol: list[dict]) -> list[dict]:
    """Predicted line spans (s from the start of line 0, firmware timing of ``ambit_chop``) with
    the LED-on span of lit lines and the dark span of pauses (``open``: no lit line after it)."""
    spans = ac.nominal_spans(protocol)
    roles = line_roles(protocol)
    rows = []
    for i, (s, sp, r) in enumerate(zip(protocol, spans, roles)):
        rows.append(dict(line=i, actinic=int(s["actinic"]), freq=int(s["freq"]), pulses=int(s["pulses"]),
                         t0=sp["t0"], t1=sp["t1"], t_edge0=sp["t_edge0"], **r))
    for i, r in enumerate(rows):
        if r["role"] == "lit":
            r["on"], r["off"] = r["t0"] + LED_ON_DELAY_S, r["t1"]
    for i, r in enumerate(rows):
        if r["role"] != "pause":
            continue
        a = rows[i - 1]["off"]
        nxt = rows[i + 1] if i + 1 < len(rows) else None
        closed = nxt is not None and nxt["role"] == "lit"
        b = nxt["on"] if closed else r["t1"] + 0.0003 + LED_ON_DELAY_S
        r.update(dark0=a, dark1=b, open=not closed)
    return rows


def window_check(protocol: list[dict], window_s: float) -> str | None:
    """Why a Felix window of ``window_s`` cannot resolve the pauses of ``protocol``, or None.

    A pause only yields a decay when a window lies fully inside its dark span whatever the window
    phase, and the LED-edge alignment needs one fully dark window per pause: both hold when
    window <= dark / 2. Seen on 2026-11-05: 1-point pauses (9.3 ms) with 10 ms windows gave
    "0 Felix windows" for every pause and lines placed from the nominal latency.
    """
    pauses = [r for r in predict(protocol) if r["role"] == "pause"]
    if not pauses:
        return None
    dark = min(r["dark1"] - r["dark0"] for r in pauses)
    if window_s <= dark / 2:
        return None
    r = min(pauses, key=lambda p: p["dark1"] - p["dark0"])
    overhead = dark - r["pulses"] / r["freq"]              # dark time around the pause's own points
    need = max(r["pulses"] + 1, int(np.ceil((2 * window_s - overhead) * r["freq"])))
    return (f"The Felix window ({window_s * 1e3:.0f} ms) is longer than half of the shortest pause "
            f"({dark * 1e3:.1f} ms dark): no window fits inside it reliably and the LED edges blur. "
            f"Use a window of at most {max(1, int(dark * 500))} ms, or pauses of at least {need} points.")


def describe(protocol: list[dict]) -> str:
    pred = predict(protocol)
    out = [f"{len(protocol)} lines (max {ag.MAX_SEGMENTS}), {sum(s['pulses'] for s in protocol)} AMBIT samples, "
           f"about {pred[-1]['t1']:.2f} s",
           "line  role    actinic  samples  freq   starts s   length ms   dark ms (pauses, predicted)"]
    for r in pred:
        dark = ""
        if r["role"] == "pause":
            dark = f"{(r['dark1'] - r['dark0']) * 1e3:6.1f}" + ("  (no light after it: stays dark)" if r["open"] else "")
        out.append(f"{r['line']:>4}  {r['role']:<6}  {r['actinic']:>7}  {r['pulses']:>7}  {r['freq']:>4}  "
                   f"{r['t0']:>9.3f}  {(r['t1'] - r['t0']) * 1e3:>10.1f}   {dark}")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------
def led_intervals(counts: np.ndarray, starts: np.ndarray, W: float, dark: float) -> list[dict]:
    """LED-on intervals in the Felix totals (counts per window): runs of windows well above the
    dark level; the edges are placed inside their windows from the fill fraction of the partly
    lit neighbours (the lit level is the median of the run's inner windows)."""
    thr = dark + max(10.0 * np.sqrt(max(dark, 1.0)), dark)
    st, en = ac._runs(counts > thr)
    out = []
    for a, b in zip(st, en):                       # windows a .. b-1 lit
        inner = counts[a + 1:b - 1] if b - a >= 3 else counts[a:b]
        level = float(np.median(inner))
        frac = lambda c: float(np.clip((c - dark) / max(level - dark, 1.0), 0.0, 1.0))
        before = frac(counts[a - 1]) if a > 0 else 0.0
        after = frac(counts[b]) if b < counts.size else 0.0
        on = starts[a] + (1.0 - frac(counts[a])) * W - before * W
        off = starts[b - 1] + frac(counts[b - 1]) * W + after * W
        out.append(dict(on=float(on), off=float(off), level_counts=level, n_windows=int(b - a)))
    return out


def place(pred: list[dict], ivs: list[dict], t_send: float, log=print) -> tuple[float, dict]:
    """Run start (Felix monotonic time of line 0's start) from the measured LED onsets matched
    in order to the lit lines; measured on/off per lit line and dark span per pause."""
    lit = [r for r in pred if r["role"] == "lit"]
    info = dict(n_intervals=len(ivs), n_lit=len(lit), mode="nominal", residuals_ms=[])
    if not lit or len(ivs) != len(lit):
        log(f"WARNING: {len(ivs)} LED-on intervals found for {len(lit)} lit lines: "
            "positions from the nominal latency")
        return t_send + ag.NOMINAL_LATENCY_S, info
    offs = np.array([iv["on"] - r["on"] for iv, r in zip(ivs, lit)])
    run0 = float(np.median(offs))
    info.update(mode="led-edges", latency_s=run0 - t_send,
                residuals_ms=[float(x) for x in (offs - run0) * 1e3],
                lit_len_error_ms=[float(((iv["off"] - iv["on"]) - (r["off"] - r["on"])) * 1e3)
                                  for iv, r in zip(ivs, lit)])
    for iv, r in zip(ivs, lit):
        r["on_meas"], r["off_meas"] = iv["on"] - run0, iv["off"] - run0
        r["level_cps"] = iv["level_counts"]
    for i, r in enumerate(pred):
        if r["role"] != "pause":
            continue
        a = pred[i - 1]["off_meas"]
        b = pred[i + 1]["on_meas"] if not r["open"] else a + (r["dark1"] - r["dark0"])
        r["dark0_meas"], r["dark1_meas"] = a, b
    return run0, info


def run_pauses(ambit: ag.AmbitLink, felix: ag.FelixLink, out_dir, protocol: list[dict],
               window_s: float = WINDOW_S, pre_s: float = 1.0, tail_s: float = 1.0, note: str = "",
               threshold_v: float = ag.THRESHOLD_V, laser_hz: float = ag.LASER_REP_RATE_HZ, log=print,
               extra: dict | None = None) -> Path:
    """One run: the Felix records ``window_s`` histograms from ``pre_s`` before the arrun to
    ``tail_s`` after it; the lines are placed from the LED edges. Saves felix.npz, ambit.csv, run.json
    (``extra``: more keys for run.json, e.g. the pre-illumination of a light series)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ag.build_arrun(protocol)
    pred = predict(protocol)
    duration = pred[-1]["t1"]
    try:
        ac.disarm_chop(ambit)        # no gaps: a chop armed and never consumed would apply here
    except Exception:
        pass
    felix.set_thresholds(threshold_v)
    rates = felix.count_rates(0.5)
    log("Felix input rates: " + ", ".join(f"ch{c}={r:,.0f}/s" for c, r in rates.items()))
    if rates.get(ag.SYNC_CH, 0) < 0.9 * laser_hz:
        log("WARNING: SYNC rate below the laser repetition rate: check the laser and the SYNC cable")
    warn = ag.rejection_warning(felix.tdc_rejection(0.5), rates.get(ag.SPAD_CH, 0))
    if warn:
        raise RuntimeError(warn)         # the run would only produce empty histograms
    warn = window_check(protocol, window_s)
    if warn:
        log("WARNING: " + warn)
    felix.configure_histogram(ag.SYNC_CH, ag.SPAD_CH, 0.0, 1.0 / laser_hz)
    n_win = int(np.ceil((pre_s + ag.NOMINAL_LATENCY_S + duration + tail_s + 0.5) / window_s))
    W = felix.start_windows(window_s, n_win)
    log(f"Felix: {n_win} windows of {W * 1e3:.2f} ms")
    time.sleep(pre_s)
    log("AMBIT: " + cmd.strip())
    data, t_send, t_done = ag.run_arrun(ambit, cmd, protocol, log=log)
    if not felix.wait_done(timeout_s=tail_s + 10.0):
        log(f"WARNING: the Felix delivered {len(felix.windows)} of {n_win} windows")
    wins = list(felix.windows)
    if not wins:
        raise RuntimeError("no Felix windows")
    bins = np.vstack([w.bins for w in wins]).astype(np.uint32)
    k = np.array([w.k for w in wins])
    win_id = np.array([w.win_id for w in wins])
    lost = int(sum(bool(w.lost) for w in wins))
    starts = felix.t_soa + k * W
    if lost or (np.ptp(win_id) > 0 and np.any(np.diff(win_id) != 1)):
        log(f"WARNING: the Felix lost windows (flag {lost}, id gaps {int(np.sum(np.diff(win_id) != 1))}): "
            "window times after a gap are wrong; try 10 ms windows")
    counts = bins.sum(1).astype(float)
    pre = starts + W <= t_send
    dark = float(np.median(counts[pre])) if pre.sum() >= 10 else float(np.percentile(counts, 10))
    ivs = led_intervals(counts, starts, W, dark)
    run0, info = place(pred, ivs, t_send, log=log)
    if info["mode"] == "led-edges":
        log(f"lines placed from {len(ivs)} LED edges: run start = send + {info['latency_s'] * 1e3:.0f} ms, "
            f"onset residuals {min(info['residuals_ms']):+.2f}..{max(info['residuals_ms']):+.2f} ms")
    # AMBIT sample times: firmware timing, each lit line shifted onto its own measured onset
    s = np.asarray(data.get("s_630", np.zeros(0)), float)
    r = np.asarray(data.get("r_630", np.zeros(0)), float)
    rows, i = [], 0
    for p, seg in zip(pred, protocol):
        shift = (p["on_meas"] - p["on"]) if "on_meas" in p else 0.0
        for kk in range(seg["pulses"]):
            if i >= min(s.size, r.size):
                break
            rows.append([p["line"], kk, f"{p['t_edge0'] + shift + kk / seg['freq']:.5f}", int(s[i]), int(r[i]),
                         f"{s[i] / r[i]:.6f}" if r[i] > 0 else "nan", seg["actinic"], seg["freq"], p["role"]])
            i += 1
    with open(out_dir / "ambit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["line", "k", "t_s", "s_630", "r_630", "ratio", "actinic", "freq", "role"])
        wr.writerows(rows)
    np.savez_compressed(out_dir / "felix.npz", bins=bins, t_start_s=starts - run0, window_s=W,
                        bin_s=felix.bin_s, left_s=felix.left_s, win_id=win_id)
    meta = dict(kind="pauses", protocol=protocol, cmd=cmd.strip(), lines=pred, align=info, window_s=W,
                bin_s=felix.bin_s, laser_hz=laser_hz, dark_counts_per_window=dark, t_send_s=t_send - run0,
                rates_before={str(c): v for c, v in rates.items()}, lost=lost,
                ambit=dict(name=ambit.name, fw=ambit.fw), felix=felix.description, note=note,
                saved=time.strftime("%Y-%m-%d %H:%M:%S"), **(extra or {}))
    (out_dir / "run.json").write_text(json.dumps(meta, indent=1, default=float), encoding="utf-8")
    log(f"saved to {out_dir}")
    return out_dir


def run_levels(ambit: ag.AmbitLink, felix: ag.FelixLink, session, label: str, levels,
               block_s: float = BLOCK_S, pause_samples=PAUSE_SAMPLES, freq: int = FREQ, rest_s: float = 0.0,
               **kw) -> list[Path]:
    """Many actinic levels, one arrun per lit level. The AMBIT takes 16 lines per arrun, so the
    levels are not chained in one command: every run is the actinic-0 reference line followed by
    one level's lit blocks and pauses (7 lines with 3 pauses). The pauses stay inside their arrun,
    so their dark time stays firmware-timed; between two runs the leaf is under the laser alone
    (Felix tail + pre time, data transfer, ``rest_s``), then gets its own reference line.

    ``levels``: actinic bytes in run order; levels <= 3 (LED off) are skipped, since every run
    starts with the reference. ``kw`` goes to ``run_pauses`` (window_s, note, threshold_v,
    laser_hz, log). Each run is saved in ``session``/pauses_<label>_a<level>_<time>. Returns the
    run folders."""
    log = kw.get("log", print)
    lit = [int(x) for x in levels if x > ag.ACTINIC_ON]
    if not lit:
        raise ValueError("no lit level (> 3) in levels")
    prots = [pause_protocol((0, lv), block_s, pause_samples, freq) for lv in lit]   # all checked before starting
    dur = predict(prots[0])[-1]["t1"]
    log(f"{len(lit)} runs of {len(prots[0])} lines, {dur:.1f} s of AMBIT each, plus Felix pre/tail and transfer")
    dirs = []
    for i, (lv, prot) in enumerate(zip(lit, prots)):
        log(f"--- actinic {lv} ({i + 1} of {len(lit)})")
        dirs.append(run_pauses(ambit, felix, Path(session) / f"pauses_{label}_a{lv:03d}_{time.strftime('%H%M%S')}",
                               prot, **kw))
        if rest_s > 0 and i < len(lit) - 1:
            time.sleep(rest_s)
    return dirs


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
def window_stats(bins: np.ndarray, bin_s: float, peak_bin: int, W: float) -> dict:
    """Per window, rotated like ``blank_dose.decay_summary``: total, laser-correlated and floor
    rates (counts/s, raw), mean arrival after the peak (only where the laser part is > 5 sigma)."""
    hr = fl._rebin_matrix(np.asarray(bins, float), bd.REBIN)
    bw = bin_s * bd.REBIN * 1e9
    p0 = int(round(bd.PEAK_AT_NS / bw))
    hr = np.roll(hr, p0 - peak_bin, axis=1)
    t = (np.arange(hr.shape[1]) + 0.5) * bw
    pre = (t > 0.3) & (t < bd.PEAK_AT_NS - 0.6)
    floor = hr[:, pre].mean(1)
    floor_err = hr[:, pre].std(1, ddof=1) / np.sqrt(pre.sum())
    rise = p0 - int(round(bd.RISE_NS / bw))
    tail = hr[:, rise:]
    lc = (tail - floor[:, None]).sum(1)
    lc_err = np.sqrt(np.maximum(tail.sum(1), 1.0) + (floor_err * tail.shape[1]) ** 2)
    sel = (t >= t[p0]) & (t < t[p0] + 6.0)
    dt = t[sel] - t[p0]
    w = hr[:, sel] - floor[:, None]
    ws = w.sum(1)
    ok = (lc > 5 * lc_err) & (ws > 0)
    wsn = np.where(ws > 0, ws, 1.0)
    tm = (w * dt[None, :]).sum(1) / wsn
    d = dt[None, :] - tm[:, None]
    var = (d ** 2 * hr[:, sel]).sum(1) + (floor_err * d.sum(1)) ** 2
    return dict(total_cps=np.asarray(bins).sum(1) / W, laser_cps=lc / W, laser_cps_err=lc_err / W,
                floor_cps=floor * hr.shape[1] / W, t_mean_ns=np.where(ok, tm, np.nan),
                t_mean_err_ns=np.where(ok, np.sqrt(var) / wsn, np.nan), significant=ok)


def load_run(run_dir) -> dict:
    run_dir = Path(run_dir)
    meta = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    z = np.load(run_dir / "felix.npz")
    rows = list(csv.DictReader(open(run_dir / "ambit.csv")))
    amb = dict(t=np.array([float(r["t_s"]) for r in rows]), ratio=np.array([float(r["ratio"]) for r in rows]),
               line=np.array([int(r["line"]) for r in rows]), actinic=np.array([int(r["actinic"]) for r in rows]),
               s_630=np.array([float(r["s_630"]) for r in rows]), r_630=np.array([float(r["r_630"]) for r in rows]))
    return dict(run_dir=run_dir, meta=meta, bins=z["bins"], ts=z["t_start_s"].astype(float),
                W=float(z["window_s"]), bin_s=float(z["bin_s"]), ambit=amb)


def _summary(d: dict, m: np.ndarray, pk: int) -> dict | None:
    """decay_summary of the windows ``m`` (raw rates), plus the IRF reconvolution fit ``rc``."""
    if not m.any():
        return None
    h = d["bins"][m].sum(0)
    s = bd.decay_summary(h, d["bin_s"], m.sum() * d["W"], peak_bin=pk, dead_time_s=0.0)
    s["n_windows"] = int(m.sum())
    s["rc"] = rc.fit(h, d["bin_s"], pk, d["irf_p"]) if d.get("irf_p") is not None and s["significant"] else None
    return s


def analyse_pauses(run_dir, irf: dict | None = None) -> dict:
    """Reference (laser only), lit levels and every pause of a saved run. ``irf``: the measured IRF
    for the reconvolution fits (``blank_dose.load_irf``; None: the first one in irf-library)."""
    d = load_run(run_dir)
    irf = irf if irf is not None else bd.load_irf()
    d["irf_p"] = rc.prepare_irf(irf, d["bins"].shape[1], d["bin_s"]) if irf is not None else None
    meta, ts, W = d["meta"], d["ts"], d["W"]
    lines = meta["lines"]
    warnings = []
    if meta["align"]["mode"] != "led-edges":
        warnings.append("lines placed from the nominal latency, not from the LED edges: pause positions uncertain")
    # reference: the level-0 dark line, else the windows before the run
    dark_lines = [r for r in lines if r["role"] == "dark"]
    if dark_lines:
        r0 = dark_lines[0]
        ref_m = (ts >= r0["t0"] + 0.02) & (ts + W <= r0["t1"] - 0.02)
        ref_name = f"actinic 0 line ({r0['t1'] - r0['t0']:.2f} s)"
    else:
        ref_m = (ts + W <= -0.02)
        ref_name = "before the run"
    own = bd.decay_summary(d["bins"][ref_m].sum(0), d["bin_s"], ref_m.sum() * W, dead_time_s=0.0)
    if not own["significant"]:
        raise ValueError("no laser peak in the reference windows: was the laser on?")
    pk = own["peak_bin"]
    ref = _summary(d, ref_m, pk)
    ref["name"] = ref_name
    ws = window_stats(d["bins"], d["bin_s"], pk, W)
    before = _summary(d, ts + W <= -0.02, pk)
    # lit levels: windows fully inside the measured LED-on spans
    lit = {}
    for r in lines:
        if r["role"] != "lit":
            continue
        a, b = r.get("on_meas", r["on"]), r.get("off_meas", r["off"])
        m = (ts >= a) & (ts + W <= b)
        lit.setdefault(r["level"], []).append(m)
    lit_sum = {}
    for lv, ms in lit.items():
        m = np.any(ms, axis=0)
        tot = ws["total_cps"][m]
        lit_sum[lv] = dict(level=lv, n_windows=int(m.sum()), total_cps=float(np.median(tot)),
                           total_min=float(tot.min()), total_max=float(tot.max()))
        if lit_sum[lv]["total_cps"] > bd.SPAD_MAX_CPS:
            warnings.append(f"actinic {lv}: SPAD at {lit_sum[lv]['total_cps'] / 1e6:.1f} Mcps while lit, "
                            f"above its rated {bd.SPAD_MAX_CPS / 1e6:.1f} Mcps")
    # pauses
    pauses = []
    for r in lines:
        if r["role"] != "pause":
            continue
        a, b = r.get("dark0_meas", r["dark0"]), r.get("dark1_meas", r["dark1"])
        m = (ts >= a) & (ts + W <= b)
        idx = np.flatnonzero(m)
        p = dict(line=r["line"], level=r["level"], pause=r["pause"], samples=r["pulses"], dark0=a, dark1=b,
                 dark_ms=(b - a) * 1e3, dark_pred_ms=(r["dark1"] - r["dark0"]) * 1e3, open=r["open"],
                 n_windows=int(m.sum()), summary=_summary(d, m, pk))
        if idx.size:
            i0 = idx[0]
            p["first"] = dict(t_after_off_ms=(ts[i0] - a) * 1e3, total_cps=float(ws["total_cps"][i0]),
                              laser_cps=float(ws["laser_cps"][i0]), floor_cps=float(ws["floor_cps"][i0]),
                              t_mean_ns=float(ws["t_mean_ns"][i0]), t_mean_err_ns=float(ws["t_mean_err_ns"][i0]))
        else:
            warnings.append(f"actinic {r['level']} pause {r['pause']}: no Felix window fully inside "
                            f"its {(b - a) * 1e3:.0f} ms")
        pauses.append(p)
    # after the run: the LED stays off after the last line
    end = max(r["t1"] for r in lines)
    tail = _summary(d, (ts >= end + 0.1) & (ts + W <= end + 1.0), pk)
    return dict(run_dir=str(d["run_dir"]), note=meta.get("note", ""), W=W, ts=ts, ws=ws, lines=lines,
                align=meta["align"], ref=ref, before=before, lit=lit_sum, pauses=pauses, tail=tail,
                ambit=d["ambit"], peak_bin=pk, bins=d["bins"], bin_s=d["bin_s"], warnings=warnings)


def report(res: dict) -> str:
    k = lambda x: f"{x / 1e3:7.1f}"
    al = res["align"]
    out = [f"{Path(res['run_dir']).name}" + (f"  ({res['note']})" if res["note"] else "")]
    if al["mode"] == "led-edges":
        out.append(f"  lines placed from {al['n_intervals']} LED edges; onsets within "
                   f"{min(al['residuals_ms']):+.2f}..{max(al['residuals_ms']):+.2f} ms of the firmware timing, "
                   f"lit lengths {min(al['lit_len_error_ms']):+.2f}..{max(al['lit_len_error_ms']):+.2f} ms")
    out += [f"  WARNING: {w}" for w in res["warnings"]]
    f = res["ref"]
    out.append(f"  reference, {f['name']}: SPAD {f['total_cps'] / 1e3:.1f} kcps total, laser-correlated "
               f"{f['laser_cps'] / 1e3:.1f} kcps, floor {f['floor_cps'] / 1e3:.1f} kcps, t_mean "
               f"{f['t_mean_ns']:.3f} +- {f['t_mean_err_ns']:.3f} ns, tau1 {f['tau1_ns']:.3f} ns")
    for name in ("before", "tail"):
        s = res[name]
        if s is not None and s["significant"]:
            lab = "before the run" if name == "before" else "0.1-1 s after the run"
            out.append(f"  {lab:<22}: SPAD {s['total_cps'] / 1e3:.1f} kcps, laser {s['laser_cps'] / 1e3:.1f}, "
                       f"floor {s['floor_cps'] / 1e3:.1f} kcps, t_mean {s['t_mean_ns']:.3f} +- {s['t_mean_err_ns']:.3f} ns")
    for lv, s in res["lit"].items():
        out.append(f"  actinic {lv:>3} lit: SPAD {s['total_cps'] / 1e6:.2f} Mcps (median of {s['n_windows']} "
                   f"windows, {s['total_min'] / 1e6:.2f}..{s['total_max'] / 1e6:.2f})")
    out.append("  pauses (kcps, raw):        dark ms  windows | first window: after off  total  laser  floor  "
               "t_mean | whole pause: total  laser  floor  t_mean +- err   tau1")
    for p in res["pauses"]:
        fw, s = p.get("first"), p["summary"]
        a = (f"{fw['t_after_off_ms']:6.1f} ms {k(fw['total_cps'])} {k(fw['laser_cps'])} {k(fw['floor_cps'])}  "
             f"{fw['t_mean_ns']:5.2f}" if fw else " " * 40)
        b = (f"{k(s['total_cps'])} {k(s['laser_cps'])} {k(s['floor_cps'])}  {s['t_mean_ns']:.3f} +- "
             f"{s['t_mean_err_ns']:.3f}  {s['tau1_ns']:.3f}" if s is not None and s["significant"] else "-")
        out.append(f"  actinic {p['level']:>3} pause {p['pause'] + 1} ({p['samples']} samples) {p['dark_ms']:6.1f} "
                   f"{p['n_windows']:>4}    | {a} | {b}" + ("  (stays dark)" if p["open"] else ""))
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def _lit_mask(res: dict) -> np.ndarray:
    """Windows that overlap any LED-on span (their laser-correlated part is meaningless)."""
    ts, W = res["ts"], res["W"]
    m = np.zeros(ts.size, bool)
    for r in res["lines"]:
        if r["role"] == "lit":
            m |= (ts + W > r.get("on_meas", r["on"])) & (ts < r.get("off_meas", r["off"]))
    return m


def _shade(ax, res, alpha=0.18):
    for r in res["lines"]:
        if r["role"] == "lit":
            ax.axvspan(r.get("on_meas", r["on"]), r.get("off_meas", r["off"]), color=_colour(r["level"]),
                       alpha=alpha, lw=0)


def ambit_fo_fm(res: dict, n: int = 50) -> dict:
    """AMBIT fluorescence relative yield (s/r): Fo = mean of the last ``n`` points before the LED
    first comes on; Fm = mean of the last ``n`` points measured with the LED on, per lit level
    (pause points have actinic 0, so they are not among them). Each as dict(mean, sem, n, t0, t1),
    t0/t1 = time of the first and last point used."""
    a = res["ambit"]
    ok = np.isfinite(a["ratio"])
    lit = [r for r in res["lines"] if r["role"] == "lit"]

    def stat(m):
        idx = np.flatnonzero(m & ok)
        idx = idx[np.argsort(a["t"][idx])][-n:]
        if idx.size == 0:
            return None
        y = a["ratio"][idx]
        return dict(mean=float(y.mean()), sem=float(y.std(ddof=1) / np.sqrt(y.size)) if y.size > 1 else float("nan"),
                    n=int(y.size), t0=float(a["t"][idx].min()), t1=float(a["t"][idx].max()))

    t_on = min((r.get("on_meas", r["on"]) for r in lit), default=np.inf)
    out = dict(fo=stat((a["t"] < t_on) & (a["actinic"] <= ag.ACTINIC_ON)), fm={})
    for lv in sorted({int(x) for x in a["actinic"] if x > ag.ACTINIC_ON}):
        out["fm"][lv] = stat(a["actinic"] == lv)
    return out


def plot_overview(res: dict, save_to=None, show: bool = True, roll_n: int = 10, n_mean: int = 50):
    """The whole run: SPAD total per 5 ms window and the laser-correlated rate in the dark windows
    (top), and the AMBIT fluorescence relative yield (s/r) of every measuring pulse (bottom), with
    its rolling mean over ``roll_n`` consecutive points (all points, pauses included) and Fo / Fm
    from ``ambit_fo_fm`` (mean of the last ``n_mean`` points before the first LED-on / with the LED
    on, per level), given in the legend."""
    import matplotlib.pyplot as plt

    ts, W, ws = res["ts"], res["W"], res["ws"]
    tc = ts + W / 2
    lit_any = _lit_mask(res)
    fig, axes = plt.subplots(2, 1, figsize=(14, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[1.3, 0.8]))
    ax = axes[0]
    ax.semilogy(tc, ws["total_cps"], "-", lw=0.6, color="k", label="SPAD total")
    pos = ws["laser_cps"] > 0
    ax.semilogy(tc[pos & ~lit_any], ws["laser_cps"][pos & ~lit_any], ".", ms=2, color="tab:red",
                label="laser-correlated (dark windows)")
    ax.axhline(bd.SPAD_MAX_CPS, color="tab:orange", ls="--", lw=0.8, label="SPAD rated max")
    _shade(ax, res)
    ax.set_ylabel(f"counts/s per {W * 1e3:.3g} ms window")
    ax.set_title(f"{Path(res['run_dir']).name}: SPAD counts (shaded: LED on, colour = actinic level)")
    ax.legend(fontsize=8, loc="upper left")
    ax = axes[1]
    a = res["ambit"]
    for lv in np.unique(a["actinic"]):
        mm = a["actinic"] == lv
        ax.plot(a["t"][mm], a["ratio"][mm], ".", ms=3, color=_colour(int(lv)), label=f"actinic {lv}")
    ok = np.isfinite(a["ratio"])
    o = np.argsort(a["t"][ok])
    tt, yy = a["t"][ok][o], a["ratio"][ok][o]
    if yy.size >= roll_n:
        k = np.ones(roll_n) / roll_n               # time of each mean = mean time of its points
        ax.plot(np.convolve(tt, k, "valid"), np.convolve(yy, k, "valid"), color="k", lw=1.0,
                label=f"rolling mean, {roll_n} points")
    ff = ambit_fo_fm(res, n_mean)
    fo = ff["fo"]
    labels = [f"Fo = {fo['mean']:.4f} ± {fo['sem']:.4f}"] if fo else []
    for lv, st in ff["fm"].items():
        if not st:
            continue
        tag = f" (actinic {lv})" if len(ff["fm"]) > 1 else ""
        labels.append(f"Fm{tag} = {st['mean']:.4f} ± {st['sem']:.4f}")
        if fo and st["mean"] > 0:                  # phi2 = (Fm - Fo) / Fm = 1 - Fo/Fm
            phi = (st["mean"] - fo["mean"]) / st["mean"]
            err = fo["mean"] / st["mean"] * np.hypot(fo["sem"] / fo["mean"], st["sem"] / st["mean"])
            labels.append(r"$\phi_2$" + f"{tag} = (Fm - Fo) / Fm = {phi:.3f} ± {err:.3f}")
    for lab in labels:                             # legend text only, no marker
        ax.plot([], [], " ", label=lab)
    _shade(ax, res)
    ax.set_ylabel("Ambit (fluo. rel. yield)")
    ax.set_xlabel("time since the start of the AMBIT run [s]")
    ax.legend(fontsize=8, loc="upper left")
    for x in axes:
        x.grid(alpha=0.3)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=110)
    if show:
        plt.show()
    return fig


def plot_raw(run_dir, irf: dict | None = None, save: bool = True, show: bool = False):
    """The raw data of a saved run, before and independent of any fit: SPAD counts per Felix window
    with the lit blocks and the pauses where the alignment put them (top), the AMBIT trace (middle)
    and the summed raw decays of the reference line and of every pause (bottom), the IRF overlaid.
    Drawn by the GUI for every run, also when the analysis fails. Saved as raw.png in ``run_dir``."""
    import matplotlib.pyplot as plt

    d = load_run(run_dir)
    meta, bins, ts, W, bin_s, amb = d["meta"], d["bins"], d["ts"], d["W"], d["bin_s"], d["ambit"]
    per = bins.sum(1)
    tc = ts + W / 2
    lines = meta.get("lines", [])
    lit = [r for r in lines if r["role"] == "lit"]
    pauses = [r for r in lines if r["role"] == "pause"]
    darks = [r for r in lines if r["role"] == "dark"]
    pre = per[ts + W <= -0.02]
    dark_level = meta.get("dark_counts_per_window") or (float(np.median(pre)) if pre.size else 0.0)
    align = meta.get("align", {})
    rates = meta.get("rates_before", {})

    fig, axes = plt.subplots(3, 1, figsize=(14, 11.5), gridspec_kw=dict(height_ratios=[1.2, 0.7, 1.1]))
    # top: counts per window with the protocol as placed
    ax = axes[0]
    ax.semilogy(tc, np.maximum(per, 0.5), "-", lw=0.6, color="k", label=f"SPAD counts per {W * 1e3:.3g} ms window")
    for r in lit:
        ax.axvspan(r.get("on_meas", r["on"]), r.get("off_meas", r["off"]), color=_colour(int(r["level"])),
                   alpha=0.15, lw=0)
    inside = []
    for k, r in enumerate(pauses):
        m = (ts >= r["dark0"]) & (ts + W <= r["dark1"])
        inside.append(int(m.sum()))
        ax.axvspan(r["dark0"], r["dark1"], color="tab:blue", alpha=0.3, lw=0, label="pause (dark span)" if k == 0 else None)
        ax.annotate(f"{int(m.sum())} win", ((r["dark0"] + r["dark1"]) / 2, max(per.max(), 1)), ha="center",
                    va="bottom", fontsize=7, color="tab:blue")
    if darks:
        ax.axvspan(darks[0]["t0"], darks[0]["t1"], color="gray", alpha=0.12, lw=0, label="reference line (laser only)")
    if dark_level:
        thr = dark_level + max(10.0 * np.sqrt(max(dark_level, 1.0)), dark_level)
        ax.axhline(thr, color="tab:red", ls="--", lw=0.8, label="LED-edge threshold")
    ax.axhline(bd.SPAD_MAX_CPS * W, color="tab:orange", ls=":", lw=0.8, label="SPAD rated max")
    ax.set_ylabel("counts per window")
    ax.set_title(f"{Path(run_dir).name}: raw Felix windows (shaded: LED on; blue: pauses placed from "
                 f"{align.get('mode', '?')}, {align.get('n_intervals', '?')}/{align.get('n_lit', '?')} LED edges found)")
    ax.legend(fontsize=8, loc="upper left")
    # middle: AMBIT trace
    ax = axes[1]
    for lv in np.unique(amb["actinic"]):
        mm = amb["actinic"] == lv
        ax.plot(amb["t"][mm], amb["ratio"][mm], ".", ms=3, color=_colour(int(lv)), label=f"actinic {lv}")
    for r in lit:
        ax.axvspan(r.get("on_meas", r["on"]), r.get("off_meas", r["off"]), color=_colour(int(r["level"])), alpha=0.15, lw=0)
    for r in pauses:
        ax.axvspan(r["dark0"], r["dark1"], color="tab:blue", alpha=0.3, lw=0)
    ax.set_ylabel("AMBIT s/r (rel. yield)")
    ax.set_xlabel("time since the start of the AMBIT run [s]")
    ax.legend(fontsize=8, loc="upper left")
    # bottom: raw decays
    ax = axes[2]
    t_ns = (np.arange(bins.shape[1]) + 0.5) * bin_s * 1e9
    if darks:
        r0 = darks[0]
        ref_m = (ts >= r0["t0"] + 0.02) & (ts + W <= r0["t1"] - 0.02)
        ref_name = f"reference: actinic 0 line, {int(ref_m.sum())} windows"
    else:
        ref_m = ts + W <= -0.02
        ref_name = f"reference: before the run, {int(ref_m.sum())} windows"
    h_ref = bins[ref_m].sum(0) if ref_m.any() else np.zeros(bins.shape[1])
    ax.semilogy(t_ns, np.maximum(h_ref, 0.5), "-", lw=0.8, color="k", label=f"{ref_name}, {int(h_ref.sum()):,} counts")
    for k, r in enumerate(pauses):
        m = (ts >= r["dark0"]) & (ts + W <= r["dark1"])
        if m.any():
            h = bins[m].sum(0)
            ax.semilogy(t_ns, np.maximum(h, 0.5), "-", lw=0.7, alpha=0.8,
                        label=f"pause {k + 1}: {int(m.sum())} windows, {int(h.sum()):,} counts")
        else:
            ax.plot([], [], " ", label=f"pause {k + 1}: no Felix window inside its {(r['dark1'] - r['dark0']) * 1e3:.1f} ms")
    if irf is not None and h_ref.sum() > 0:
        try:
            t_i, c_i = bd._irf_curve(irf, bins.shape[1] * bin_s)
            shift = t_ns[int(np.argmax(h_ref))] - bd.PEAK_AT_NS           # put the IRF peak on the data peak
            t_i = (t_i + shift) % (bins.shape[1] * bin_s * 1e9)
            o = np.argsort(t_i)
            ax.semilogy(t_i[o], c_i[o] * h_ref.max(), "--", lw=0.8, color="tab:red", label=f"IRF {irf.get('name', '')} (scaled)")
        except Exception:
            pass
    ax.set_xlabel("time in the laser period [ns]")
    ax.set_ylabel("counts per bin (raw)")
    ax.set_ylim(bottom=0.5)
    ax.legend(fontsize=8, loc="upper right")
    for x in axes:
        x.grid(alpha=0.3)
    info = (f"{meta.get('saved', '')}  note: {meta.get('note', '') or '-'}   |   {len(per)} windows of {W * 1e3:.3g} ms, "
            f"{int(per.sum()):,} counts   |   laser-only level {dark_level / W / 1e3:.0f} kcps, max {per.max() / W / 1e6:.1f} Mcps"
            f"   |   before the run: SYNC {float(rates.get('0', 0)) / 1e6:.1f} MHz, SPAD {float(rates.get('1', 0)) / 1e3:.0f} kcps"
            f"   |   windows inside the pauses: {inside}   |   IRF: {irf.get('name') if irf else 'NONE'}")
    fig.text(0.01, 0.004, info, fontsize=7.5, family="monospace")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    if save:
        fig.savefig(Path(run_dir) / "raw.png", dpi=110)
    if show:
        plt.show()
    return fig


def plot_pause_zoom(res: dict, before_ms: float = 40.0, after_ms: float = 40.0, save_to=None, show: bool = True,
                    ymax: float = 1e8):
    """Every pause against the time since LED-off: SPAD total per 5 ms window (black, left axis,
    log, up to ``ymax`` counts/s) and the fluorescence relative yield (AMBIT s/r) of every measuring
    pulse (dots, right axis, coloured by the actinic level of their line: grey in the pause)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    levels = sorted({p["level"] for p in res["pauses"]})
    n_p = max(p["pause"] for p in res["pauses"]) + 1
    fig, axes = plt.subplots(len(levels), n_p, figsize=(4.8 * n_p, 3.6 * len(levels)), squeeze=False)
    ts, W, ws, f = res["ts"], res["W"], res["ws"], res["ref"]
    a = res["ambit"]
    sr_max = 1.15 * float(np.nanmax(a["ratio"]))
    for p in res["pauses"]:
        ax = axes[levels.index(p["level"]), p["pause"]]
        t0, t1 = p["dark0"], p["dark1"]
        lo, hi = t0 - before_ms * 1e-3, t1 + after_ms * 1e-3
        m = (ts + W > lo) & (ts < hi)
        ax.step((ts[m] + W / 2 - t0) * 1e3, ws["total_cps"][m], where="mid", color="k", lw=0.9)
        ax.axvspan(0, (t1 - t0) * 1e3, color="0.85", lw=0)
        ax.set_yscale("log")
        ax.set_ylim(0.1 * f["total_cps"], ymax)
        ax2 = ax.twinx()
        ma = (a["t"] >= lo) & (a["t"] <= hi)
        for lv in np.unique(a["actinic"][ma]):
            mm = ma & (a["actinic"] == lv)
            ax2.plot((a["t"][mm] - t0) * 1e3, a["ratio"][mm], "o", ms=4, color=_colour(int(lv)), mec="k", mew=0.4)
        ax2.set_ylim(0, sr_max)
        ax.set_title(f"actinic {p['level']}, pause {p['pause'] + 1}: {p['dark_ms']:.0f} ms dark"
                     + (" (stays dark)" if p["open"] else ""), fontsize=9)
        ax.set_xlabel("time since LED off [ms] (grey: dark)")
        ax.grid(alpha=0.3)
        if p["pause"] == 0:
            ax.set_ylabel("SPAD total, counts/s per window")
        if p["pause"] == n_p - 1:
            ax2.set_ylabel("fluo. rel. yield (AMBIT s/r)")
        if p["pause"] == 0 and levels.index(p["level"]) == 0:
            ax.legend(handles=[Line2D([], [], color="k", lw=0.9, label="SPAD total"),
                               Line2D([], [], color="0.5", marker="o", ls="", mec="k", label="fluo. rel. yield (colour: actinic)")],
                      fontsize=7, loc="center right")
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=110)
    if show:
        plt.show()
    return fig


def plot_pause_decays(res: dict, save_to=None, show: bool = True):
    """Summed decay of every pause (fully dark windows) against the reference: raw counts per
    second per 47 ps bin, one panel per pause."""
    import matplotlib.pyplot as plt

    n_p = max(p["pause"] for p in res["pauses"]) + 1
    fig, axes = plt.subplots(1, n_p, figsize=(4.8 * n_p, 4.4), squeeze=False)
    f = res["ref"]
    for j in range(n_p):
        ax = axes[0, j]
        curves = [("reference", f, "0.3")] + [(f"actinic {p['level']}", p["summary"], _colour(p["level"]))
                                              for p in res["pauses"] if p["pause"] == j]
        dark_ms = [p["dark_ms"] for p in res["pauses"] if p["pause"] == j]
        for name, s, c in curves:
            if s is None or not s["significant"]:
                continue
            t = (np.arange(s["h"].size) + 0.5) * s["bin_ns"]
            lab = (f"{name}: {s['laser_cps'] / 1e3:.0f} kcps, {s['n_windows']} windows, "
                   f"t_mean {s['t_mean_ns']:.3f} +- {s['t_mean_err_ns']:.3f}")
            ax.semilogy(t, s["h"] / s["seconds"], color=c, lw=0.8, label=lab)
        ax.set_title(f"pause {j + 1} (~{np.mean(dark_ms):.0f} ms)", fontsize=9)
        ax.legend(fontsize=7)
        ax.set_xlim(0, 12.5)
        ax.set_xlabel("micro-time [ns] (laser peak at 2.5 ns)")
        ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel("counts/s per 47 ps bin")
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=110)
    if show:
        plt.show()
    return fig


def fit_groups(res: dict) -> list[tuple[str, list]]:
    """The decays of a run in figure groups: laser only (reference, before and after the run),
    then the pauses of each lit level. Each decay is (name, decay_summary, colour)."""
    laser = [(f"reference: {res['ref']['name']}", res["ref"], "0.3")]
    if res["before"] is not None:
        laser.append(("before the run (laser only)", res["before"], "tab:blue"))
    if res["tail"] is not None:
        laser.append(("0.1-1 s after the run (laser only)", res["tail"], "tab:cyan"))
    groups = [("laser", laser)]
    for lv in sorted({p["level"] for p in res["pauses"]}):
        ps = sorted((p for p in res["pauses"] if p["level"] == lv), key=lambda p: p["pause"])
        groups.append((f"actinic{lv}", [(f"actinic {lv}, pause {p['pause'] + 1}: {p['dark_ms']:.0f} ms dark"
                                          + (" (stays dark)" if p["open"] else ""), p["summary"], _colour(lv))
                                         for p in ps]))
    return groups


def fit_report(res: dict) -> str:
    """Main fit: 2-exponential IRF reconvolution (``reconv.fit``) of every decay, the whole decay
    from 1.5 ns before the peak to 7 ns after it. Per component: the lifetime with its error and
    the share of the photons (``ph``). A lifetime marked ``*`` is below half the IRF FWHM: the IRF
    cannot separate such a component from scatter, so it has no measurable amplitude. ``t_amp``:
    amplitude-weighted lifetime of the components with photons (proportional to the yield; the
    amplitude of an unresolved component is poorly determined). ``t_int``: photon-weighted mean
    lifetime, the model's counterpart of the model-free ``t_mean``. Then the IRF shift and the
    reduced chi2. The amplitudes in counts per bin are in the figures (``plot_fits``)."""
    out = [f"{'decay':<44} {'laser cts':>9} {'t_mean ns':>15} {'tau1 ns':>13} {'ph1':>4} {'tau2 ns':>13} "
           f"{'ph2':>4} {'t_amp ns':>15} {'t_int ns':>15} {'shift ps':>8} {'chi2r':>6}"]
    for _, group in fit_groups(res):
        for name, s, _c in group:
            if s is None or not s["significant"]:
                out.append(f"{name:<44} no significant laser peak")
                continue
            n = s["laser_cps"] * s["seconds"] / s["dead_time_factor"]
            tm = f"{s['t_mean_ns']:.3f} +- {s['t_mean_err_ns']:.3f}"
            r = s.get("rc")
            if r is None or not r.get("ok"):
                out.append(f"{name:<44} {n:>9,.0f} {tm:>15}   no reconvolution fit (no IRF, or too few counts)")
                continue
            comp = [f"{'none':^19} " if r["empty"][i] else
                    f"{r['taus'][i]:>6.3f}+-{r['taus_err'][i]:.3f}{' ' if r['resolved'][i] else '*'} "
                    f"{100 * r['share'][i]:>3.0f}% " for i in range(2)]
            out.append(f"{name:<44} {n:>9,.0f} {tm:>15} {comp[0]}{comp[1]}"
                       f"{r['t_avg_all']:>6.3f} +- {r['t_avg_all_err']:.3f} "
                       f"{r['tau_int']:>6.3f} +- {r['tau_int_err']:.3f} {r['shift_ns'] * 1e3:>+8.1f} "
                       f"{r['chi2r']:>6.2f}" + ("  (parameter at a bound)" if r["at_bound"] else ""))
    r0 = next((s["rc"] for _, g in fit_groups(res) for _, s, _c in g if s is not None and s.get("rc")), None)
    if r0 is not None:
        out.append(f"(IRF {r0['irf_name']}, FWHM {r0['irf_fwhm_ps']:.0f} ps; fit window {r0['window'][0]:.1f}-"
                   f"{r0['window'][1]:.1f} ns, laser peak at 2.5 ns. *: lifetime below {r0['resolved_ns']:.2f} ns = "
                   f"half the IRF FWHM, unresolved: no amplitude. none: a component without photons, the decay "
                   f"is a single exponential. ph: share of the photons)")
    return "\n".join(out)


def levels_report(run_dirs, irf: dict | None = None) -> str:
    """One table over the runs of ``run_levels``: per level, its reference line (t_amp of the
    reconvolution fit, mean arrival, fluorescence relative yield = AMBIT s/r), the lit block
    (SPAD Mcps, yield) and every pause (t_amp +- error, mean arrival, yield). Each run is analysed
    with ``analyse_pauses``; to look at one level in full, analyse and plot its folder."""
    head = (f"{'level':>5} | {'ref t_amp ns':>14} {'t_mean':>6} {'yield':>6} | {'lit Mcps':>8} {'yield':>6} |"
            + "".join(f" {'pause ' + str(j + 1) + ': ms, t_amp ns, t_mean, yield':<40}|" for j in range(3)))
    out = [head]
    for d in run_dirs:
        res = analyse_pauses(d, irf=irf)
        a = res["ambit"]
        lv = max(res["lit"]) if res["lit"] else 0
        ref = res["ref"]
        dark = [r for r in res["lines"] if r["role"] == "dark"]
        y_ref = float(np.median(a["ratio"][a["line"] == dark[0]["line"]])) if dark else float("nan")
        lit_lines = [r["line"] for r in res["lines"] if r["role"] == "lit"]
        y_lit = float(np.median(a["ratio"][np.isin(a["line"], lit_lines)])) if lit_lines else float("nan")
        r0 = ref.get("rc") or {}
        tau = f"{r0['t_avg_all']:.3f}+-{r0['t_avg_all_err']:.3f}" if r0.get("ok") else "-"
        lit = res["lit"].get(lv)
        row = (f"{lv:>5} | {tau:>14} {ref['t_mean_ns']:>6.3f} {y_ref:>6.3f} | "
               f"{(lit['total_cps'] / 1e6 if lit else float('nan')):>8.2f} {y_lit:>6.3f} |")
        for p in sorted(res["pauses"], key=lambda q: q["pause"]):
            sm, rc_ = p["summary"], (p["summary"] or {}).get("rc") or {}
            yp = a["ratio"][(a["t"] > p["dark0"]) & (a["t"] < p["dark1"])]
            cell = (f"{p['dark_ms']:>4.0f} "
                    + (f"{rc_['t_avg_all']:.3f}+-{rc_['t_avg_all_err']:.3f} " if rc_.get("ok") else f"{'-':>12} ")
                    + (f"{sm['t_mean_ns']:.3f} " if sm is not None and sm["significant"] else "  -   ")
                    + (f"{yp.mean():.3f}" if yp.size else "  -  "))
            row += f" {cell:<40}|"
        out.append(row)
        for w in res["warnings"]:
            out.append(f"      {w}")
    out.append("(t_amp: amplitude-weighted lifetime of the 2-exp IRF reconvolution fit; yield: AMBIT s/r, "
               "median for ref and lit, mean of the samples in a pause)")
    return "\n".join(out)


def plot_fits_vs_ref(res: dict, ref: str = "before", save_to=None, show: bool = True, ylim=None,
                     yscale: str = "linear"):
    """Every pause against one laser-only decay, one panel per pause (rows: levels, columns: pauses).
    Each panel: the pause ("Fm n", n = pause number) as points with its reconvolution fit, and behind
    it the laser-only decay ("Fo") as grey points with its fit; no components. Both in counts per
    second per 47 ps bin (raw counts divided by the acquisition time, not normalised), since a pause
    lasts 25-90 ms and Fo about 1 s. The legend gives t_amp (amplitude-weighted lifetime of the
    reconvolution fit) and t_int (photon-weighted) of both, and phi2_lifetime = (t_amp,Fm -
    t_amp,Fo) / t_amp,Fm with its error. ``ref`` (the Fo decay): "before"
    (before the run, laser only), "reference" (the actinic-0 line) or "after" (0.1-1 s after the run)."""
    import matplotlib.pyplot as plt

    names = {"before": ("before the run", res.get("before")), "reference": ("reference line", res.get("ref")),
             "after": ("after the run", res.get("tail"))}
    ref_name, rs = names[ref]
    if rs is None or not (rs.get("rc") or {}).get("ok"):
        ref_name, rs = "reference line", res["ref"]
    rr = rs["rc"]
    levels = sorted({p["level"] for p in res["pauses"]})
    n_p = max(p["pause"] for p in res["pauses"]) + 1
    fig, axes = plt.subplots(len(levels), n_p, figsize=(4.9 * n_p, 4.2 * len(levels)), squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    t = rr["t"]
    win = rr["fit_mask"]
    for p in res["pauses"]:
        ax = axes[levels.index(p["level"]), p["pause"]]
        ax.set_visible(True)
        s = p["summary"]
        r = (s or {}).get("rc") or {}
        to, eo = rr["t_avg_all"], rr["t_avg_all_err"]
        ax.plot(t, rr["h"] / rs["seconds"], "o", ms=2.2, color="0.72", mec="none", label="Fo")
        ax.plot(t[win], rr["model"][win] / rs["seconds"], color="0.5", lw=1.4,
                label=r"Fo: fit, $t_\mathrm{amp}$" + f" = {to:.3f} \u00b1 {eo:.3f} ns, "
                      + r"$t_\mathrm{int}$" + f" = {rr['tau_int']:.3f} ns")
        top = (rr["h"][win] / rs["seconds"]).max()
        title = f"actinic {p['level']}, pause {p['pause'] + 1}: {p['dark_ms']:.0f} ms dark"
        k = p["pause"] + 1
        if r.get("ok"):
            tm, em = r["t_avg_all"], r["t_avg_all_err"]
            ax.plot(t, r["h"] / s["seconds"], "o", ms=2.2, color=_colour(p["level"]), mec="none", label=f"Fm {k}")
            ax.plot(t[win], r["model"][win] / s["seconds"], color="k", lw=1.3,
                    label=f"Fm {k}: fit, " + r"$t_\mathrm{amp}$" + f" = {tm:.3f} \u00b1 {em:.3f} ns, "
                          + r"$t_\mathrm{int}$" + f" = {r['tau_int']:.3f} ns")
            top = max(top, (r["h"][win] / s["seconds"]).max())
            phi, phi_err = _phi2(to, eo, tm, em)
            ax.legend(fontsize=7, loc="upper right" if yscale != "log" else "lower left", title_fontsize=7.5, markerscale=2.5,
                      title=r"$\phi_{2,\mathrm{lifetime}} = (t_{\mathrm{amp,Fm}} - t_{\mathrm{amp,Fo}})\,/\,"
                            r"t_{\mathrm{amp,Fm}}$" + f" = {phi:.3f} \u00b1 {phi_err:.3f}")
        else:
            ax.text(0.5, 0.5, "no fit (too few counts)", transform=ax.transAxes, ha="center", fontsize=9)
            ax.legend(fontsize=7, loc="upper right")
        if yscale == "log":
            ax.set_yscale("log")
            ax.set_ylim(*(ylim or (0.5 * min(rs["floor_per_bin"] / rs["seconds"], top / 1e3), 2 * top)))
        else:
            ax.set_ylim(0, ylim[1] if ylim else 1.12 * top)
        ax.set_xlim(0, t[-1] + t[0])
        ax.set_title(title + (" (stays dark)" if p["open"] else ""), fontsize=9)
        ax.set_xlabel("micro-time [ns] (laser peak at 2.5 ns)")
        if p["pause"] == 0:
            ax.set_ylabel("counts/s per 47 ps bin")
        ax.grid(alpha=0.3)
    fig.suptitle(f"{Path(res['run_dir']).name}: every pause (Fm n) against Fo = {ref_name} (grey), raw rates",
                 fontsize=10)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=110)
    if show:
        plt.show()
    return fig


def plot_fits(res: dict, irf: dict | None = None, save: bool = True, show: bool = True, ylim=None,
              yscale: str = "linear"):
    """Main fit figures: raw decay and the 2-exponential IRF reconvolution fit (``reconv.plot``) of
    every decay, each component convolved with the IRF under the data, with its height at t0 (the
    laser peak) and photon share, and t_int; one figure per group of ``fit_groups``, saved as
    pauses_fits_<group>.png. (Every pause against Fo, with phi2_lifetime, is the third figure of
    ``plot_pauses``: ``plot_fits_vs_ref``.) ``yscale``: "linear" (component heights at t0 add up by eye) or "log" (the tails). ``ylim``: the
    same count range in every panel; None: each panel its own. ``irf`` is not needed: the fits use the
    IRF given to ``analyse_pauses``."""
    d = Path(res["run_dir"])
    figs = []
    for key, group in fit_groups(res):
        decays = [(name, s.get("rc") if s is not None else None, c) for name, s, c in group]
        figs.append(rc.plot(decays, save_to=d / f"pauses_fits_{key}.png" if save else None, show=show, ylim=ylim,
                            yscale=yscale))
    return figs


def plot_tail_fits(res: dict, irf: dict | None = None, save: bool = True, show: bool = True, ylim=(1, 1e5)):
    """The old 2-exponential tail fits (``blank_dose.plot_fits``, from 0.25 ns after the peak), for
    comparison: saved as pauses_tailfits_<group>.png."""
    d = Path(res["run_dir"])
    return [bd.plot_fits(group, irf, save_to=d / f"pauses_tailfits_{key}.png" if save else None, show=show,
                         models=(2,), ylim=ylim)
            for key, group in fit_groups(res)]


def plot_pauses(res: dict, save: bool = True, show: bool = True, zoom_ymax: float = 1e8, ref: str = "before",
                yscale: str = "linear"):
    """The three figures: overview, pause zoom, and every pause (Fm n) against Fo with its fit and
    phi2_lifetime (``plot_fits_vs_ref``); saved in the run folder as pauses_overview.png,
    pauses_zoom.png, pauses_fits_vs_<ref>.png. ``zoom_ymax``: top of the SPAD axis in the pause zoom,
    counts/s. ``ref``: the Fo decay ("before", "reference" or "after"); ``yscale``: of the third
    figure. The raw pause decays without fits: ``plot_pause_decays``."""
    d = Path(res["run_dir"])
    sv = (lambda n: d / f"pauses_{n}.png") if save else (lambda n: None)
    figs = [plot_overview(res, sv("overview"), show), plot_pause_zoom(res, save_to=sv("zoom"), show=show, ymax=zoom_ymax)]
    if res["pauses"]:
        figs.append(plot_fits_vs_ref(res, ref, save_to=sv(f"fits_vs_{ref}"), show=show, yscale=yscale))
    return figs


# ---------------------------------------------------------------------------
# Light series: pre-illumination at several levels, one pauses run per level
# ---------------------------------------------------------------------------
PRE_LEVELS = (0, 10, 20, 40, 80, 160, 255)   # pre-illumination actinic bytes, in run order
T_WAIT_S = 120.0                             # LED on, AMBIT idle, before each measurement
FO_POINTS = 50                               # AMBIT Fo: mean of the last n dark points before the train
FM_POINTS = 50                               # AMBIT Fm: mean of the last n lit points of the train
SERIES_MARGIN_S = 0.2                        # Felix margin before and after the arrun (not analysed)


def light_on(ambit: ag.AmbitLink, actinic: int) -> None:
    """A one-sample arrun at ``actinic`` with the persist flag: one measuring flash, then the LED
    stays on at that level while the AMBIT is idle (no more flashes until the next command). An
    ``actinic`` <= 3 leaves the LED off (the firmware switches it off at the start of such a line)."""
    prot = [ac.line(0.1, 10, int(actinic))]
    ag.run_arrun(ambit, ag.build_arrun(prot, persist=1), prot, log=lambda *a: None)


def light_off(ambit: ag.AmbitLink) -> None:
    """Switch a persistent LED off: a one-sample dark arrun without persist."""
    ac.led_off(ambit)


def _write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float), encoding="utf-8")


def run_light_series(ambit: ag.AmbitLink, felix: ag.FelixLink, session, label: str, levels=PRE_LEVELS,
                     t_wait_s: float = T_WAIT_S, protocol: list[dict] | None = None,
                     pre_s: float = SERIES_MARGIN_S, tail_s: float = SERIES_MARGIN_S, **kw) -> list[Path]:
    """Per level X of ``levels``, in order: ``light_on(X)``, ``t_wait_s`` of waiting with the AMBIT
    idle, ``light_off``, then one ``run_pauses`` of ``protocol`` (the train of the measurement cell:
    dark line, lit blocks with their pauses, dark line) with only ``pre_s``/``tail_s`` of Felix
    margin, so the leaf is in the dark for the 1 s dark line (plus ~0.5 s of set-up) before the
    first 255 block and for the dark line after the last one. Then straight on to the next level.
    The SPAD rate is logged at the start and the end of every wait: with the LED on it must not
    drop. Each run goes to ``session``/pauses_<label>_pre<X>_<time> with its pre-illumination in
    run.json; ``session``/series.json lists the runs (updated after each one, so an aborted series
    can still be analysed). ``kw`` goes to ``run_pauses`` (window_s, note, threshold_v, laser_hz,
    log). Whatever happens, the LED is switched off at the end."""
    log = kw.get("log", print)
    if not protocol:
        raise ValueError("protocol: the arrun lines of one measurement (train_protocol)")
    ag.build_arrun(protocol)                                   # validated before the first wait
    levels = [int(x) for x in levels]
    session = Path(session)
    session.mkdir(parents=True, exist_ok=True)
    meta = dict(kind="light_series", label=label, levels=levels, t_wait_s=float(t_wait_s), protocol=protocol,
                pre_s=pre_s, tail_s=tail_s, runs=[], started=time.strftime("%Y-%m-%d %H:%M:%S"))
    dirs: list[Path] = []
    try:
        for i, x in enumerate(levels):
            log(f"--- pre-illumination {x} ({i + 1} of {len(levels)}): LED {'on' if x > ag.ACTINIC_ON else 'off'}, "
                f"waiting {t_wait_s:.0f} s")
            light_on(ambit, x)
            r0 = felix.count_rates(0.5).get(ag.SPAD_CH, float("nan"))
            t_end = time.monotonic() + t_wait_s - 1.0              # the two rate readings take 0.5 s each
            while time.monotonic() < t_end:
                time.sleep(min(1.0, t_end - time.monotonic()))
            r1 = felix.count_rates(0.5).get(ag.SPAD_CH, float("nan"))
            log(f"    SPAD {r0 / 1e3:,.0f} kcps at the start of the wait, {r1 / 1e3:,.0f} kcps at the end")
            if x > ag.ACTINIC_ON and r1 < 0.5 * r0:
                log("    WARNING: the SPAD rate fell during the wait: did the LED stay on?")
            light_off(ambit)
            d = run_pauses(ambit, felix, session / f"pauses_{label}_pre{x:03d}_{time.strftime('%H%M%S')}", protocol,
                           pre_s=pre_s, tail_s=tail_s,
                           extra=dict(pre_illumination=dict(actinic=x, wait_s=float(t_wait_s),
                                                            spad_kcps_start=r0 / 1e3, spad_kcps_end=r1 / 1e3)), **kw)
            dirs.append(d)
            meta["runs"].append(d.name)
            _write_json(session / "series.json", meta)
    finally:
        try:
            light_off(ambit)
        except Exception as exc:                                   # the LED must not stay on
            log(f"could not switch the LED off: {exc}")
    meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    _write_json(session / "series.json", meta)
    return dirs


def _mean_se(v) -> tuple[float, float]:
    """Mean and standard error of the finite values (se nan with fewer than 2)."""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan")
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else float("nan")


def _phi2(tau_o: float, err_o: float, tau_m: float, err_m: float) -> tuple[float, float]:
    """phi2 = 1 - tau_o / tau_m and its error from the two lifetime errors."""
    if not (np.isfinite(tau_o) and np.isfinite(tau_m)) or tau_o <= 0 or tau_m <= 0:
        return float("nan"), float("nan")
    return 1.0 - tau_o / tau_m, (tau_o / tau_m) * float(np.hypot(err_o / tau_o, err_m / tau_m))


AMBIT_FLOOR_COUNTS = 5.0     # s_630 median below this: the AMBIT fluorescence signal sits at its floor


def _fit_numbers(s: dict | None) -> dict | None:
    """The numbers of one decay's reconvolution fit for a table: laser-correlated counts, both
    lifetimes with errors, resolved and empty flags and photon shares, t_amp, t_int, IRF shift,
    chi2r. None without a fit."""
    r = (s or {}).get("rc") or {}
    if not r.get("ok"):
        return None
    return dict(counts=float(s["laser_cps"] * s["seconds"] / s["dead_time_factor"]), seconds=float(s["seconds"]),
                tau1=float(r["taus"][0]), tau1_err=float(r["taus_err"][0]), resolved1=bool(r["resolved"][0]),
                empty1=bool(r["empty"][0]), share1=float(r["share"][0]),
                tau2=float(r["taus"][1]), tau2_err=float(r["taus_err"][1]), resolved2=bool(r["resolved"][1]),
                empty2=bool(r["empty"][1]), share2=float(r["share"][1]),
                t_amp=float(r["t_avg_all"]), t_amp_err=float(r["t_avg_all_err"]), t_int=float(r["tau_int"]),
                t_int_err=float(r["tau_int_err"]), shift_ps=float(r["shift_ns"] * 1e3), chi2r=float(r["chi2r"]),
                at_bound=bool(r["at_bound"]), t_mean=float(s["t_mean_ns"]), t_mean_err=float(s["t_mean_err_ns"]))


def run_summary(res: dict, fo_points: int = FO_POINTS, fm_points: int = FM_POINTS) -> dict:
    """phi2 of one analysed run (``analyse_pauses``), from the AMBIT and from the lifetimes, with
    every division guarded and a list of ``warnings`` in plain words instead of inf, nan or a
    silently wrong number.

    AMBIT: Fo = mean s/r of the last ``fo_points`` points of the first dark line, Fm = mean of the
    last ``fm_points`` lit points (actinic > 3) of the train, phi2 = (Fm - Fo) / Fm, error from the
    scatter of the points. Warned: no points, Fm = 0 (phi2 not computed), Fo = 0 (phi2 = 1, no
    error), Fm not above Fo (phi2 <= 0), the raw s_630 signal at its floor (phi2 meaningless).
    Lifetimes: t_amp (the amplitude-weighted lifetime of the reconvolution fit, proportional to
    the yield) of the reference line (Fo) and of every pause (Fm), phi2_tau = (t_amp,Fm -
    t_amp,Fo) / t_amp,Fm per pause, mean and SD; t_int (photon-weighted) and tau2 (the slow
    component) are kept for the table but not used. Warned: no fit of the reference or of a
    pause, a pause with t_amp not above the reference. ``fits`` holds the fit numbers of the
    reference and the pauses for a table (``_fit_numbers``)."""
    a = res["ambit"]
    warnings: list[str] = []
    dark0 = next((r for r in res["lines"] if r["role"] == "dark"), None)
    fo_idx = np.flatnonzero(a["line"] == dark0["line"])[-int(fo_points):] if dark0 else np.zeros(0, int)
    fm_idx = np.flatnonzero(a["actinic"] > ag.ACTINIC_ON)[-int(fm_points):]
    if dark0 is None:
        warnings.append("no dark line before the train: AMBIT Fo not defined")
    elif fo_idx.size == 0:
        warnings.append("AMBIT Fo: no points in the dark line")
    if fm_idx.size == 0:
        warnings.append("AMBIT Fm: no lit points in the run")
    fo, se_o = _mean_se(a["ratio"][fo_idx])
    fm, se_m = _mean_se(a["ratio"][fm_idx])
    phi_a, phi_a_err = float("nan"), float("nan")
    if not np.isfinite(fm) or fm <= 0:
        if fm_idx.size:
            warnings.append(f"AMBIT Fm is {fm:.4g}: phi2 = (Fm - Fo) / Fm cannot be computed (division by zero)")
    elif not np.isfinite(fo) or fo <= 0:
        phi_a = 1.0
        warnings.append(f"AMBIT Fo is {fo:.4g}: phi2 set to 1, its error is not defined")
    else:
        phi_a = 1.0 - fo / fm
        phi_a_err = (fo / fm) * float(np.hypot(se_o / fo, se_m / fm))
        if fm <= fo:
            warnings.append(f"AMBIT Fm ({fm:.4f}) is not above Fo ({fo:.4f}): phi2 = {phi_a:.2f} is not meaningful")
    s_med = {}
    if "s_630" in a:
        for name, idx in (("Fo", fo_idx), ("Fm", fm_idx)):
            if idx.size:
                s_med[name] = float(np.median(a["s_630"][idx]))
                if s_med[name] < AMBIT_FLOOR_COUNTS:
                    warnings.append(f"AMBIT fluorescence signal at its floor in {name} (s_630 median "
                                    f"{s_med[name]:.0f} counts): the AMBIT phi2 is meaningless")
    ref_fit = _fit_numbers(res["ref"])
    if ref_fit is None:
        warnings.append("no reconvolution fit of the reference decay (no IRF, or too few counts): "
                        "phi2 from the lifetimes not defined")
    nan = float("nan")
    to, eo = (ref_fit["t_amp"], ref_fit["t_amp_err"]) if ref_fit else (nan, nan)       # phi2 uses t_amp
    ti_o, ei_o = (ref_fit["t_int"], ref_fit["t_int_err"]) if ref_fit else (nan, nan)
    t2_o, e2_o = (ref_fit["tau2"], ref_fit["tau2_err"]) if ref_fit else (nan, nan)
    row = dict(run=Path(res["run_dir"]).name, note=res.get("note", ""), ref_name=res["ref"]["name"],
               fo_ambit=fo, fo_ambit_se=se_o, n_fo=int(fo_idx.size), fm_ambit=fm, fm_ambit_se=se_m,
               n_fm=int(fm_idx.size), s_630_fo_median=s_med.get("Fo", float("nan")),
               s_630_fm_median=s_med.get("Fm", float("nan")), phi2_ambit=phi_a, phi2_ambit_err=phi_a_err,
               tau2_fo=t2_o, tau2_fo_err=e2_o, t_amp_fo=to, t_amp_fo_err=eo, t_int_fo=ti_o, t_int_fo_err=ei_o,
               laser_fo_kcps=res["ref"]["laser_cps"] / 1e3,
               lit_mcps=max((s["total_cps"] for s in res["lit"].values()), default=float("nan")) / 1e6,
               tau2_fm=[], tau2_fm_err=[], t_amp_fm=[], t_amp_fm_err=[], t_int_fm=[], t_int_fm_err=[],
               phi2_tau=[], phi2_tau_err=[],
               dark_ms=[], fits=dict(reference=ref_fit, pauses=[]))
    for p in sorted(res["pauses"], key=lambda q: q["pause"]):
        k = p["pause"] + 1
        f = _fit_numbers(p["summary"])
        row["fits"]["pauses"].append(f)
        tm, em = (f["t_amp"], f["t_amp_err"]) if f else (nan, nan)
        if f is None:
            warnings.append(f"pause {k}: no fit ({p['n_windows']} Felix windows in {p['dark_ms']:.0f} ms dark)")
        phi, perr = _phi2(to, eo, tm, em)
        if f is not None and ref_fit is not None and np.isfinite(phi) and phi <= 0:
            warnings.append(f"pause {k}: t_amp Fm ({tm:.3f} ns) is not above t_amp Fo ({to:.3f} ns): phi2 = {phi:.2f}")
        row["tau2_fm"].append(f["tau2"] if f else nan)
        row["tau2_fm_err"].append(f["tau2_err"] if f else nan)
        row["t_amp_fm"].append(tm)
        row["t_amp_fm_err"].append(em)
        row["t_int_fm"].append(f["t_int"] if f else nan)
        row["t_int_fm_err"].append(f["t_int_err"] if f else nan)
        row["phi2_tau"].append(phi)
        row["phi2_tau_err"].append(perr)
        row["dark_ms"].append(p["dark_ms"])
    if not res["pauses"]:
        warnings.append("no pause in the run: no Fm lifetime")
    v = np.array(row["phi2_tau"], float)
    v = v[np.isfinite(v)]
    row["phi2_tau_mean"] = float(v.mean()) if v.size else float("nan")
    row["phi2_tau_sd"] = float(v.std(ddof=1)) if v.size > 1 else float("nan")
    t2 = np.array(row["tau2_fm"], float)
    t2 = t2[np.isfinite(t2)]
    row["tau2_fm_mean"] = float(t2.mean()) if t2.size else float("nan")
    row["tau2_fm_sd"] = float(t2.std(ddof=1)) if t2.size > 1 else float("nan")
    for key in ("t_amp_fm", "t_int_fm"):
        ti = np.array(row[key], float)
        ti = ti[np.isfinite(ti)]
        row[key + "_mean"] = float(ti.mean()) if ti.size else float("nan")
        row[key + "_sd"] = float(ti.std(ddof=1)) if ti.size > 1 else float("nan")
    row["warnings"] = warnings + list(res["warnings"])
    return row


def analyse_light_series(src, irf: dict | None = None, fo_points: int = FO_POINTS, fm_points: int = FM_POINTS,
                         log=print) -> dict:
    """``src``: a session folder with series.json, or a list of run folders. Every run is analysed
    with ``analyse_pauses``; per run: AMBIT Fo = mean s/r of the last ``fo_points`` points of the
    first dark line (the 1 s before the train), Fm = mean of the last ``fm_points`` lit points
    (actinic > 3) of the train, phi2 = (Fm - Fo) / Fm with the error from the scatter of the points;
    t_amp, t_int and tau2 of the reference line (Fo) and of every pause (Fm); phi2_tau = (t_amp,Fm -
    t_amp,Fo) / t_amp,Fm per pause, with their mean and SD. Rows are in run order with the
    pre-illumination from run.json (nan if the run has none). Saves series.csv and
    series_results.json in the session folder. Returns dict(session, rows, results)."""
    if isinstance(src, (str, Path)) and (Path(src) / "series.json").exists():
        session = Path(src)
        dirs = [session / n for n in json.loads((session / "series.json").read_text(encoding="utf-8"))["runs"]]
    else:
        dirs = [Path(d) for d in src]
        session = dirs[0].parent if dirs else None
    dirs = [d for d in dirs if (d / "felix.npz").exists()]
    if not dirs:
        raise ValueError("no run with data")
    irf = irf if irf is not None else bd.load_irf()
    rows, results = [], []
    for d in dirs:
        meta = json.loads((d / "run.json").read_text(encoding="utf-8"))
        pre = meta.get("pre_illumination") or {}
        res = analyse_pauses(d, irf=irf)
        results.append(res)
        row = dict(actinic=int(pre["actinic"]) if "actinic" in pre else float("nan"),
                   wait_s=float(pre.get("wait_s", float("nan"))))
        row.update(run_summary(res, fo_points, fm_points))
        rows.append(row)
        log(f"{d.name}: pre {row['actinic']}, Fo {row['fo_ambit']:.4f} Fm {row['fm_ambit']:.4f} phi2 "
            f"{row['phi2_ambit']:.3f}; t_amp Fo {row['t_amp_fo']:.3f} ns, Fm {np.round(row['t_amp_fm'], 3).tolist()}, "
            f"phi2_tau {row['phi2_tau_mean']:.3f} +- {row['phi2_tau_sd']:.3f}"
            + "".join(f"\n    WARNING: {w}" for w in row["warnings"]))
    ser = dict(session=str(session) if session else "", rows=rows, results=results,
               fo_points=int(fo_points), fm_points=int(fm_points))
    if session is not None:
        _write_json(session / "series_results.json",
                    dict(session=ser["session"], fo_points=fo_points, fm_points=fm_points, rows=rows))
        n_p = max((len(r["phi2_tau"]) for r in rows), default=0)
        scalar = [k for k in rows[0] if not isinstance(rows[0][k], (list, dict))]
        lists = ("tau2_fm", "tau2_fm_err", "t_amp_fm", "t_amp_fm_err", "t_int_fm", "t_int_fm_err", "phi2_tau",
                 "phi2_tau_err", "dark_ms")
        with open(session / "series.csv", "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(scalar + [f"{k}_{j + 1}" for k in lists for j in range(n_p)])
            for r in rows:
                wr.writerow([r[k] for k in scalar]
                            + [r[k][j] if j < len(r[k]) else "" for k in lists for j in range(n_p)])
    return ser


def series_report(ser: dict) -> str:
    """One line per level: AMBIT Fo, Fm, phi2; t_amp and t_int of the reference line; t_amp and t_int
    of every pause; phi2 from the lifetimes per pause and their mean +- SD."""
    out = [f"{'pre':>4} {'wait s':>6} {'Fo_A':>7} {'Fm_A':>7} {'phi2_A':>13} | {'t_amp Fo ns':>13} {'t_int Fo':>8} | "
           f"{'t_amp Fm ns, pauses':<30} {'t_int Fm':<24} | {'phi2_tau, pauses':<28} {'mean +- sd':>14}"]
    for r in ser["rows"]:
        pre = f"{r['actinic']:>4}" if np.isfinite(r["actinic"]) else "   -"
        out.append(f"{pre} {r['wait_s']:>6.0f} {r['fo_ambit']:7.4f} {r['fm_ambit']:7.4f} "
                   f"{r['phi2_ambit']:6.3f}+-{r['phi2_ambit_err']:.3f} | {r['t_amp_fo']:6.3f}+-{r['t_amp_fo_err']:.3f} "
                   f"{r['t_int_fo']:8.3f} | {' '.join(f'{t:.3f}' for t in r['t_amp_fm']):<30} "
                   f"{' '.join(f'{t:.3f}' for t in r['t_int_fm']):<24} | "
                   f"{' '.join(f'{p:.3f}' for p in r['phi2_tau']):<28} {r['phi2_tau_mean']:6.3f} +- {r['phi2_tau_sd']:.3f}"
                   + (f"   {len(r['warnings'])} warning(s)" if r["warnings"] else ""))
    out.append(f"(Fo_A: mean s/r of the last {ser['fo_points']} dark points before the train; Fm_A: mean of the last "
               f"{ser['fm_points']} lit points; phi2_A = (Fm - Fo) / Fm, error from the scatter of the points. "
               "t_amp: amplitude-weighted lifetime of the reconvolution fit, Fo from the reference line, Fm per "
               "pause; phi2_tau = (t_amp,Fm - t_amp,Fo) / t_amp,Fm. t_int: photon-weighted mean lifetime, for "
               "comparison)")
    return "\n".join(out)


def plot_light_series(ser: dict, save_to=None, show: bool = True):
    """Three panels over the pre-illumination level: phi2 (AMBIT, and the mean of the 4 lifetime
    values); t_amp of Fo (left axis) and of every Fm pause (right axis, the mean joined); and
    phi2 from the lifetimes against phi2 from the AMBIT, the mean of the 4 pauses as the marker
    (bar: their SD) with the 4 values as faint points. Saved as series_phi2.png in the session
    folder unless ``save_to`` says otherwise (False: not saved)."""
    import matplotlib.pyplot as plt

    rows = ser["rows"]
    x = np.array([r["actinic"] for r in rows], float)
    xlabel = "pre-illumination (actinic byte)"
    if not np.all(np.isfinite(x)):
        x = np.arange(len(rows), dtype=float)
        xlabel = "run"
    n_p = max((len(r["phi2_tau"]) for r in rows), default=0)
    pad = lambda key: np.array([r[key] + [np.nan] * (n_p - len(r[key])) for r in rows], float)
    phi_a = np.array([r["phi2_ambit"] for r in rows])
    phi_a_err = np.array([r["phi2_ambit_err"] for r in rows])
    phi_t, phi_t_err = pad("phi2_tau"), pad("phi2_tau_err")
    phi_t_mean = np.array([r["phi2_tau_mean"] for r in rows])
    phi_t_sd = np.array([r["phi2_tau_sd"] for r in rows])
    t2o, t2o_err = np.array([r["t_amp_fo"] for r in rows]), np.array([r["t_amp_fo_err"] for r in rows])
    t2m, t2m_err = pad("t_amp_fm"), pad("t_amp_fm_err")
    t2m_mean = np.array([r["t_amp_fm_mean"] for r in rows])
    t2m_sd = np.array([r["t_amp_fm_sd"] for r in rows])
    dx = 0.012 * (np.ptp(x) if np.ptp(x) > 0 else 1.0)

    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    ax = axs[0]
    ax.errorbar(x, phi_a, yerr=phi_a_err, fmt="o-", color="tab:blue", capsize=2, label="AMBIT: (Fm - Fo) / Fm")
    ax.errorbar(x, phi_t_mean, yerr=phi_t_sd, fmt="s--", color="tab:red", mfc="none", capsize=2,
                label="lifetime: mean of the 4 pauses (bar: SD)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("phi2")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("phi2 against the pre-illumination", fontsize=10)

    ax = axs[1]
    ax2 = ax.twinx()
    h1 = ax.errorbar(x, t2o, yerr=t2o_err, fmt="o-", color="tab:blue", capsize=2, label="t_amp Fo (reference line)")
    for j in range(n_p):
        ax2.errorbar(x + (j - (n_p - 1) / 2) * dx, t2m[:, j], yerr=t2m_err[:, j], fmt=".", color="tab:red",
                     alpha=0.45, ms=7, capsize=0)
    h2 = ax2.errorbar(x, t2m_mean, yerr=t2m_sd, fmt="s-", color="tab:red", mfc="none", capsize=2,
                      label="t_amp Fm (mean of the pauses; dots: each pause)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("t_amp Fo [ns]", color="tab:blue")
    ax2.set_ylabel("t_amp Fm [ns]", color="tab:red")
    ax.tick_params(axis="y", colors="tab:blue")
    ax2.tick_params(axis="y", colors="tab:red")
    ax.grid(alpha=0.3)
    ax.legend(handles=[h1, h2], fontsize=8, loc="upper left")
    ax.set_title("lifetimes against the pre-illumination", fontsize=10)

    ax = axs[2]
    for i in range(len(rows)):
        ax.plot(np.full(n_p, phi_a[i]), phi_t[i], ".", color="tab:red", alpha=0.35, ms=8)
    ax.errorbar(phi_a, phi_t_mean, xerr=phi_a_err, yerr=phi_t_sd, fmt="s", color="tab:red", ms=7, capsize=2,
                label="mean of the 4 pauses (bar: SD); dots: each pause")
    for i, r in enumerate(rows):
        if np.isfinite(r["actinic"]) and np.isfinite(phi_a[i]) and np.isfinite(phi_t_mean[i]):
            ax.annotate(str(int(r["actinic"])), (phi_a[i], phi_t_mean[i]), fontsize=7, xytext=(4, 4),
                        textcoords="offset points")
    ax.plot([0, 1], [0, 1], "k--", lw=0.7, label="1:1")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.5)
    ax.set_xlabel("phi2 AMBIT")
    ax.set_ylabel("phi2 from the lifetimes, (t_amp,Fm - t_amp,Fo) / t_amp,Fm")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("phi2: lifetime against AMBIT (labels: pre-illumination)", fontsize=10)
    fig.suptitle(Path(ser["session"]).name if ser["session"] else "light series", fontsize=10)
    fig.tight_layout()
    if save_to is None and ser["session"]:
        save_to = Path(ser["session"]) / "series_phi2.png"
    if save_to:
        fig.savefig(save_to, dpi=110)
    if show:
        plt.show()
    return fig
