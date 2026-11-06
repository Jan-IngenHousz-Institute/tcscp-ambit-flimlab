"""AMBIT chopped actinic + Felix time tags: TCSPC lifetimes of Fo, Fs and Fm' in dark gaps.

Why
---
While the AMBIT's actinic light or saturating pulse is on, the leaf fluorescence that the
LED excites reaches the SPAD at about 40 times the laser-excited fluorescence (runs of
3 Nov 2026). An ND filter attenuates both by the same factor, so it only trades SPAD
saturation for photon starvation. Instead, the AMBIT firmware (``chop`` verb, branch
``feat/chop-gaps`` of the ambit repo) cuts short dark gaps into every lit line, and only
the SPAD photons inside those gaps are kept. In a gap the floor is the stray-light leak
plus a little delayed fluorescence, and the laser can run near 1 Mcps.

How
---
1. ``arm_chop`` arms the gaps for the next ``arrun`` (one-shot, on the AMBIT).
2. The Felix time-tags every SPAD event (TTM module, binary file), plus the laser SYNC
   when the TTM can deliver it (``selftest`` finds out which mode works).
3. ``analyse`` reads the file in chunks, places the protocol lines on the Felix clock,
   finds the gaps in the SPAD rate (one matched-filter dip per gap), rebuilds every
   photon's micro-time (from SYNC tags, or from the photons' own laser periodicity when
   SYNC is absent), and builds decays of Fo (dark), Fs (actinic line) and Fm' (pulse).
   Fits reuse ``fvfm_lifetime.fit_window`` and ``phi2_lifetime.global_fit``.

A lit line whose gaps cannot be found reliably (low actinic: the LED-excited floor is then
comparable to the laser signal anyway) is analysed with all its photons and the flat floor
subtracted; the report says which lines were gated. No cable joins the AMBIT and the Felix.

Command line::

    python ambit_chop.py --check                 # AMBIT chop verb, Felix, TTM module
    python ambit_chop.py --analyse runs/<dir>    # re-analyse a saved run
    python ambit_chop.py --simulate runs/sim     # synthetic run through the same analysis
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import ambit_gate as ag  # noqa: E402  AMBIT link, arrun builder, FelixLink
import fvfm_lifetime as fl  # noqa: E402  single-histogram tail fits
import phi2_lifetime as p2  # noqa: E402  shared-lifetime fit of two histograms

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
GAP_US = 150  # requested dark time per gap: firmware chop,<gap_us>,<period_us>
CHOP_PERIOD_US = 1000  # gap spacing inside one ADPD period
WRITE_US = 118.0  # one AS7341 LED-register write (AmbitV003, 3 Nov 2026: dark 265-272 us at gap 150)
GUARD_START_US = 10.0  # skipped after the LED goes dark: LED tail, SPAD recovery, afterpulses
GUARD_END_US = 8.0  # skipped before the on-write lands
DEAD_TIME_S = 200e-9  # SPAD non-paralyzable dead time fitted on 3 Nov 2026
N_MICRO = 5333  # micro-time bins per laser period: 2.34 ps at 80 MHz, as the hardware histogram
PEAK_AT_NS = 2.5  # every decay is rotated so the laser peak sits here (clean floor before it)
LINE_EDGE_S = 0.010  # ignored at both ends of every line
DARK_AFTER_LIT_S = 0.005  # ignored in addition at the start of a dark line that follows light
BIN_US = 1.0  # rate trace resolution
CHUNK = 4_000_000  # TTM records per chunk
REBIN = fl.REBIN  # 20 x 2.34 ps = 47 ps for kinetics and fits
TTM_MAGIC = 0x69B58C9FF09A8CE2
REC = np.dtype([("ce", "u1"), ("t", "<u8")])  # 9-byte TTM record: edge<<7 | channel, timestamp

# Firmware timing (ambit src/PAM.cpp run_arr_trigger) used to predict line positions and
# gap counts. Only predictions: the analysis measures the real ones.
LINE_SETUP_S = 0.0062  # line start -> first ADPD edge: SPI/I2C set-up + 5 ms arm settle
FIRST_LINE_WARMUP_S = 0.012  # three dark warm-up sequences before the first line of a run
LINE_END_S = 0.0006  # last read-out -> LED off
LINE_TRANSITION_S = 0.0013  # LED dark between two lit lines (AS_LED_OFF ... AS_LED_ON)
GAP_SETTLE_US = 50.0  # kChopSettleUs
GAP_EDGE_GUARD_US = 100.0  # kChopEdgeGuardUs
GAP_WRITE_PLAN_US = 120.0  # kChopWritePlanUs
MIN_REGULAR_FRACTION = 0.6  # a lit line is gated only if this fraction of its planned gaps is found


# ---------------------------------------------------------------------------
# Protocols
# ---------------------------------------------------------------------------
def line(seconds: float, freq: int, actinic: int) -> dict:
    """One arrun line: type 2 (no 730 nm reference pulses), ``seconds`` at ``freq`` Hz."""
    return dict(
        type=2,
        far_red=0,
        pulses=max(1, int(round(seconds * freq))),
        freq=int(freq),
        actinic=int(actinic),
        subsampling=0,
    )


def fvfm_protocol(dark_s=3.0, pulse_s=0.8, post_s=1.0, pulse=255, f_dark=10, f_pulse=100) -> list[dict]:
    """Dark-adapted leaf: Fo in the dark, Fm in the chopped saturating pulse, relaxation."""
    return [line(dark_s, f_dark, 0), line(pulse_s, f_pulse, pulse), line(post_s, f_dark, 0)]


def phi2_protocol(
    actinic: int,
    fs_s=3.0,
    pulse_s=0.8,
    after_s=1.0,
    dark_s=0.0,
    post_dark_s=0.5,
    pulse=255,
    f_lit=10,
    f_pulse=100,
    f_dark=10,
) -> list[dict]:
    """Light-adapted leaf: Fs in the chopped actinic line, Fm' in the chopped pulse."""
    p = []
    if dark_s > 0:
        p.append(line(dark_s, f_dark, 0))
    p += [line(fs_s, f_lit, actinic), line(pulse_s, f_pulse, pulse)]
    if after_s > 0:
        p.append(line(after_s, f_lit, actinic))
    if post_dark_s > 0:
        p.append(line(post_dark_s, f_dark, 0))
    return p


def ready_after_edge_s(freq: float) -> float:
    """ADPD edge -> end of the read-out: the core sleeps or WFIs through the sequence first."""
    period_us = 1e6 / freq
    if period_us >= 1100:
        return 905e-6  # light sleep, 820-880 us actual, plus the read-out
    if period_us >= 600:
        return 510e-6  # WFI 450 us plus the read-out
    return 470e-6


def plan_period_gaps(freq: float, gap_us: float, period_us: float, write_us: float = WRITE_US):
    """Dark-gap starts (s after an ADPD edge) and dark length, as chop_period_gaps() plans them."""
    P = 1.0 / freq
    gap, per, w = gap_us * 1e-6, period_us * 1e-6, write_us * 1e-6
    w_plan = max(GAP_WRITE_PLAN_US * 1e-6, w)
    deadline = P - GAP_EDGE_GUARD_US * 1e-6
    slot = ready_after_edge_s(freq) + GAP_SETTLE_US * 1e-6
    now, starts = slot, []
    while True:
        start = max(slot, now)
        if start + w_plan + gap + w_plan > deadline:
            break
        t_dark = start + w
        starts.append(t_dark)
        now = t_dark + gap + w
        slot += per
    return np.array(starts), gap + w


def chop_plan(protocol: list[dict], gap_us=GAP_US, period_us=CHOP_PERIOD_US, write_us=WRITE_US) -> list[dict]:
    """Per line: planned gaps, dark-time fraction (light lost by the leaf), counting time."""
    out = []
    for i, s in enumerate(protocol):
        lit = s["actinic"] > ag.ACTINIC_ON
        starts, dark = plan_period_gaps(s["freq"], gap_us, period_us, write_us) if lit else (np.zeros(0), 0.0)
        n = len(starts) * max(s["pulses"] - 1, 0)
        dur = s["pulses"] / s["freq"]
        inner = max(dark - (GUARD_START_US + GUARD_END_US) * 1e-6, 0.0)
        out.append(
            dict(
                line=i,
                actinic=s["actinic"],
                freq=s["freq"],
                seconds=dur,
                lit=lit,
                gaps_per_period=len(starts),
                gaps=n,
                dark_us=dark * 1e6,
                dark_fraction=n * dark / dur if dur else 0.0,
                counting_s=n * inner if lit else dur,
            )
        )
    return out


def describe_plan(protocol, gap_us=GAP_US, period_us=CHOP_PERIOD_US) -> str:
    rows = ["line  actinic  freq   seconds  gaps/period  gaps   dark %   counting s"]
    for r in chop_plan(protocol, gap_us, period_us):
        rows.append(
            f"{r['line']:>4}  {r['actinic']:>7}  {r['freq']:>4}  {r['seconds']:>8.2f}  "
            f"{r['gaps_per_period'] if r['lit'] else '-':>11}  {r['gaps'] if r['lit'] else '-':>5}  "
            f"{100 * r['dark_fraction'] if r['lit'] else 0:>6.1f}  {r['counting_s']:>10.3f}"
        )
    return "\n".join(rows)


# ---------------------------------------------------------------------------
# AMBIT: chop verb
# ---------------------------------------------------------------------------
CHOP_RE = re.compile(r"chop:\s*(?:armed gap_us=(?P<gap>\d+) period_us=(?P<per>\d+)|(?P<off>off))")
KV_RE = re.compile(r"(\w+)=(-?\d+)")


def _console(link: ag.AmbitLink, command: str, prefix: str, wait_s: float = 1.0) -> str | None:
    """Send one console line, return the first reply line that starts with ``prefix`` or ERROR."""
    ser = link.ser
    ser.reset_input_buffer()
    ser.write((command + "\n").encode("ascii"))
    ser.flush()
    deadline = time.monotonic() + wait_s
    while time.monotonic() < deadline:
        s = ser.readline().decode("ascii", errors="replace").strip()
        if s.startswith(prefix) or s.startswith("ERROR"):
            return s
    return None


def arm_chop(link: ag.AmbitLink, gap_us: int = GAP_US, period_us: int = CHOP_PERIOD_US) -> dict:
    """Arm dark gaps for the NEXT arrun only. Raises if the firmware has no chop verb."""
    # every console argument is read up to a ',' terminator, so the line ends with a comma
    reply = _console(link, f"chop,{int(gap_us)},{int(period_us)},", "chop:")
    if reply is None:
        raise RuntimeError(
            "the AMBIT did not answer 'chop': this firmware has no chop verb. "
            "Build and flash branch feat/chop-gaps of the ambit repo."
        )
    if reply.startswith("ERROR"):
        raise ValueError(
            f"the AMBIT refused chop,{gap_us},{period_us} ({reply}). Allowed: gap 20..2000 us, period >= gap + 300 us."
        )
    m = CHOP_RE.search(reply)
    if m is None or m["off"]:
        return dict(gap_us=0, period_us=0, reply=reply)
    return dict(gap_us=int(m["gap"]), period_us=int(m["per"]), reply=reply)


def disarm_chop(link: ag.AmbitLink) -> None:
    _console(link, "chop,0,0,", "chop:")


def chop_stats(link: ag.AmbitLink) -> dict:
    """What the last run did with its gaps (firmware chopstat)."""
    reply = _console(link, "chopstat", "chopstat:")
    return {k: int(v) for k, v in KV_RE.findall(reply)} if reply else {}


# ---------------------------------------------------------------------------
# Felix: time tags
# ---------------------------------------------------------------------------
def felix_clock(felix: ag.FelixLink, integration_s: float = 1.0) -> dict:
    """Count rates plus what the TDC says about its LSB and the SYNC period."""
    rates = felix.count_rates(integration_s)
    info = dict(
        rates={int(k): float(v) for k, v in rates.items()}, f_sync_hz=float(rates.get(ag.SYNC_CH, float("nan")))
    )
    try:
        tdc = felix.backend.getTDC()[0]
        info["tdc_lsb"] = float(tdc.LSB)
        info["sync_period_lsb"] = int(tdc.ch(ag.SYNC_CH).period)
    except Exception as exc:  # informative only
        info["tdc_info_error"] = str(exc)
    return info


class TTMRecorder:
    """Felix time-tagging module writing a binary file (format: docs/binaryFormatGuide.pdf)."""

    def __init__(self, felix: ag.FelixLink, log=print):
        self.felix, self.log = felix, log
        mods = list(felix.backend.getTTMModule())
        if not mods:
            raise RuntimeError("this Felix reports no time-tagging module")
        self.ttm = mods[0]
        self.path: Path | None = None
        self.t_start = float("nan")
        self.settings: dict = {}

    def _get(self, name):
        try:
            return bool(getattr(self.ttm, name))
        except Exception:
            return None

    def start(self, path, sync: bool = True, reconstruct: bool | None = None, collapse: bool | None = None) -> None:
        tdc = self.felix.tdc
        E = tdc.TDCChannel.edgeTrigger_type
        # pylibtdc 1.21.8 exposes these read-only (reconstruction on, collapse off); set them
        # only where a library version allows it
        for name, val in (("sync_enable_reconstruction", reconstruct), ("sync_collapse", collapse)):
            if val is not None:
                try:
                    setattr(self.ttm, name, bool(val))
                except AttributeError:
                    pass
        edges = {ag.SPAD_CH: E.RISING}  # the SPAD pulse is 10 ns: one edge per photon
        if sync:
            edges[ag.SYNC_CH] = E.RISING
        ff = getattr(tdc.TTMModule, "file_format", None)  # pybind enum: never test it for truth
        fmt = ff.BINARY if ff is not None and hasattr(ff, "BINARY") else tdc.TTMModule.BINARY
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            path.unlink()
        self.ttm.enable_to_file(
            filename_base=str(path),
            file_format=fmt,
            append=False,
            max_file_count=1,
            max_file_size=0,
            unit=tdc.TTMModule.HW_LSB,
            channel_edges=edges,
            default_edges=E.NONE,
        )
        self.t_start = time.monotonic()
        self.path = path
        self.settings = dict(
            sync=bool(sync), reconstruct=self._get("sync_enable_reconstruction"), collapse=self._get("sync_collapse")
        )

    def stop(self, timeout_s: float = 60.0) -> dict:
        self.ttm.disable()
        try:
            empty = bool(self.ttm.queue_waitEmpty(float(timeout_s)))
        except Exception:
            empty = None
        seconds = time.monotonic() - self.t_start
        out = dict(seconds=seconds, queue_empty=empty, **self.settings)
        for k in ("captured_events", "lost_events"):
            try:
                out[k.split("_")[0]] = int(getattr(self.ttm, k)())
            except Exception:
                out[k.split("_")[0]] = None
        out["bytes"] = self.path.stat().st_size if self.path and self.path.exists() else 0
        if out.get("lost"):
            self.log(
                f"WARNING: the TTM lost {out['lost']:,} events (bandwidth). Drop SYNC from the "
                f"stream or lower the SPAD rate."
            )
        return out


def read_ttm(path) -> tuple[dict, np.ndarray]:
    """Header dict and a read-only memmap of the 9-byte records."""
    path = Path(path)
    with open(path, "rb") as fh:
        head = np.frombuffer(fh.read(16), "<u8")
        if head.size < 2 or int(head[0]) != TTM_MAGIC:
            raise ValueError(f"{path} is not a Felix TTM binary file")
        n_words = int(head[1])
        fh.seek(0)
        head = np.frombuffer(fh.read(8 * n_words), "<u8")
    period_fs, conv_a, conv_b = float(head[4]), int(head[5]), int(head[6])
    lsb_fs = period_fs / 2.0**conv_b
    if conv_a:
        lsb_fs *= 2.0**64 / conv_a
    n_rec = (path.stat().st_size - 8 * n_words) // REC.itemsize
    rec = np.memmap(path, dtype=REC, mode="r", offset=8 * n_words, shape=(int(n_rec),)) if n_rec else np.zeros(0, REC)
    last = int(head[8])
    info = dict(
        lsb_s=lsb_fs * 1e-15,
        n_channels=int(head[7]),
        unix_ms=int(head[2]),
        records=int(n_rec),
        last_file=last,
        lost=int(head[9]) if last else None,
    )
    return info, rec


def _chunks(rec: np.ndarray, chunk: int = CHUNK):
    for i in range(0, rec.shape[0], chunk):
        a = np.asarray(rec[i : i + chunk])
        ch = (a["ce"] & 0x7F).astype(np.int8)
        t = a["t"].astype(np.int64)
        if t.size > 1 and np.any(np.diff(t) < 0):  # channels may interleave out of order
            o = np.argsort(t, kind="stable")
            t, ch = t[o], ch[o]
        yield ch, t


# ---------------------------------------------------------------------------
# Laser micro-time
# ---------------------------------------------------------------------------
def _sync_period(ts: np.ndarray, T0: float) -> tuple[float, float, int]:
    """Sums for the mean laser period from consecutive SYNC tags (any density)."""
    if ts.size < 2:
        return 0.0, 0.0, 0
    d = np.diff(ts).astype(np.float64)
    m = np.rint(d / T0)
    ok = (m >= 1) & (m < 2e5) & (np.abs(d - m * T0) < 0.25 * T0)
    return float(d[ok].sum()), float(m[ok].sum()), int(ok.sum())


def micro_from_sync(t_ph: np.ndarray, t_sync: np.ndarray, T: float) -> np.ndarray:
    """Micro-time (same units as the tags) of each photon from the nearest preceding SYNC.

    With sparse SYNC tags the local period is taken from the bracketing pair, so a slow
    drift between the laser and the Felix clocks is interpolated, not extrapolated.
    """
    out = np.full(t_ph.size, np.nan)
    if t_sync.size == 0 or t_ph.size == 0:
        return out
    k = np.searchsorted(t_sync, t_ph, side="right") - 1
    ok = k >= 0
    k = np.clip(k, 0, t_sync.size - 1)
    dt = (t_ph - t_sync[k]).astype(np.float64)
    k2 = np.minimum(k + 1, t_sync.size - 1)
    d = (t_sync[k2] - t_sync[k]).astype(np.float64)
    m = np.rint(d / T)
    good = (m >= 1) & (np.abs(d - m * T) < 0.25 * T)
    Tl = np.where(good, d / np.maximum(m, 1), T)
    out[ok] = np.mod(dt[ok], Tl[ok])
    return out


def _chi2_score(x: np.ndarray, nb: int) -> float:
    """How far the folded phases are from uniform (all bins count, unlike a max-bin score)."""
    h = np.bincount((x * nb).astype(np.int64) % nb, minlength=nb).astype(float)
    m = h.mean()
    return float(((h - m) ** 2).sum() / max(m, 1e-12))


def _refine_period(t: np.ndarray, T0: float, log=print) -> float:
    """Epoch folding: the period that makes the folded decay sharpest, coarse to fine.

    A chunk of length L cannot tell periods apart that drift by less than about the decay
    width (~0.5 ns) over L, so every stage spans the flat region of the previous one.
    """
    T = T0
    t = t - t[0]
    # beyond ~0.3 s the laser's own phase wander (ns over seconds) dominates the folding
    # first stage +-20 ppm: enough to start from the nominal 80 MHz (this laser runs ~7 ppm off)
    for L, span, n in ((0.004, 2e-5, 801), (0.02, 2e-7, 81), (0.1, 3e-8, 41), (0.3, 6e-9, 41)):
        seg = t[t < L]
        if seg.size < 500 or t[-1] < 0.5 * L:
            break
        deltas = np.linspace(-span, span, n)
        scores = [_chi2_score((seg / (T * (1 + d))) % 1.0, 256) for d in deltas]
        T = T * (1 + deltas[int(np.argmax(scores))])
    return T


# The laser and the Felix run on separate oscillators. Measured on 3 Nov 2026 (SYNC tags,
# 1.5 s): the laser phase wanders smoothly by 2.5 ns peak to peak against a constant-period
# model, while 5 ms phase blocks interpolated linearly follow it to 7 ps rms. Tagging SYNC
# is no way out on this Felix: pylibtdc rebuilds all 80 M SYNC tags per second into the
# file (1.1 GB for 1.5 s, 14 s to write), and the collapse switch is read-only. So the
# phase is tracked from the SPAD photons themselves, all of them, lit ones included: in the
# light ~2.5 % of the ~3.6 Mcps are laser photons, i.e. ~450 per 5 ms block.
FOLD_BLOCK_S = 0.005
FOLD_NB = 512  # folded-histogram bins per period (24 ps) for the phase tracking
FOLD_HALF_WINDOW_S = 0.02  # local quadratic smoothing of the block phases
FOLD_MIN_Z = 6.0  # cross-correlation peak significance for a usable block


@dataclass
class LaserPhase:
    """micro = ((t - t_ref) / T - phase(t)) mod 1, times T; phase from 5 ms blocks."""

    T: float
    t_ref: float
    centres: np.ndarray
    phase: np.ndarray  # fraction of a period, smoothed
    valid_fraction: float
    phase_noise_ps: float  # rms of block phases around the smoothed curve
    wander_ps: float  # peak-to-peak of the smoothed phase

    def micro(self, t_s: np.ndarray) -> np.ndarray:
        x = (np.asarray(t_s, np.float64) - self.t_ref) / self.T
        return np.mod(x - np.interp(t_s, self.centres, self.phase), 1.0) * self.T


def fold_hist(t_s: np.ndarray, T: float, t_ref: float, n_blocks: int, block_s=FOLD_BLOCK_S, nb=FOLD_NB) -> np.ndarray:
    """(n_blocks x nb) folded histograms of photons in seconds; accumulate over chunks."""
    blk = np.floor(t_s / block_s).astype(np.int64)
    ok = (blk >= 0) & (blk < n_blocks)
    x = ((t_s[ok] - t_ref) / T) % 1.0
    return np.bincount(blk[ok] * nb + (x * nb).astype(np.int64) % nb, minlength=n_blocks * nb).reshape(n_blocks, nb)


def _block_shifts(h: np.ndarray, template: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Shift (fraction of a period) of each row against the template, and its significance."""
    nb = h.shape[1]
    hf = h - h.mean(1, keepdims=True)
    tf = template - template.mean()
    cc = np.fft.irfft(np.fft.rfft(hf, axis=1) * np.conj(np.fft.rfft(tf))[None, :], n=nb, axis=1)
    k = np.argmax(cc, axis=1)
    rows = np.arange(h.shape[0])
    y0, y1, y2 = cc[rows, (k - 1) % nb], cc[rows, k], cc[rows, (k + 1) % nb]
    den = y0 - 2 * y1 + y2
    frac = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
    s = (k + np.clip(frac, -0.5, 0.5)) / nb
    s = np.where(s > 0.5, s - 1.0, s)
    # significance against a flat (laser-free) block: var(cc) = mean(h) * sum(template^2)
    noise = np.sqrt(np.maximum(h.mean(1), 1e-12) * (tf**2).sum())
    return s, y1 / noise


def _smooth_phase(c: np.ndarray, s: np.ndarray, w: np.ndarray, valid: np.ndarray, hw: float) -> np.ndarray:
    """Local weighted quadratic of the block phases (linear with few neighbours)."""
    out = np.full(c.size, np.nan)
    cv, sv, wv = c[valid], s[valid], w[valid]
    if cv.size == 0:
        return np.zeros(c.size)
    lo = np.searchsorted(cv, c - hw)
    hi = np.searchsorted(cv, c + hw)
    for i in range(c.size):
        a, b = lo[i], hi[i]
        if b - a >= 5:
            p = np.polyfit(cv[a:b] - c[i], sv[a:b], 2, w=np.sqrt(wv[a:b]))
            out[i] = p[-1]
        elif b - a >= 2:
            p = np.polyfit(cv[a:b] - c[i], sv[a:b], 1, w=np.sqrt(wv[a:b]))
            out[i] = p[-1]
    ok = np.isfinite(out)
    return np.interp(c, c[ok], out[ok]) if ok.any() else np.zeros(c.size)


def track_phase(h: np.ndarray, T: float, t_ref: float, block_s=FOLD_BLOCK_S, log=print) -> LaserPhase:
    """Laser phase per block from the folded block histograms (two template iterations)."""
    nb = h.shape[1]
    n = h.shape[0]
    centres = (np.arange(n) + 0.5) * block_s
    counts = h.sum(1)
    used = counts > 0
    first = np.flatnonzero(used)
    # first template from the first 30 ms: short enough that a residual period error cannot blur it
    template = h[first[: max(1, int(0.03 / block_s))]].sum(0).astype(float)
    phase = np.zeros(n)
    for _ in range(3):
        s, z = _block_shifts(h.astype(float), template)
        valid = used & (z > FOLD_MIN_Z)
        if valid.sum() < 3:
            raise ValueError(
                "the laser peak is not visible in the photon stream: too few laser photons "
                "(need ~100 kcps; use less ND or check that the laser is on)"
            )
        su = s.copy()
        su[valid] = np.unwrap(2 * np.pi * s[valid]) / (2 * np.pi)
        phase = _smooth_phase(centres, su, counts.astype(float), valid, FOLD_HALF_WINDOW_S)
        # new template: every block rolled onto the common phase
        sh = np.rint(phase * nb).astype(np.int64)
        idx = (np.arange(nb)[None, :] + sh[:, None]) % nb
        template = np.take_along_axis(h, idx, axis=1)[used].sum(0).astype(float)
    noise = (su - phase)[valid] * T * 1e12
    return LaserPhase(
        T=T,
        t_ref=t_ref,
        centres=centres,
        phase=phase,
        valid_fraction=float(valid[used].mean()),
        phase_noise_ps=float(np.std(noise)),
        wander_ps=float(np.ptp(phase[used]) * T * 1e12),
    )


def fold_microtimes(t_s: np.ndarray, T0_s: float, log=print) -> tuple[np.ndarray, dict]:
    """Micro-times (s) of in-memory photons without SYNC tags (see LaserPhase)."""
    t = np.asarray(t_s, np.float64)
    if t.size < 500:
        raise ValueError("too few photons to fold the laser period")
    t = t - t[0]
    T = _refine_period(t, T0_s, log)
    n_blocks = int(t[-1] / FOLD_BLOCK_S) + 1
    lp = track_phase(fold_hist(t, T, 0.0, n_blocks), T, 0.0, log=log)
    return lp.micro(t), dict(
        mode="fold",
        period_s=T,
        valid_blocks=lp.valid_fraction,
        phase_noise_ps=lp.phase_noise_ps,
        wander_ps=lp.wander_ps,
    )


def rotate_to_peak(micro: np.ndarray, T: float, peak_at=PEAK_AT_NS * 1e-9) -> tuple[np.ndarray, float]:
    """Rotate micro-times (s) so the laser peak sits at ``peak_at``. Returns (micro, shift)."""
    ok = np.isfinite(micro)
    h = np.bincount((micro[ok] / T * N_MICRO).astype(np.int64) % N_MICRO, minlength=N_MICRO).astype(float)
    k = 21
    hs = np.convolve(np.concatenate([h[-k:], h, h[:k]]), np.ones(k) / k, mode="same")[k:-k]
    t_peak = (np.argmax(hs) + 0.5) * T / N_MICRO
    shift = peak_at - t_peak
    return np.mod(micro + shift, T), shift


# ---------------------------------------------------------------------------
# Pass 1: rate trace and SYNC statistics
# ---------------------------------------------------------------------------
@dataclass
class Scan:
    t0: int  # first tag (LSB); every time below is relative to it
    lsb_s: float
    duration_s: float
    counts: np.ndarray  # SPAD events per BIN_US
    n_spad: int
    n_sync: int
    sync_period_lsb: float  # mean laser period from the SYNC tags (nan without them)
    other: dict = field(default_factory=dict)
    fold: np.ndarray | None = None  # folded block histograms of every SPAD photon
    fold_T: float = float("nan")  # epoch-folded period (s) they were built with


def scan_ttm(rec: np.ndarray, lsb_s: float, T0_lsb: float, fold="auto", log=print) -> Scan:
    """Pass 1: rate trace, SYNC period and (for SPAD-only files) the folded block histograms.

    ``fold="auto"`` folds only when the file has no SYNC tags.
    """
    tmin, tmax, n_ch = None, None, np.zeros(128, np.int64)
    first_sp = None
    for ch, t in _chunks(rec):
        if t.size == 0:
            continue
        tmin = int(t[0]) if tmin is None else min(tmin, int(t[0]))
        tmax = int(t[-1]) if tmax is None else max(tmax, int(t[-1]))
        n_ch += np.bincount(ch.astype(np.int64), minlength=128)
        if first_sp is None:
            first_sp = t[ch == ag.SPAD_CH]
    if tmin is None:
        raise ValueError("the TTM file holds no events")
    per_bin = lsb_s / (BIN_US * 1e-6)
    n_bins = int((tmax - tmin) * per_bin) + 2
    duration = (tmax - tmin) * lsb_s
    counts = np.zeros(n_bins, np.int32)
    T_fold, fh, n_blocks = float("nan"), None, int(duration / FOLD_BLOCK_S) + 1
    if fold == "auto":
        fold = n_ch[ag.SYNC_CH] < 1000
    if fold and first_sp is not None and first_sp.size > 500:
        T_fold = _refine_period((first_sp - tmin) * lsb_s, T0_lsb * lsb_s, log)
        fh = np.zeros((n_blocks, FOLD_NB), np.int64)
    sd, sm, last_sync = 0.0, 0.0, None
    for ch, t in _chunks(rec):
        t = t - tmin
        sp = t[ch == ag.SPAD_CH]
        if sp.size:
            b = (sp * per_bin).astype(np.int64)
            lo = int(b.min())
            counts[lo : lo + int(b.max()) - lo + 1] += np.bincount(b - lo).astype(np.int32)
            if fh is not None:
                fh += fold_hist(sp * lsb_s, T_fold, 0.0, n_blocks)
        sy = t[ch == ag.SYNC_CH]
        if sy.size:
            if last_sync is not None:
                sy = np.concatenate([[last_sync], sy])
            d, m, _ = _sync_period(sy, T0_lsb)
            sd, sm, last_sync = sd + d, sm + m, int(sy[-1])
    others = {int(c): int(n) for c, n in enumerate(n_ch) if n and c not in (ag.SPAD_CH, ag.SYNC_CH)}
    return Scan(
        t0=tmin,
        lsb_s=lsb_s,
        duration_s=duration,
        counts=counts,
        n_spad=int(n_ch[ag.SPAD_CH]),
        n_sync=int(n_ch[ag.SYNC_CH]),
        sync_period_lsb=sd / sm if sm else float("nan"),
        other=others,
        fold=fh,
        fold_T=T_fold,
    )


# ---------------------------------------------------------------------------
# Lines and gaps on the Felix clock
# ---------------------------------------------------------------------------
def _runs(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Start and (exclusive) end indices of the True runs."""
    m = np.concatenate([[False], np.asarray(mask, bool), [False]])
    d = np.diff(m.astype(np.int8))
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def _window_sums(c: np.ndarray, w: int) -> np.ndarray:
    """S[j] = c[j:j+w].sum()."""
    C = np.concatenate([[0.0], np.cumsum(c, dtype=np.float64)])
    return C[w:] - C[:-w]


def nominal_spans(protocol: list[dict]) -> list[dict]:
    """Predicted start/end (s from run start) of each line and of its LED-on part."""
    out, t = [], 0.0
    for i, s in enumerate(protocol):
        f, n = s["freq"], s["pulses"]
        first = FIRST_LINE_WARMUP_S if i == 0 else 0.0
        t_edge0 = t + LINE_SETUP_S + first
        t_end = t_edge0 + (n - 1) / f + ready_after_edge_s(f) + LINE_END_S
        out.append(dict(t0=t, t_edge0=t_edge0, t1=t_end, lit=s["actinic"] > ag.ACTINIC_ON))
        t = t_end + 0.0003
    return out


def _otsu(v: np.ndarray) -> float:
    h, e = np.histogram(v, bins=128)
    c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h)
    w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1)
    m1 = (np.sum(h * c) - np.cumsum(h * c)) / np.maximum(w1, 1)
    return float(c[int(np.argmax(w0 * w1 * (m0 - m1) ** 2))])


def led_on_intervals(counts: np.ndarray, log=print) -> tuple[list[tuple[int, int]], dict]:
    """Intervals (bin indices) where the AMBIT LED is on, split at dark runs >= 0.6 ms.

    Chop gaps (<= ~0.5 ms) never darken a whole 0.6 ms window; line transitions (~1.3 ms),
    dark lines and the recording margins do.
    """
    k = int(round(1000 / BIN_US))
    n1 = counts.size // k
    r1 = counts[: n1 * k].reshape(n1, k).sum(1).astype(float)
    lr = np.log10(r1 + 1.0)
    thr = _otsu(lr)
    hi, lo = lr > thr, lr <= thr
    if hi.sum() == 0 or lo.sum() == 0 or 10 ** (np.median(lr[hi]) - np.median(lr[lo])) < 1.5:
        return [], dict(lit_level=float("nan"), dark_level=float(np.median(r1) / k))
    lit_level = float(np.median(r1[hi])) / k  # counts per bin while lit (gaps included)
    dark_level = float(np.median(r1[lo])) / k
    w = int(round(600 / BIN_US))
    S = _window_sums(counts.astype(float), w)
    # A 0.6 ms window holding one chop gap (<= ~0.3 ms dark) keeps >= ~55 % of the local
    # light; a line transition (~1.3 ms dark) holds only dark counts. The dark level itself
    # follows the yield (a transition right after the saturating pulse is ~2x a dark line),
    # so the threshold is local: 35 % of the brightest window within +-3 ms, and at least
    # 4 sigma above the dark-line level.
    from scipy.ndimage import maximum_filter1d

    t_dark = w * dark_level
    floor_thr = t_dark + 4 * np.sqrt(max(t_dark, 1))
    step = int(round(100 / BIN_US))
    local = maximum_filter1d(S[::step], size=2 * int(round(3000 / 100)) + 1, mode="nearest")
    thr_s = np.maximum(0.35 * np.repeat(local, step)[: S.size], floor_thr)
    dark = S < thr_s
    st, en = _runs(dark)
    # a dark run of window positions [st, en) means bins [st, en + w - 1) are dark; runs
    # shorter than 0.2 ms of positions are noise dips, not a 1.3 ms transition
    dark_bins = np.zeros(counts.size, bool)
    for a, b in zip(st, en):
        if b - a >= int(round(200 / BIN_US)):
            dark_bins[a : b + w - 1] = True
    on_st, on_en = _runs(~dark_bins)
    min_len = int(round(3000 / BIN_US))
    ivs = [(int(a), int(b)) for a, b in zip(on_st, on_en) if b - a >= min_len]
    return ivs, dict(lit_level=lit_level, dark_level=dark_level, floor_threshold=float(floor_thr))


def place_lines(counts: np.ndarray, protocol: list[dict], log=print) -> tuple[list[dict], dict]:
    """Start/end (s, Felix clock relative to the first tag) of every protocol line.

    The lit line with the highest actinic (the saturating pulse) is the anchor; every other
    line is predicted from the protocol durations and snapped to a measured LED-on interval
    when one matches.
    """
    spans = nominal_spans(protocol)
    ivs, lv = led_on_intervals(counts, log)
    to_s = BIN_US * 1e-6
    lit_idx = [i for i, s in enumerate(protocol) if s["actinic"] > ag.ACTINIC_ON]
    info = dict(levels=lv, intervals=[(a * to_s, b * to_s) for a, b in ivs], anchor=None, notes=[])
    if not lit_idx or not ivs:
        info["notes"].append("no LED-on interval found: whole recording treated as one dark line")
        dur = counts.size * to_s
        return [
            dict(line=i, t0=0.0, t1=dur if i == 0 else 0.0, matched=False, lit=False) for i in range(len(protocol))
        ], info
    led_dur = {i: spans[i]["t1"] - spans[i]["t0"] - LINE_SETUP_S + 0.005 for i in lit_idx}
    anchor = max(lit_idx, key=lambda i: (protocol[i]["actinic"], led_dur[i]))
    best, best_err = None, np.inf
    for a, b in ivs:
        err = abs((b - a) * to_s - led_dur[anchor]) / led_dur[anchor]
        if err < best_err:
            best, best_err = (a, b), err
    if best_err > 0.2:
        info["notes"].append(
            f"no LED-on interval matches the anchor line {anchor} ({led_dur[anchor]:.3f} s); positions are nominal"
        )
    # LED-on of the anchor starts ~1 ms after its line start; run start on the Felix clock
    run0 = best[0] * to_s - (spans[anchor]["t0"] + 0.001)
    info["anchor"] = anchor
    lines = []
    for i, sp in enumerate(spans):
        t0, t1, matched = run0 + sp["t0"], run0 + sp["t1"], False
        if sp["lit"]:
            pred0 = t0 + 0.001
            for a, b in ivs:
                a_s, b_s = a * to_s, b * to_s
                if abs(a_s - pred0) < 0.03 and abs((b_s - a_s) - led_dur[i]) < 0.1 * led_dur[i] + 0.01:
                    t0, t1, matched = a_s, b_s, True
                    break
        lines.append(
            dict(
                line=i,
                t0=t0,
                t1=t1,
                matched=matched,
                lit=sp["lit"],
                t_edge0=(t0 - 0.001 if matched else t0) + (sp["t_edge0"] - sp["t0"]),
            )
        )
    # dark lines end where the next lit line starts and start where the previous one ends
    for i, ln in enumerate(lines):
        if ln["lit"]:
            continue
        if i > 0 and lines[i - 1]["lit"]:
            ln["t0"] = lines[i - 1]["t1"]
        if i + 1 < len(lines) and lines[i + 1]["lit"]:
            ln["t1"] = lines[i + 1]["t0"]
    for i in lit_idx:
        if not lines[i]["matched"]:
            info["notes"].append(f"line {i}: LED-on interval not found, nominal position used")
    return lines, info


def find_gaps(
    counts: np.ndarray, i0: int, i1: int, dark_bins: float, period_bins: float, adpd_period_bins: float
) -> dict:
    """Dark gaps inside one lit line [i0, i1): one matched-filter dip per gap.

    The filter is a box of the measured dark length; its minimum marks the gap start to a
    few bins even when the gap is only ~2x darker than the light. A dip is kept only if a
    neighbour sits one chop period away (gaps of one ADPD period) or one ADPD period away
    (one gap per period); isolated dips are noise.
    """
    W = max(int(round(dark_bins)), 3)
    seg = counts[i0:i1].astype(np.float64)
    if seg.size < 4 * W:
        return dict(starts=np.zeros(0, np.int64), lit_per_bin=float("nan"), threshold=float("nan"))
    S = _window_sums(seg, W)
    # up to half of the window positions overlap a gap; the upper part of the distribution is the light
    s_lit = float(np.percentile(S, 80))
    thr = min(s_lit - 5.0 * np.sqrt(max(s_lit, 1.0)), 0.85 * s_lit)
    st, en = _runs(S < thr)
    cand = np.array([a + int(np.argmin(S[a:b])) for a, b in zip(st, en)], np.int64)
    if cand.size == 0:
        return dict(starts=cand, lit_per_bin=s_lit / W, threshold=thr)
    # one dip per half chop period: keep the deeper of two close candidates
    keep = [cand[0]]
    for c in cand[1:]:
        if c - keep[-1] < 0.5 * period_bins:
            if S[c] < S[keep[-1]]:
                keep[-1] = c
        else:
            keep.append(c)
    cand = np.array(keep, np.int64)
    tol = max(15.0 / BIN_US, 0.02 * period_bins)
    d = np.diff(cand)
    # across an ADPD edge the spacing follows the read-out, which jitters by tens of us
    near = (np.abs(d - period_bins) < tol) | (
        np.abs(d - adpd_period_bins) < max(80.0 / BIN_US, 0.02 * adpd_period_bins)
    )
    regular = np.zeros(cand.size, bool)
    regular[:-1] |= near
    regular[1:] |= near
    return dict(
        starts=cand[regular] + i0,
        lit_per_bin=s_lit / W,
        threshold=thr,
        candidates=int(cand.size),
        depth=float(np.median(S[cand[regular]]) / W) if regular.any() else float("nan"),
    )


def _fold_profile(seg: np.ndarray, P: float) -> np.ndarray:
    n = int(np.ceil(P))
    idx = np.floor(np.mod(np.arange(seg.size, dtype=np.float64), P)).astype(np.int64)
    w = np.bincount(idx, weights=seg, minlength=n)
    c = np.bincount(idx, minlength=n)
    return w / np.maximum(c, 1)


def find_gaps_folded(
    counts: np.ndarray,
    i0: int,
    i1: int,
    dark_bins: float,
    period_bins: float,
    adpd_period_bins: float,
    clock_ppm: float | None = None,
    max_per_period: int | None = None,
) -> dict:
    """Dark gaps inside one lit line [i0, i1), using the firmware's regularity.

    Inside one ADPD period the gaps follow each other at exactly the chop period (busy-wait
    timing, late by <= 3 us on the bench); only the start of each period's train moves, with
    the read-out (light-sleep wake, tens of us). So:
      1. fold the line at the ADPD period (the AMBIT clock is searched within +-200 ppm of the
         Felix clock) -> a high-statistics profile of where the gaps sit in a period;
      2. place each period's train by sliding that template by up to +-100 us;
      3. refine each gap with a box filter of the dark length within +-25 us.
    Works when single gaps are buried in shot noise (a few tens of counts per gap).
    ``clock_ppm``: AMBIT clock against the Felix, measured on a brighter line of the same run
    (then only +-6 ppm are searched). ``max_per_period``: the firmware's planned gaps per
    ADPD period; extra template dips are noise. The placed gaps are verified afterwards by
    stacking them (``verify_gaps``).
    """
    out = dict(starts=np.zeros(0, np.int64), lit_per_bin=float("nan"), template=0, periods=0)
    W = max(int(round(dark_bins)), 3)
    seg = counts[i0:i1].astype(np.float64)
    if seg.size < 3 * adpd_period_bins or seg.size < 10 * W:
        return out
    S = _window_sums(seg, W)
    s_lit = float(np.percentile(S, 80))
    out["lit_per_bin"] = s_lit / W
    k = max(1, int(round(W / 4)))
    # 1. ADPD period on the Felix clock: the one that makes the folded profile sharpest
    best_P, best_score = adpd_period_bins, -np.inf
    stages = ((200e-6, 10e-6), (12e-6, 1e-6)) if clock_ppm is None else ((6e-6, 1e-6),)
    if clock_ppm is not None:
        best_P = adpd_period_bins * (1 + clock_ppm * 1e-6)
    for span, step in stages:
        centre = best_P
        for e in np.arange(-span, span + step / 2, step):
            P = centre * (1 + e)
            prof = _fold_profile(seg, P)
            sc = float(np.var(np.convolve(prof, np.ones(k) / k, mode="valid")))
            if sc > best_score:
                best_P, best_score = P, sc
    prof = _fold_profile(seg, best_P)
    n_prof = prof.size
    n_periods = int(seg.size // best_P)
    # template: dips of the profile, one per chop slot
    ext = np.concatenate([prof, prof[:W]])
    Sp = _window_sums(ext, W)[:n_prof]
    lit_p = float(np.percentile(Sp, 80))
    noise_p = np.sqrt(max(lit_p, 1e-9) / max(n_periods, 1))
    low = Sp < min(lit_p - 6 * noise_p, 0.9 * lit_p)
    st, en = _runs(low)
    offs = []
    for a, b in zip(st, en):
        j = a + int(np.argmin(Sp[a:b]))
        if not offs or j - offs[-1] > 0.5 * period_bins:
            offs.append(j)
        elif Sp[j] < Sp[offs[-1]]:
            offs[-1] = j
    offs = np.array(offs, np.int64)
    if max_per_period and offs.size > max_per_period:  # keep the deepest, in time order
        offs = np.sort(offs[np.argsort(Sp[offs])[:max_per_period]])
    out.update(
        template=int(offs.size),
        periods=n_periods,
        adpd_period_bins=best_P,
        clock_ppm=float((best_P / adpd_period_bins - 1) * 1e6),
    )
    if offs.size == 0:
        return out
    # per-gap refinement only where a single gap stands well out of the shot noise; below
    # that the template position (period shift only) is more accurate than a noisy minimum
    gap_snr = (lit_p - float(Sp[offs].mean())) / np.sqrt(max(lit_p, 1e-9))
    refine = gap_snr > 8.0
    out.update(gap_snr=float(gap_snr), refined=bool(refine))
    # 2. each period's train: slide the template
    shifts = np.arange(-100, 101, 2) / BIN_US
    starts = []
    for p in range(n_periods + 1):  # + the partial period at the end
        base = p * best_P
        pos = np.rint(base + offs[None, :] + shifts[:, None]).astype(np.int64)
        okp = (pos >= 0) & (pos < S.size)
        tot = np.where(okp, S[np.clip(pos, 0, S.size - 1)], s_lit).sum(1)
        b = int(np.argmin(tot))
        n_in = int(okp[b].sum())
        if n_in == 0:
            continue
        depth = (s_lit * n_in - tot[b]) / np.sqrt(s_lit * n_in)
        if depth < 6.0:  # no gaps in this period (the last one, or a hiccup)
            continue
        for q in pos[b][okp[b]]:
            if not refine:
                starts.append(q)
                continue
            # 3. one gap: box filter within +-25 us
            a_, b_ = max(q - int(25 / BIN_US), 0), min(q + int(25 / BIN_US) + 1, S.size)
            starts.append(a_ + int(np.argmin(S[a_:b_])))
    out["starts"] = np.unique(np.array(starts, np.int64)) + i0
    return out


def verify_gaps(counts: np.ndarray, starts: np.ndarray, dark_bins: float) -> dict:
    """Stack the counts around every placed gap start: the dark step must sit at 0 and at
    the dark length. Misplaced trains (a wrong clock, a noise template) smear or shift it."""
    D = int(round(dark_bins))
    o = np.arange(-100, D + 101)
    idx = starts[:, None] + o[None, :]
    ok = (idx >= 0) & (idx < counts.size)
    prof = np.where(ok, counts[np.clip(idx, 0, counts.size - 1)], 0).sum(0) / np.maximum(ok.sum(0), 1)
    inside = (o >= 15) & (o <= D - 15)
    outside = ((o >= -100) & (o <= -20)) | ((o >= D + 20) & (o <= D + 100))
    p_in, p_out = float(prof[inside].mean()), float(prof[outside].mean())
    sig = np.sqrt(max(p_out, 1e-12) / max(starts.size, 1) / inside.sum())
    mid = 0.5 * (p_in + p_out)
    ps = np.convolve(prof, np.ones(3) / 3, mode="same")
    a = np.flatnonzero((o >= -50) & (o <= 50) & (ps < mid))
    b = np.flatnonzero((o >= D - 50) & (o <= D + 50) & (ps < mid))
    e0 = float(o[a[0]]) if a.size else float("nan")
    e1 = float(o[b[-1]] + 1) if b.size else float("nan")
    good = bool(p_in < p_out - 5 * sig and abs(e0) <= 8 and abs(e1 - D) <= 8)
    return dict(
        ok=good,
        inside_per_bin=p_in,
        outside_per_bin=p_out,
        edge_start_us=e0 * BIN_US,
        edge_end_us=e1 * BIN_US,
        contrast=p_out / p_in if p_in > 0 else float("inf"),
        profile=prof,
    )


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------
def run_chopped(
    ambit: ag.AmbitLink,
    felix: ag.FelixLink,
    protocol: list[dict],
    out_dir,
    gap_us: int = GAP_US,
    period_us: int = CHOP_PERIOD_US,
    pre_s: float = 0.5,
    post_s: float = 0.4,
    ttm_mode: dict | None = None,
    persist: int = 0,
    threshold_v: float = ag.THRESHOLD_V,
    laser_hz: float = ag.LASER_REP_RATE_HZ,
    log=print,
) -> Path:
    """One chopped arrun while the Felix time-tags the SPAD (and SYNC). Saves a run folder."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ttm_mode = dict(ttm_mode or dict(sync=True))
    cmd = ag.build_arrun(protocol, persist=persist)
    duration = ag.protocol_duration(protocol)
    felix.set_thresholds(threshold_v)
    clock = felix_clock(felix, 1.0)
    log("Felix input rates: " + ", ".join(f"ch{c}={r:,.0f}/s" for c, r in clock["rates"].items()))
    if clock["f_sync_hz"] < 0.9 * laser_hz:
        log("WARNING: SYNC rate below the laser repetition rate. Check the SYNC cable and laser mode.")
    armed = arm_chop(ambit, gap_us, period_us)
    log(f"AMBIT: {armed['reply']}")
    rec = TTMRecorder(felix, log=log)
    rec.start(out_dir / "ttm.bin", **ttm_mode)
    time.sleep(pre_s)
    result: dict = {}

    def worker():
        try:
            result["data"], result["t_send"], result["t_done"] = ag.run_arrun(ambit, cmd, protocol, log=log)
        except Exception as exc:  # re-raised in the main thread
            result["error"] = exc

    th = threading.Thread(target=worker, daemon=True)
    log("AMBIT: " + cmd.strip())
    t_cmd = time.monotonic()
    th.start()
    # the AMBIT dumps its arrays only after the run: stop the tags at the end of the light
    n_lines = len(protocol)
    t_stop = t_cmd + ag.NOMINAL_LATENCY_S + duration + n_lines * 0.008 + FIRST_LINE_WARMUP_S + post_s
    while time.monotonic() < t_stop and th.is_alive():
        time.sleep(0.02)
    while time.monotonic() < t_stop:
        time.sleep(0.02)
    ttm = rec.stop(timeout_s=600.0)  # with SYNC tags the file takes ~10x the run to write
    th.join(timeout=duration + 30)
    if "error" in result:
        raise result["error"]
    if th.is_alive() or "data" not in result:
        raise TimeoutError("no AMBIT data after the run")
    stats = chop_stats(ambit)
    log(f"TTM: {ttm['captured']:,} events in {ttm['seconds']:.2f} s, lost {ttm['lost']}, {ttm['bytes'] / 1e6:.0f} MB")
    if stats:
        log(
            f"AMBIT chopstat: {stats.get('gaps', 0):,} gaps in {stats.get('lines', 0)} lines, dark "
            f"{stats.get('dark_min_us', 0)}-{stats.get('dark_max_us', 0)} us, write max "
            f"{stats.get('write_max_us', 0)} us, errors {stats.get('write_errors', 0)}"
        )
        if stats.get("result", 0) != 0:
            log(f"WARNING: the AMBIT run ended with result {stats['result']}")
    save_run_meta(
        out_dir,
        protocol=protocol,
        cmd=cmd,
        data=result["data"],
        chop=armed,
        chopstat=stats,
        clock=clock,
        ttm=ttm,
        ttm_mode=ttm_mode,
        laser_hz=laser_hz,
        pre_s=pre_s,
        post_s=post_s,
        persist=persist,
        ambit=dict(name=ambit.name, fw=ambit.fw),
        felix=felix.description,
    )
    log(f"saved to {out_dir}")
    return out_dir


def preilluminate(
    ambit: ag.AmbitLink, actinic: int, seconds: float, gap_us=GAP_US, period_us=CHOP_PERIOD_US, log=print
) -> None:
    """Chopped actinic at ``actinic`` for ``seconds``; the LED stays on afterwards (persist)."""
    freq = int(max(1, min(10, (ag.MAX_POINTS - 1) // max(seconds, 1))))
    prot = [line(seconds, freq, actinic)]
    arm_chop(ambit, gap_us, period_us)
    log(f"pre-illumination: actinic {actinic} for {seconds:.0f} s (chopped)")
    ag.run_arrun(ambit, ag.build_arrun(prot, persist=1), prot, log=lambda *a: None)


def led_off(ambit: ag.AmbitLink) -> None:
    """A one-sample dark run with persist 0 switches a persistent LED off."""
    prot = [line(0.1, 10, 0)]
    ag.run_arrun(ambit, ag.build_arrun(prot), prot, log=lambda *a: None)


def save_run_meta(out_dir: Path, protocol, cmd, data, **meta) -> None:
    s = np.asarray(data.get("s_630", np.zeros(0)), float)
    r = np.asarray(data.get("r_630", np.zeros(0)), float)
    n = min(s.size, r.size)
    with open(out_dir / "ambit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["line", "k", "s_630", "r_630", "ratio", "actinic", "freq"])
        i = 0
        for li, seg in enumerate(protocol):
            for k in range(seg["pulses"]):
                if i >= n:
                    break
                wr.writerow(
                    [
                        li,
                        k,
                        int(s[i]),
                        int(r[i]),
                        f"{s[i] / r[i]:.6f}" if r[i] > 0 else "nan",
                        seg["actinic"],
                        seg["freq"],
                    ]
                )
                i += 1
    meta = dict(
        protocol=protocol,
        cmd=cmd.strip(),
        ambit_tags={k: int(len(v)) for k, v in data.items()},
        saved=time.strftime("%Y-%m-%d %H:%M:%S"),
        **meta,
    )
    (out_dir / "run.json").write_text(json.dumps(meta, indent=1, default=str), encoding="utf-8")


def load_ambit(run_dir: Path) -> dict:
    rows = list(csv.DictReader(open(Path(run_dir) / "ambit.csv")))
    return dict(
        line=np.array([int(r["line"]) for r in rows], int),
        k=np.array([int(r["k"]) for r in rows], int),
        ratio=np.array([float(r["ratio"]) for r in rows]),
        actinic=np.array([int(r["actinic"]) for r in rows]),
    )


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
def _hist_summary(h: np.ndarray, bin_ns: float) -> dict:
    """Floor, laser counts, peak FWHM and mean arrival of one rotated histogram."""
    t = (np.arange(h.size) + 0.5) * bin_ns
    peak = int(np.argmax(np.convolve(h, np.ones(5) / 5, mode="same")))
    pre = (t > 0.3) & (t < t[peak] - 0.6)
    floor = float(h[pre].mean()) if pre.sum() > 5 else 0.0
    y = np.convolve(h - floor, np.ones(5) / 5, mode="same")
    half = y[peak] / 2
    lo = peak
    while lo > 0 and y[lo] > half:
        lo -= 1
    hi = peak
    while hi < y.size - 1 and y[hi] > half:
        hi += 1
    fwhm = (hi - lo) * bin_ns * 1e3
    sel = (t >= t[peak]) & (t < t[peak] + 6.0)
    wts = np.clip(h[sel] - floor, 0, None)
    t_mean = float(((t[sel] - t[peak]) * wts).sum() / wts.sum()) if wts.sum() > 0 else float("nan")
    laser = float((h - floor).sum())
    return dict(
        floor_per_bin=floor,
        laser_counts=laser,
        counts=float(h.sum()),
        fwhm_ps=float(fwhm),
        t_mean_ns=t_mean,
        peak_ns=float(t[peak]),
    )


def _fit(h: np.ndarray, bin_s: float) -> dict:
    """fvfm_lifetime.fit_window, reduced to JSON-friendly numbers."""
    if h.sum() < 300:
        return dict(ok=False, counts=float(h.sum()))
    f = fl.fit_window(h, bin_s)
    if not f.get("ok"):
        return dict(ok=False, counts=float(h.sum()))
    return dict(
        ok=True,
        counts=f["n_counts"],
        laser_counts=f["laser_counts"],
        floor=f["floor"],
        t_mean_ns=f["t_mean"],
        tau1_ns=f["tau1"],
        tau_amp2_ns=f["tau_amp2"],
        tau_int2_ns=f["tau_int2"],
        taus2_ns=[float(x) for x in f["taus2"]],
        fractions2=[float(x) for x in f["fractions2"]],
        chi2_1=f["chi2_1"],
        chi2_2=f["chi2_2"],
        degenerate=f["fit2_degenerate"],
    )


def _ratio(a, b) -> float:
    return float(1.0 - a / b) if (b and np.isfinite(a) and np.isfinite(b) and b != 0) else float("nan")


def analyse(
    run_dir,
    fs_last_s: float = 1.0,
    fm_skip_s: float = 0.25,
    fo_last_s: float = 1.0,
    guard_start_us: float = GUARD_START_US,
    guard_end_us: float = GUARD_END_US,
    dead_time_s: float = DEAD_TIME_S,
    clock: str = "auto",
    kin_bin_s: float = 0.05,
    tig_bin_us: float = 10.0,
    early_us: float = 80.0,
    log=print,
) -> dict:
    """Decays, lifetimes, intensities and Phi2 / Fv/Fm of one chopped run folder.

    ``early_us``: besides the whole counted part of every gap, the estimates are repeated
    with only its first ``early_us`` (the ``*_early`` keys). If the yield relaxes while the
    LED is dark, the early values are the closer ones to the lit state; a clear difference
    between the two flags that relaxation.
    """
    run_dir = Path(run_dir)
    meta = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    protocol = meta["protocol"]
    chop = meta.get("chop") or {}
    gap_us = float(chop.get("gap_us") or GAP_US)
    period_us = float(chop.get("period_us") or CHOP_PERIOD_US)
    stats = meta.get("chopstat") or {}
    if stats.get("dark_min_us") and stats.get("dark_max_us"):
        dark_us = 0.5 * (stats["dark_min_us"] + stats["dark_max_us"])
        # counting stops relative to the SHORTEST dark time: at ~3 Mcps of light, 5 % of a lit
        # microsecond in the counted window is already ~150 kcps of floor
        dark_lo_us = float(stats["dark_min_us"])
    else:
        dark_us = gap_us + WRITE_US
        dark_lo_us = dark_us - 6.0
    laser_hz = float(meta.get("laser_hz", ag.LASER_REP_RATE_HZ))
    head, rec = read_ttm(run_dir / "ttm.bin")
    lsb = head["lsb_s"]
    f_sync = float((meta.get("clock") or {}).get("f_sync_hz") or laser_hz)
    T0_lsb = 1.0 / f_sync / lsb
    log(f"TTM file: {head['records']:,} records, LSB {lsb * 1e15:.2f} fs")
    scan = scan_ttm(rec, lsb, T0_lsb, fold={"auto": "auto", "sync": False, "fold": True}[clock], log=log)
    log(
        f"   {scan.duration_s:.2f} s, SPAD {scan.n_spad:,} ({scan.n_spad / scan.duration_s / 1e6:.2f} Mcps mean), "
        f"SYNC {scan.n_sync:,}"
    )
    use_sync = scan.n_sync > 1000 and clock in ("auto", "sync")
    lphase, phase_error = None, None
    if not use_sync:
        try:
            lphase = track_phase(scan.fold, scan.fold_T, 0.0, log=log)
            log(
                f"   laser phase from the photons: {lphase.valid_fraction:.0%} of 5 ms blocks usable, block noise "
                f"{lphase.phase_noise_ps:.1f} ps rms, wander {lphase.wander_ps:.0f} ps p-p"
            )
        except ValueError as exc:
            phase_error = str(exc)
            log(f"   NO MICRO-TIMES: {exc}. Lines, gaps and rates are still analysed.")

    # -- lines and gaps ------------------------------------------------------
    lines, place = place_lines(scan.counts, protocol, log)
    for n in place["notes"]:
        log("   note: " + n)
    to_bin = 1e6 / BIN_US
    sel_iv = []  # (t0, t1, line, gap id, gap start) selection intervals, seconds
    gaps_all = []
    plan = chop_plan(protocol, gap_us, period_us)
    # the AMBIT clock is measured on the brightest lit line and reused for the dimmer ones
    order = sorted(range(len(lines)), key=lambda j: -(protocol[j]["actinic"] if lines[j]["lit"] else -1))
    run_clock_ppm = None
    gap_found = {}
    for i in order:
        if not lines[i]["lit"]:
            continue
        a = lines[i]["t0"] + LINE_EDGE_S
        b = lines[i]["t1"] - LINE_EDGE_S
        i0, i1 = int(max(a, 0) * to_bin), int(max(b, 0) * to_bin)
        g = find_gaps_folded(
            scan.counts,
            i0,
            i1,
            dark_us / BIN_US,
            period_us / BIN_US,
            1e6 / protocol[i]["freq"] / BIN_US,
            clock_ppm=run_clock_ppm,
            max_per_period=plan[i]["gaps_per_period"] or None,
        )
        g["verify"] = verify_gaps(scan.counts, g["starts"], dark_us / BIN_US) if g["starts"].size else dict(ok=False)
        if run_clock_ppm is None and g["verify"]["ok"] and g.get("clock_ppm") is not None:
            run_clock_ppm = g["clock_ppm"]
        gap_found[i] = g
    for i, ln in enumerate(lines):
        seg = protocol[i]
        ln.update(
            actinic=seg["actinic"],
            freq=seg["freq"],
            pulses=seg["pulses"],
            gated=False,
            gaps=0,
            planned_gaps=plan[i]["gaps"],
        )
        a = ln["t0"] + LINE_EDGE_S
        b = ln["t1"] - LINE_EDGE_S
        if not ln["lit"]:
            if i > 0 and lines[i - 1]["lit"]:
                a += DARK_AFTER_LIT_S
            if b > a:
                sel_iv.append((a, b, i, -1, a))
            ln["mode"] = "dark"
            continue
        g = gap_found[i]
        v = g["verify"]
        ln["clock_ppm"] = g.get("clock_ppm")
        ln["gap_snr"] = g.get("gap_snr")
        ln["gap_check"] = {k: v.get(k) for k in ("ok", "contrast", "edge_start_us", "edge_end_us")}
        starts = g["starts"] / to_bin
        ln["lit_rate_cps"] = g["lit_per_bin"] * to_bin
        frac = starts.size / max(plan[i]["gaps"] * (b - a) / max(ln["t1"] - ln["t0"], 1e-9), 1)
        ln["found_fraction"] = float(frac)
        if starts.size and frac >= MIN_REGULAR_FRACTION and v.get("ok"):
            ln["gated"], ln["mode"], ln["gaps"] = True, "gaps", int(starts.size)
            g0 = starts + guard_start_us * 1e-6
            g1 = starts + (dark_lo_us - guard_end_us) * 1e-6
            for k, (u, v, s0) in enumerate(zip(g0, g1, starts)):
                sel_iv.append((u, v, i, len(gaps_all) + k, s0))
            gaps_all.extend(starts.tolist())
        else:
            ln["mode"] = "all photons"
            sel_iv.append((a, b, i, -1, a))
            if plan[i]["gaps"]:
                why = "stacked gaps not aligned" if starts.size and frac >= MIN_REGULAR_FRACTION else "too few found"
                log(
                    f"   line {i} (actinic {seg['actinic']}): {starts.size} of ~{plan[i]['gaps']} gaps ({why}), "
                    f"analysed ungated with the floor subtracted"
                )
    sel_iv.sort()
    S0 = np.array([x[0] for x in sel_iv])
    S1 = np.array([x[1] for x in sel_iv])
    S_line = np.array([x[2] for x in sel_iv], int)
    S_gap = np.array([x[3] for x in sel_iv], int)
    S_gs = np.array([x[4] for x in sel_iv])
    # counting time at bin resolution
    count_mask = np.zeros(scan.counts.size, bool)
    for u, v in zip(S0, S1):
        count_mask[int(np.ceil(u * to_bin)) : int(np.floor(v * to_bin))] = True

    # -- pass 2: selected photons and their micro-times ------------------------
    T_lsb = scan.sync_period_lsb if (use_sync and np.isfinite(scan.sync_period_lsb)) else T0_lsb
    t_sel, m_sel, l_sel, gid_sel, tig_sel = [], [], [], [], []
    carry = np.zeros(0, np.int64)
    for ch, t in _chunks(rec):
        t = t - scan.t0
        sp = t[ch == ag.SPAD_CH]
        if use_sync:
            sy = t[ch == ag.SYNC_CH]
            sy = np.concatenate([carry, sy])
            carry = sy[-2:]
        ts = sp * lsb
        j = np.searchsorted(S0, ts, side="right") - 1
        ok = j >= 0
        ok[ok] = ts[ok] < S1[j[ok]]
        if not ok.any():
            continue
        jj = j[ok]
        t_sel.append(ts[ok])
        l_sel.append(S_line[jj])
        gid_sel.append(S_gap[jj])
        tig_sel.append((ts[ok] - S_gs[jj]) * 1e6)
        if use_sync:
            m_sel.append(micro_from_sync(sp[ok], sy, T_lsb) * lsb)
        elif lphase is not None:
            m_sel.append(lphase.micro(ts[ok]))
        else:
            m_sel.append(np.full(int(ok.sum()), np.nan))
    t_sel = np.concatenate(t_sel) if t_sel else np.zeros(0)
    if t_sel.size < 1000:
        raise ValueError(f"only {t_sel.size} photons selected: check the run")
    l_sel = np.concatenate(l_sel)
    gid_sel = np.concatenate(gid_sel)
    tig_sel = np.concatenate(tig_sel)
    micro = np.concatenate(m_sel)
    if use_sync:
        T = T_lsb * lsb
        clk = dict(mode="sync", period_s=T, sync_tags=scan.n_sync)
        bad = ~np.isfinite(micro)
        if bad.mean() > 0.01:
            log(f"   {bad.mean():.1%} of the photons had no SYNC before them")
    elif lphase is not None:
        T = lphase.T
        clk = dict(
            mode="fold",
            period_s=T,
            valid_blocks=lphase.valid_fraction,
            phase_noise_ps=lphase.phase_noise_ps,
            wander_ps=lphase.wander_ps,
        )
    else:
        return _rates_only(
            run_dir,
            protocol,
            lines,
            place,
            scan,
            stats,
            chop,
            gaps_all,
            dark_us,
            S0,
            S1,
            S_line,
            S_gap,
            t_sel,
            l_sel,
            count_mask,
            meta,
            phase_error,
            guard_start_us,
            guard_end_us,
            log,
        )
    micro, shift = rotate_to_peak(micro, T)
    good = np.isfinite(micro)
    mb = np.minimum((micro[good] / T * N_MICRO).astype(np.int64), N_MICRO - 1)
    t_g, l_g, tig_g = t_sel[good], l_sel[good], tig_sel[good]
    bin_s = T / N_MICRO
    h_all = np.bincount(mb, minlength=N_MICRO).astype(float)
    clk.update(
        rotation_ns=shift * 1e9,
        **{k: v for k, v in _hist_summary(h_all, bin_s * 1e9).items() if k in ("fwhm_ps", "peak_ns")},
    )
    log(
        f"   micro-time from {clk['mode']}: period {T * 1e9:.6f} ns, peak FWHM {clk['fwhm_ps']:.0f} ps, "
        f"{t_g.size:,} photons selected"
    )

    # -- windows ---------------------------------------------------------------
    run0 = lines[0]["t0"]
    lit_idx = [i for i, ln in enumerate(lines) if ln["lit"]]
    pulse_val = max((lines[i]["actinic"] for i in lit_idx), default=0)
    pulse_lines = [i for i in lit_idx if lines[i]["actinic"] == pulse_val and pulse_val >= 150]
    windows = {}

    def add_window(name, i, a, b):
        if b - a < 0.02:
            return
        m = (l_g == i) & (t_g >= a) & (t_g < b)
        h = np.bincount(mb[m], minlength=N_MICRO).astype(float)
        cnt_s = float(count_mask[int(a * to_bin) : int(b * to_bin)].sum()) * BIN_US * 1e-6
        windows[name] = dict(
            line=i,
            t0=a - run0,
            t1=b - run0,
            counting_s=cnt_s,
            photons=int(m.sum()),
            mode=lines[i]["mode"],
            h=h,
            tig=tig_g[m] if lines[i]["gated"] else None,
        )
        if lines[i]["gated"]:
            me = m & (tig_g < guard_start_us + early_us)
            n_g = int(((S_line == i) & (S_gap >= 0) & (S_gs >= a) & (S_gs < b)).sum())
            span = min(early_us, dark_lo_us - guard_start_us - guard_end_us) * 1e-6
            windows[name]["early"] = dict(
                h=np.bincount(mb[me], minlength=N_MICRO).astype(float), counting_s=n_g * span, photons=int(me.sum())
            )

    if pulse_lines:
        p = pulse_lines[0]
        lp = lines[p]
        add_window("Fm", p, lp["t0"] + fm_skip_s, lp["t1"] - LINE_EDGE_S)
        prev = p - 1
        if prev >= 0 and lines[prev]["lit"]:
            add_window("Fs", prev, lines[prev]["t1"] - LINE_EDGE_S - fs_last_s, lines[prev]["t1"] - LINE_EDGE_S)
        first_lit = lit_idx[0]
        if first_lit > 0 and not lines[first_lit - 1]["lit"] and first_lit == p:
            d = lines[first_lit - 1]
            add_window(
                "Fo",
                first_lit - 1,
                max(d["t1"] - LINE_EDGE_S - fo_last_s, d["t0"] + LINE_EDGE_S),
                d["t1"] - LINE_EDGE_S,
            )
        nxt = p + 1
        if nxt < len(lines):
            n_ = lines[nxt]
            add_window(
                "after", nxt, n_["t0"] + LINE_EDGE_S + (0 if n_["lit"] else DARK_AFTER_LIT_S), n_["t1"] - LINE_EDGE_S
            )
    # fits, intensities
    for w in windows.values():
        w["fit"] = _fit(w["h"], bin_s)
        w["summary"] = _hist_summary(w["h"], bin_s * 1e9)
        m_rate = w["photons"] / w["counting_s"] if w["counting_s"] > 0 else float("nan")
        if lines[w["line"]]["mode"] == "all photons" and lines[w["line"]]["lit"]:
            m_rate = lines[w["line"]].get("lit_rate_cps", m_rate)  # dominated by the lit parts
        corr = 1.0 / (1.0 - m_rate * dead_time_s) if m_rate * dead_time_s < 0.9 else float("nan")
        lc = w["fit"].get("laser_counts", w["summary"]["laser_counts"])
        w["measured_cps"] = m_rate
        w["dead_time_factor"] = corr
        w["laser_cps"] = lc / w["counting_s"] * corr if w["counting_s"] > 0 else float("nan")
        w["floor_fraction"] = 1.0 - lc / max(w["photons"], 1)

    # -- AMBIT trace on the same clock -------------------------------------------
    amb = load_ambit(run_dir) if (run_dir / "ambit.csv").exists() else None
    amb_t = None
    if amb is not None and amb["ratio"].size:
        amb_t = np.array([lines[li]["t_edge0"] + k / protocol[li]["freq"] for li, k in zip(amb["line"], amb["k"])])
        for w in windows.values():
            m = (amb["line"] == w["line"]) & (amb_t - run0 >= w["t0"]) & (amb_t - run0 < w["t1"])
            w["ambit_ratio"] = float(np.nanmean(amb["ratio"][m])) if m.any() else float("nan")

    # -- estimators --------------------------------------------------------------
    est = {}

    def tau(name, key):
        return windows.get(name, {}).get("fit", {}).get(key, float("nan"))

    def amb_(name):
        return windows.get(name, {}).get("ambit_ratio", float("nan"))

    for low, label in (("Fs", "phi2"), ("Fo", "fvfm")):
        if low not in windows or "Fm" not in windows:
            continue
        e = dict(
            ambit=_ratio(amb_(low), amb_("Fm")),
            tau_2exp=_ratio(tau(low, "tau_amp2_ns"), tau("Fm", "tau_amp2_ns")),
            tau_1exp=_ratio(tau(low, "tau1_ns"), tau("Fm", "tau1_ns")),
            tau_mean=_ratio(windows[low]["summary"]["t_mean_ns"], windows["Fm"]["summary"]["t_mean_ns"]),
            intensity=_ratio(windows[low]["laser_cps"], windows["Fm"]["laser_cps"]),
        )
        if windows[low]["h"].sum() > 300 and windows["Fm"]["h"].sum() > 300:
            try:
                g = p2.global_fit(windows[low]["h"], windows["Fm"]["h"], bin_s)
                e["tau_shared"] = _ratio(g["tau_amp_s"], g["tau_amp_m"])
                e["shared_fit"] = dict(
                    taus_ns=[float(x) for x in g["taus"]],
                    tau_amp_low=g["tau_amp_s"],
                    tau_amp_fm=g["tau_amp_m"],
                    chi2r=g["chi2r"],
                )
            except Exception as exc:
                log(f"   shared fit failed: {exc}")
        # first part of each gap only (a dark Fo window has nothing to relax: used whole)
        lo_w, fm_w = windows[low], windows["Fm"]
        if "early" in fm_w:
            lo_h = lo_w["early"]["h"] if "early" in lo_w else lo_w["h"]
            fm_h = fm_w["early"]["h"]
            lo_r = (
                _hist_summary(lo_h, bin_s * 1e9)["laser_counts"]
                / lo_w["early"]["counting_s"]
                * lo_w["dead_time_factor"]
                if "early" in lo_w and lo_w["early"]["counting_s"] > 0
                else lo_w["laser_cps"]
            )
            fm_r = (
                _hist_summary(fm_h, bin_s * 1e9)["laser_counts"]
                / fm_w["early"]["counting_s"]
                * fm_w["dead_time_factor"]
                if fm_w["early"]["counting_s"] > 0
                else float("nan")
            )
            e["intensity_early"] = _ratio(lo_r, fm_r)
            # (a shared lifetime fit on a third of the photons scattered too much in simulations)
            e["t_mean_early"] = _ratio(
                _hist_summary(lo_h, bin_s * 1e9)["t_mean_ns"], _hist_summary(fm_h, bin_s * 1e9)["t_mean_ns"]
            )
        est[label] = e

    # -- kinetics and gap profile ----------------------------------------------------
    nk = int(np.ceil((lines[-1]["t1"] - run0) / kin_bin_s)) + 1
    kb = np.floor((t_g - run0) / kin_bin_s).astype(np.int64)
    okk = (kb >= 0) & (kb < nk)
    nmr = N_MICRO // REBIN
    mbr = np.minimum(mb // REBIN, nmr - 1)
    kin_h = np.bincount(kb[okk] * nmr + mbr[okk], minlength=nk * nmr).reshape(nk, nmr).astype(float)
    # counting time per kinetic bin; line 0 may start before the first tag (run0 < 0)
    cmc = np.concatenate([[0], np.cumsum(count_mask, dtype=np.int64)])
    per = kin_bin_s * to_bin
    lo_i = np.clip(np.round(run0 * to_bin + np.arange(nk) * per).astype(np.int64), 0, count_mask.size)
    hi_i = np.clip(np.round(run0 * to_bin + (np.arange(nk) + 1) * per).astype(np.int64), 0, count_mask.size)
    kin_cnt = (cmc[hi_i] - cmc[lo_i]) * BIN_US * 1e-6
    kin = _kinetics(kin_h, kin_cnt, bin_s * REBIN * 1e9, dead_time_s, kin_bin_s)
    profile = {}
    for name in ("Fs", "Fm"):
        w = windows.get(name)
        if not w or w["tig"] is None:
            continue
        m = (l_g == w["line"]) & (t_g - run0 >= w["t0"]) & (t_g - run0 < w["t1"])
        tb = np.floor(tig_g[m] / tig_bin_us).astype(np.int64)
        nt = int(np.ceil(dark_us / tig_bin_us)) + 1
        okt = (tb >= 0) & (tb < nt)
        ph = np.bincount(tb[okt] * nmr + mbr[m][okt], minlength=nt * nmr).reshape(nt, nmr).astype(float)
        # counting time per time-in-gap bin: inner windows of the gaps in this window
        gm = (S_line == w["line"]) & (S_gap >= 0) & (S_gs - run0 >= w["t0"]) & (S_gs - run0 < w["t1"])
        lo_ = (S0[gm] - S_gs[gm]) * 1e6
        hi_ = (S1[gm] - S_gs[gm]) * 1e6
        edges = np.arange(nt + 1) * tig_bin_us
        cnt = (
            np.array(
                [np.clip(np.minimum(hi_, edges[k + 1]) - np.maximum(lo_, edges[k]), 0, None).sum() for k in range(nt)]
            )
            * 1e-6
        )
        profile[name] = _kinetics(ph, cnt, bin_s * REBIN * 1e9, dead_time_s, None)
        profile[name]["t_us"] = (edges[:-1] + edges[1:]) / 2

    gap_info = dict(
        found=len(gaps_all),
        firmware=stats.get("gaps"),
        dark_us=dark_us,
        dark_lo_us=dark_lo_us,
        dark_source="chopstat" if stats.get("dark_min_us") else "nominal",
        guard_start_us=guard_start_us,
        guard_end_us=guard_end_us,
    )
    res = dict(
        run_dir=str(run_dir),
        protocol=protocol,
        lines=lines,
        place=place,
        clock=clk,
        gaps=gap_info,
        windows=windows,
        estimates=est,
        kinetics=kin,
        profile=profile,
        bin_s=bin_s,
        period_s=T,
        scan=dict(duration_s=scan.duration_s, n_spad=scan.n_spad, n_sync=scan.n_sync, other=scan.other),
        rate_1ms=_rate_trace(scan.counts, 1000),
        run0=run0,
        chop=chop,
        chopstat=stats,
        ambit=dict(
            t=None if amb_t is None else amb_t - run0,
            ratio=None if amb is None else amb["ratio"],
            line=None if amb is None else amb["line"],
        ),
        gap_starts=np.array(gaps_all) - run0,
        dark_us=dark_us,
        counts=scan.counts,
        meta=meta,
        truth=meta.get("simulated"),
    )
    log(report(res))
    return res


def _rates_only(
    run_dir,
    protocol,
    lines,
    place,
    scan,
    stats,
    chop,
    gaps_all,
    dark_us,
    S0,
    S1,
    S_line,
    S_gap,
    t_sel,
    l_sel,
    count_mask,
    meta,
    phase_error,
    guard_start_us,
    guard_end_us,
    log=print,
) -> dict:
    """Result without micro-times: line placement, gaps and photon rates in the counted parts."""
    run0 = lines[0]["t0"]
    to_bin = 1e6 / BIN_US
    for ln in lines:
        a, b = int(max(ln["t0"], 0) * to_bin), int(max(ln["t1"], 0) * to_bin)
        cnt_s = float(count_mask[a:b].sum()) * BIN_US * 1e-6
        n = int((l_sel == ln["line"]).sum())
        ln["counted_s"] = cnt_s
        ln["counted_cps"] = n / cnt_s if cnt_s > 0 else float("nan")
    res = dict(
        run_dir=str(run_dir),
        protocol=protocol,
        lines=lines,
        place=place,
        clock=dict(mode="none", error=phase_error),
        gaps=dict(
            found=len(gaps_all),
            firmware=stats.get("gaps"),
            dark_us=dark_us,
            dark_lo_us=(
                float(np.median((S1 - S0)[S_gap >= 0])) * 1e6 + guard_start_us + guard_end_us
                if (S_gap >= 0).any()
                else dark_us
            ),
            dark_source="chopstat" if stats.get("dark_min_us") else "nominal",
            guard_start_us=guard_start_us,
            guard_end_us=guard_end_us,
        ),
        windows={},
        estimates={},
        kinetics=None,
        profile={},
        run0=run0,
        chop=chop,
        chopstat=stats,
        rate_1ms=_rate_trace(scan.counts, 1000),
        gap_starts=np.array(gaps_all) - run0,
        dark_us=dark_us,
        counts=scan.counts,
        meta=meta,
        truth=meta.get("simulated"),
        ambit=dict(t=None, ratio=None, line=None),
        scan=dict(duration_s=scan.duration_s, n_spad=scan.n_spad, n_sync=scan.n_sync, other=scan.other),
    )
    log(report(res))
    return res


def _rate_trace(counts: np.ndarray, per_us: int) -> np.ndarray:
    k = int(round(per_us / BIN_US))
    n = counts.size // k
    return counts[: n * k].reshape(n, k).sum(1) / (k * BIN_US * 1e-6)


def _kinetics(h2: np.ndarray, counting_s: np.ndarray, bin_ns: float, dead_time_s: float, step_s) -> dict:
    """Laser-correlated rate and model-free mean arrival per row of a (rows x micro) matrix."""
    t = (np.arange(h2.shape[1]) + 0.5) * bin_ns
    peak = int(np.argmin(np.abs(t - PEAK_AT_NS)))
    pre = (t > 0.3) & (t < PEAK_AT_NS - 0.6)
    floor = h2[:, pre].mean(1)
    lc = h2.sum(1) - floor * h2.shape[1]
    sel = (t >= t[peak]) & (t < t[peak] + 6.0)
    w = np.clip(h2[:, sel] - floor[:, None], 0, None)
    with np.errstate(invalid="ignore", divide="ignore"):
        t_mean = (w * (t[sel] - t[peak])[None, :]).sum(1) / w.sum(1)
        m = h2.sum(1) / counting_s
        corr = 1.0 / (1.0 - m * dead_time_s)
        rate = lc / counting_s * corr
        floor_rate = floor * h2.shape[1] / counting_s * corr
    ok = (counting_s > 0) & (h2.sum(1) > 50)
    out = dict(
        laser_cps=np.where(ok, rate, np.nan),
        floor_cps=np.where(ok, floor_rate, np.nan),
        t_mean_ns=np.where(ok, t_mean, np.nan),
        counting_s=counting_s,
        counts=h2.sum(1),
    )
    if step_s:
        out["t_s"] = (np.arange(h2.shape[0]) + 0.5) * step_s
    return out


def report(res: dict) -> str:
    out = [
        f"run {res['run_dir']}",
        "",
        "line  actinic  freq   start s    end s  placed   mode          gaps  found/plan  lit kcps",
    ]
    for ln in res["lines"]:
        lit_k = f"{ln.get('lit_rate_cps', float('nan')) / 1e3:8.0f}" if ln["lit"] else "       -"
        fr = f"{ln.get('found_fraction', float('nan')):.2f}" if ln["lit"] else "-"
        out.append(
            f"{ln['line']:>4}  {ln['actinic']:>7}  {ln['freq']:>4}  {ln['t0'] - res['run0']:>7.3f}  "
            f"{ln['t1'] - res['run0']:>7.3f}  {'measured' if ln.get('matched') else 'nominal ':>8}   "
            f"{ln['mode']:<12} {ln['gaps']:>5}  {fr:>10}  {lit_k}"
        )
    g = res["gaps"]
    clk = res["clock"]
    out += [
        "",
        f"gaps: {g['found']:,} found, firmware reported {g['firmware']}, dark {g['dark_us']:.0f} us "
        f"({g['dark_source']}), counted {g.get('dark_lo_us', g['dark_us']) - g['guard_start_us'] - g['guard_end_us']:.0f} us of each",
    ]
    if clk["mode"] == "none":
        out += [f"micro-time: none ({clk.get('error')})", "", "line  counted s   photons/s in the counted parts"]
        for ln in res["lines"]:
            out.append(
                f"{ln['line']:>4}  {ln.get('counted_s', float('nan')):>9.3f}   {ln.get('counted_cps', float('nan')):>12,.0f}"
            )
        return "\n".join(out)
    out += [
        f"micro-time: {clk['mode']}, peak FWHM {clk['fwhm_ps']:.0f} ps",
        "",
        "window  line  mode          counting s   photons   laser kcps  floor %   tau1 ns  tau_amp2 ns  "
        "t_mean ns  AMBIT ratio",
    ]
    for name, w in res["windows"].items():
        f = w["fit"]
        out.append(
            f"{name:<6}  {w['line']:>4}  {w['mode']:<12}  {w['counting_s']:>10.3f}  {w['photons']:>8}  "
            f"{w['laser_cps'] / 1e3:>11.1f}  {100 * w['floor_fraction']:>7.1f}  "
            f"{f.get('tau1_ns', float('nan')):>8.3f}  {f.get('tau_amp2_ns', float('nan')):>11.3f}  "
            f"{w['summary']['t_mean_ns']:>9.3f}  {w.get('ambit_ratio', float('nan')):>11.4f}"
        )
    for label, e in res["estimates"].items():
        out += ["", f"{'Phi2' if label == 'phi2' else 'Fv/Fm'}:"]
        for k in (
            "ambit",
            "tau_shared",
            "intensity",
            "intensity_early",
            "tau_mean",
            "t_mean_early",
            "tau_2exp",
            "tau_1exp",
        ):
            if k in e:
                out.append(f"   {k:<11} {e[k]:7.3f}")
    if res.get("truth"):
        tr = res["truth"]
        out += [
            "",
            "simulation truth: " + ", ".join(f"{k} {v:.3f}" for k, v in tr.items() if isinstance(v, (int, float))),
        ]
    return "\n".join(out)


def save(res: dict, out_dir=None) -> Path:
    out_dir = Path(out_dir or res["run_dir"])
    windows = {
        k: {
            kk: (vv if kk != "early" else {a: b for a, b in vv.items() if a != "h"})
            for kk, vv in w.items()
            if kk not in ("h", "tig")
        }
        for k, w in res["windows"].items()
    }
    summary = dict(
        run_dir=res["run_dir"],
        lines=res["lines"],
        clock=res["clock"],
        gaps=res["gaps"],
        windows=windows,
        estimates=res["estimates"],
        chop=res["chop"],
        chopstat=res["chopstat"],
        place_notes=res["place"]["notes"],
        truth=res.get("truth"),
    )
    (out_dir / "chop_results.json").write_text(json.dumps(summary, indent=1, default=_json_default), encoding="utf-8")
    if res.get("kinetics") is None:
        return out_dir
    np.savez_compressed(
        out_dir / "chop_decays.npz",
        bin_s=res["bin_s"],
        **{f"h_{k}": w["h"] for k, w in res["windows"].items()},
        kin_t_s=res["kinetics"].get("t_s"),
        kin_laser_cps=res["kinetics"]["laser_cps"],
        kin_t_mean_ns=res["kinetics"]["t_mean_ns"],
        rate_1ms=res["rate_1ms"],
    )
    return out_dir


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def plot(res: dict, save_to=None, show: bool = False):
    import matplotlib.pyplot as plt

    run0 = res["run0"]
    rates_only = res.get("kinetics") is None  # no micro-times: rate trace and gaps only
    fig, axes = plt.subplots(1 if rates_only else 3, 2, figsize=(13, 4 if rates_only else 11), squeeze=False)
    ax = axes[0, 0]
    r = res["rate_1ms"]
    tt = (np.arange(r.size) + 0.5) * 1e-3 - run0
    ax.semilogy(tt, np.maximum(r, 1), lw=0.6, color="0.3", label="SPAD rate (1 ms)")
    for ln in res["lines"]:
        if ln["lit"]:
            ax.axvspan(ln["t0"] - run0, ln["t1"] - run0, color="orange", alpha=0.2 if ln["gated"] else 0.08, lw=0)
    ax.set_xlabel("time since line 0 [s]")
    ax.set_ylabel("counts / s")
    ax.set_title("SPAD rate, orange = lit lines (darker: gated)")
    am = res["ambit"]
    if am["t"] is not None:
        ax2 = ax.twinx()
        ax2.plot(am["t"], am["ratio"], ".", ms=2, color="tab:blue")
        ax2.set_ylabel("AMBIT s/r", color="tab:blue")
    # zoom into gaps of the pulse
    ax = axes[0, 1]
    w = res["windows"].get("Fm") or res["windows"].get("Fs")
    if w is None:  # rates only: middle of the brightest gated line
        gated = [ln for ln in res["lines"] if ln.get("gated")]
        if gated:
            ln = max(gated, key=lambda x: x["actinic"])
            w = dict(t0=ln["t0"] - run0, t1=ln["t1"] - run0)
    c = res["counts"]
    if w is not None and res["gap_starts"].size:
        mid = (w["t0"] + w["t1"]) / 2 + run0
        i0 = int((mid - 0.0025) * 1e6 / BIN_US)
        seg = c[i0 : i0 + int(0.005 * 1e6 / BIN_US)].astype(float)
        k = 5
        sm = np.convolve(seg, np.ones(k) / k, mode="same")
        tz = (np.arange(seg.size) * BIN_US + i0 * BIN_US) * 1e-6 - run0
        ax.plot((tz - tz[0]) * 1e3, sm / (BIN_US * 1e-6) / 1e6, lw=0.6, color="0.3")
        gs = res["gap_starts"][(res["gap_starts"] >= tz[0]) & (res["gap_starts"] < tz[-1])]
        g = res["gaps"]
        for s in gs:
            a = (s - tz[0]) * 1e3
            ax.axvspan(a, a + res["dark_us"] * 1e-3, color="tab:blue", alpha=0.15, lw=0)
            ax.axvspan(
                a + g["guard_start_us"] * 1e-3,
                a + (g.get("dark_lo_us", res["dark_us"]) - g["guard_end_us"]) * 1e-3,
                color="tab:green",
                alpha=0.3,
                lw=0,
            )
        ax.set_xlabel("ms")
        ax.set_ylabel("SPAD Mcps (5 us boxcar)")
        ax.set_title("5 ms inside the window: gaps (blue), counted part (green)")
    if rates_only:
        for a in axes.flat:
            a.grid(alpha=0.3)
        fig.tight_layout()
        if save_to:
            fig.savefig(save_to, dpi=120)
        if show:
            plt.show()
        return fig
    # kinetics
    ax = axes[1, 0]
    kin = res["kinetics"]
    ax.plot(kin["t_s"], kin["laser_cps"] / 1e3, ".-", ms=3, lw=0.7, color="tab:red")
    ax.set_ylabel("laser-correlated kcps (dead-time corr.)", color="tab:red")
    ax.set_xlabel("time since line 0 [s]")
    ax3 = ax.twinx()
    ax3.plot(kin["t_s"], kin["t_mean_ns"], ".", ms=3, color="tab:purple")
    ax3.set_ylabel("mean arrival after peak [ns]", color="tab:purple")
    for ln in res["lines"]:
        if ln["lit"]:
            ax.axvspan(ln["t0"] - run0, ln["t1"] - run0, color="orange", alpha=0.12, lw=0)
    ax.set_title("yield proxies in the counted photons (50 ms bins)")
    # decays
    ax = axes[1, 1]
    bin_ns = res["bin_s"] * 1e9
    colors = dict(Fo="tab:blue", Fs="tab:green", Fm="tab:red", after="tab:gray")
    for name, w in res["windows"].items():
        h = fl._rebin(w["h"], REBIN)
        t = (np.arange(h.size) + 0.5) * bin_ns * REBIN
        fl_ = w["summary"]["floor_per_bin"] * REBIN
        y = (h - fl_) / max((h - fl_).max(), 1)
        ax.semilogy(
            t,
            np.clip(y, 1e-4, None),
            color=colors.get(name, "k"),
            lw=0.8,
            label=f"{name}: tau_amp {w['fit'].get('tau_amp2_ns', float('nan')):.2f} ns, "
            f"floor {100 * w['floor_fraction']:.0f} %",
        )
    ax.set_xlabel("micro-time [ns]")
    ax.set_ylabel("floor-subtracted, normalised")
    ax.set_ylim(1e-3, 1.5)
    ax.legend(fontsize=8)
    ax.set_title("decays")
    # gap profiles
    for col, name in enumerate(("Fs", "Fm")):
        ax = axes[2, col]
        p = res["profile"].get(name)
        if p is None:
            ax.set_visible(False)
            continue
        ok = p["counting_s"] > 0
        ax.plot(p["t_us"][ok], p["laser_cps"][ok] / 1e3, "o-", ms=3, color="tab:red", label="laser-correlated")
        ax.plot(p["t_us"][ok], p["floor_cps"][ok] / 1e3, "s-", ms=3, color="0.5", label="flat floor")
        ax.set_xlabel("time since the LED went dark [us]")
        ax.set_ylabel("kcps")
        ax4 = ax.twinx()
        ax4.plot(p["t_us"][ok], p["t_mean_ns"][ok], "^", ms=4, color="tab:purple")
        ax4.set_ylabel("mean arrival [ns]", color="tab:purple")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_title(f"{name}: inside the gap (relaxation and delayed fluorescence check)")
    for a in axes.flat:
        a.grid(alpha=0.3)
    fig.tight_layout()
    if save_to:
        fig.savefig(save_to, dpi=120)
    if show:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Self-test: time tags against the hardware histogram (laser on, AMBIT LED off)
# ---------------------------------------------------------------------------
def _ttm_photons(path):
    """SPAD and SYNC tags (int64, relative to the first tag) and the LSB of a TTM file."""
    head, rec = read_ttm(path)
    sp, sy = [], []
    for ch, t in _chunks(rec):
        sp.append(t[ch == ag.SPAD_CH])
        sy.append(t[ch == ag.SYNC_CH])
    sp = np.concatenate(sp) if sp else np.zeros(0, np.int64)
    sy = np.concatenate(sy) if sy else np.zeros(0, np.int64)
    t0 = min([int(a[0]) for a in (sp, sy) if a.size] or [0])
    return sp - t0, sy - t0, head["lsb_s"]


def _decay_of(micro: np.ndarray, T: float) -> dict:
    micro, _ = rotate_to_peak(micro, T)
    ok = np.isfinite(micro)
    h = np.bincount(np.minimum((micro[ok] / T * N_MICRO).astype(np.int64), N_MICRO - 1), minlength=N_MICRO).astype(
        float
    )
    f = _fit(h, T / N_MICRO)
    return dict(
        h=h,
        bin_s=T / N_MICRO,
        **_hist_summary(h, T / N_MICRO * 1e9),
        tau1_ns=f.get("tau1_ns", float("nan")),
        tau_amp2_ns=f.get("tau_amp2_ns", float("nan")),
    )


def decay_from_ttm(path, f_sync_hz: float, clock: str = "auto", log=print) -> dict:
    """Decay of every SPAD photon in a TTM file: from SYNC tags if present, else folded."""
    sp, sy, lsb = _ttm_photons(path)
    T0 = 1.0 / f_sync_hz / lsb
    if sp.size < 500:
        return dict(ok=False, n_spad=int(sp.size), n_sync=int(sy.size), error="too few SPAD events")
    dur = max((sp[-1] - sp[0]) * lsb, 1e-9)
    out = dict(ok=True, n_spad=int(sp.size), n_sync=int(sy.size), spad_cps=sp.size / dur, sync_eps=sy.size / dur)
    if sy.size > 1000 and clock in ("auto", "sync"):
        d, m, _ = _sync_period(sy, T0)
        T = d / m if m else T0
        out.update(mode="sync", **_decay_of(micro_from_sync(sp, sy, T) * lsb, T * lsb))
    else:
        micro, info = fold_microtimes(sp * lsb, T0 * lsb, log=log)
        out.update(**info, **_decay_of(micro, info["period_s"]))
    return out


def fold_vs_sync(path, f_sync_hz: float, log=print) -> dict:
    """Photon by photon: micro-times folded from the SPAD stream against the SYNC truth."""
    sp, sy, lsb = _ttm_photons(path)
    if sy.size < 1000 or sp.size < 500:
        return dict(ok=False, error="needs a file with SYNC tags and SPAD photons")
    T0 = 1.0 / f_sync_hz / lsb
    d, m, _ = _sync_period(sy, T0)
    T = d / m
    Ts = T * lsb
    m_sync = micro_from_sync(sp, sy, T) * lsb
    m_fold, info = fold_microtimes(sp * lsb, T0 * lsb, log=log)
    diff = np.mod(m_fold - m_sync + Ts / 2, Ts) - Ts / 2
    ok = np.isfinite(diff)
    # judge on laser photons only: those within 0.3 ns of the SYNC-referenced peak
    ref = _decay_of(m_sync, Ts)
    mr, _ = rotate_to_peak(m_sync, Ts)
    near = ok & (np.abs(mr - PEAK_AT_NS * 1e-9) < 0.3e-9)
    off = np.median(diff[near]) if near.sum() > 50 else np.median(diff[ok])
    dd = np.mod(diff - off + Ts / 2, Ts) - Ts / 2
    core = near & (np.abs(dd) < 1e-9)
    return dict(
        ok=True,
        photons=int(sp.size),
        laser_like=int(near.sum()),
        rms_ps=float(np.std(dd[core]) * 1e12) if core.sum() > 20 else float("nan"),
        outliers=float(np.mean(np.abs(dd[near]) > 0.2e-9)) if near.sum() else float("nan"),
        fwhm_sync_ps=ref["fwhm_ps"],
        fold=info,
    )


def hardware_decay(felix: ag.FelixLink, seconds: float, laser_hz=ag.LASER_REP_RATE_HZ) -> dict:
    felix.configure_histogram(ag.SYNC_CH, ag.SPAD_CH, 0.0, 1.0 / laser_hz)
    n = max(1, int(round(seconds / 0.1)))
    felix.start_windows(0.1, n)
    felix.wait_done(seconds + 5.0)
    h = np.sum([w.bins for w in felix.windows], axis=0).astype(float)
    k = int(np.argmax(np.convolve(h, np.ones(21) / 21, mode="same")))
    h = np.roll(h, int(round(PEAK_AT_NS * 1e-9 / felix.bin_s)) - k)
    s = _hist_summary(h, felix.bin_s * 1e9)
    f = _fit(h, felix.bin_s)
    return dict(
        ok=True,
        mode="hardware",
        h=h,
        bin_s=felix.bin_s,
        **s,
        tau1_ns=f.get("tau1_ns", float("nan")),
        tau_amp2_ns=f.get("tau_amp2_ns", float("nan")),
    )


def selftest(felix: ag.FelixLink, out_dir, seconds: float = 2.0, sync_seconds: float = 0.3, log=print) -> dict:
    """Laser on, AMBIT LED off. Checks the photon-folded micro-times on the real stream.

    1. hardware histogram (SYNC-referenced) as the reference decay;
    2. ``sync_seconds`` of SPAD + SYNC tags: folded micro-times against the SYNC truth,
       photon by photon (the file is ~0.7 GB per second of SYNC, so keep it short);
    3. ``seconds`` of SPAD-only tags, the measurement mode: decay compared with (1).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    clock = felix_clock(felix, 1.0)
    rates = clock["rates"]
    log(f"Felix: SYNC {rates.get(0, 0) / 1e6:.6f} MHz, SPAD {rates.get(1, 0) / 1e3:.1f} kcps")
    res = dict(clock=clock, seconds=seconds, hardware=hardware_decay(felix, seconds))
    hw = res["hardware"]
    log(
        f"hardware histogram: laser {hw['laser_counts'] / seconds / 1e3:.1f} kcps of "
        f"{hw['counts'] / seconds / 1e3:.1f} kcps, FWHM {hw['fwhm_ps']:.0f} ps, tau1 {hw['tau1_ns']:.3f} ns"
    )
    p_sync = out_dir / "selftest_sync.bin"
    rec = TTMRecorder(felix, log=log)
    rec.start(p_sync, sync=True)
    time.sleep(sync_seconds)
    res["sync_ttm"] = rec.stop(timeout_s=300)
    try:
        res["fold_vs_sync"] = fold_vs_sync(p_sync, clock["f_sync_hz"], log=log)
    except Exception as exc:
        res["fold_vs_sync"] = dict(ok=False, error=str(exc))
    finally:
        p_sync.unlink(missing_ok=True)  # hundreds of MB of SYNC tags
    p_spad = out_dir / "selftest_spad.bin"
    rec = TTMRecorder(felix, log=log)
    rec.start(p_spad, sync=False)
    time.sleep(seconds)
    res["spad_ttm"] = rec.stop()
    try:
        res["spad_only"] = decay_from_ttm(p_spad, clock["f_sync_hz"], clock="fold", log=log)
    except Exception as exc:
        res["spad_only"] = dict(ok=False, error=str(exc))
    res["usable"] = selftest_usable(res)
    log(selftest_report(res))
    return res


def selftest_usable(res: dict) -> bool:
    hw, fv, so = res["hardware"], res.get("fold_vs_sync", {}), res.get("spad_only", {})
    return bool(so.get("ok") and so["fwhm_ps"] <= 1.3 * hw["fwhm_ps"] + 15 and (not fv.get("ok") or fv["rms_ps"] < 30))


def selftest_report(res: dict) -> str:
    hw, fv, so = res["hardware"], res.get("fold_vs_sync", {}), res.get("spad_only", {})
    sec = res.get("seconds", 1.0)
    out = [
        "                      FWHM ps   tau1 ns   t_mean ns   laser kcps",
        f"hardware histogram    {hw['fwhm_ps']:>7.0f}   {hw['tau1_ns']:>7.3f}   {hw['t_mean_ns']:>9.3f}   "
        f"{hw['laser_counts'] / sec / 1e3:>10.1f}",
    ]
    if so.get("ok"):
        dur = so["n_spad"] / so["spad_cps"]
        out.append(
            f"time tags, folded     {so['fwhm_ps']:>7.0f}   {so['tau1_ns']:>7.3f}   {so['t_mean_ns']:>9.3f}   "
            f"{so['laser_counts'] / dur / 1e3:>10.1f}"
        )
        out.append(
            f"   phase blocks usable {so.get('valid_blocks', float('nan')):.0%}, block noise "
            f"{so.get('phase_noise_ps', float('nan')):.1f} ps, wander {so.get('wander_ps', float('nan')):.0f} ps"
        )
    else:
        out.append(f"time tags, folded     failed: {so.get('error')}")
    if fv.get("ok"):
        out.append(
            f"folded vs SYNC tags: rms {fv['rms_ps']:.1f} ps on {fv['laser_like']:,} laser photons, "
            f"{100 * fv['outliers']:.1f} % off by > 200 ps"
        )
    else:
        out.append(f"folded vs SYNC tags: {fv.get('error')}")
    out.append(
        "-> folded micro-times are usable"
        if res.get("usable", selftest_usable(res))
        else "-> NOT usable: more laser photons are needed (less ND), or the laser is off"
    )
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="detect the AMBIT chop verb, the Felix and its TTM")
    ap.add_argument("--analyse", metavar="RUN_DIR", help="analyse a saved run folder")
    ap.add_argument("--simulate", metavar="OUT_DIR", help="write a synthetic run and analyse it")
    ap.add_argument("--port", help="AMBIT serial port")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)
    if args.analyse or args.simulate:
        run_dir = Path(args.analyse) if args.analyse else None
        if args.simulate:
            import chop_sim

            run_dir = chop_sim.simulate_run(Path(args.simulate), phi2_protocol(60, fs_s=1.5))
        res = analyse(run_dir)
        save(res)
        if not args.no_plot:
            plot(res, save_to=Path(run_dir) / "chop_summary.png", show=False)
        return 0
    if args.check:
        ambit = felix = None
        try:
            ambit = ag.find_ambit(args.port)
            print(arm_chop(ambit, GAP_US, CHOP_PERIOD_US)["reply"])
            disarm_chop(ambit)
            felix = ag.FelixLink().open()
            felix.set_thresholds(ag.THRESHOLD_V)
            print(felix_clock(felix, 0.5))
            TTMRecorder(felix)
            print("Felix TTM module present")
            return 0
        except Exception as exc:
            print(f"ERROR: {exc}")
            return 1
        finally:
            if felix is not None:
                felix.close()
            if ambit is not None:
                ambit.close()
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
