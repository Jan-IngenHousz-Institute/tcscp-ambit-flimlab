"""Synthetic chopped AMBIT runs in the Felix TTM binary format, for testing ambit_chop.

The model follows what was measured on the bench (3 Nov 2026) where it could be:
  * firmware gap schedule of ambit src/PAM.cpp (gaps after each read-out, every chop period,
    dark = gap + one I2C write of ~118 us), AMBIT clock -28 ppm against the Felix;
  * laser at 12.4998954 ns with a smooth phase wander of ~2.5 ns peak to peak over ~1.5 s;
  * SPAD non-paralyzable dead time 200 ns, so the LED-lit stream saturates near 3.5 Mcps;
  * LED-excited fluorescence ~40x the laser-excited signal at Fm' under the saturating pulse;
  * stray-light leak, delayed fluorescence after every LED switch-off, optional yield
    relaxation inside the gaps.
Truth (amplitude-weighted lifetimes, Phi2, Fv/Fm) is written to run.json["simulated"], and
ambit_chop.report() prints it next to the estimates.

    import chop_sim, ambit_chop as ac
    run = chop_sim.simulate_run("runs/sim_phi2", ac.phi2_protocol(100))
    res = ac.analyse(run)
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

import ambit_chop as ac
import ambit_gate as ag

LSB_S = 2.4e-9 / 2**16  # Felix: TDC period 2.4 ns / 2^16 = 36.62 fs
T_LASER_S = 12.4998954e-9  # measured laser period on the Felix clock
TAU1_NS, TAU2_NS = 0.30, 1.80  # two lifetime components; their mix sets the yield


@dataclass
class Leaf:
    tau_o: float = 0.55  # amplitude-weighted lifetimes (ns): dark-adapted Fo ...
    tau_m: float = 1.55  # ... Fm (pulse after darkness)
    tau_s: float = 0.80  # Fs in the actinic light
    tau_mp: float = 1.10  # Fm' (pulse after light)
    rise_s: float = 0.03  # yield rise in the pulse
    relax_light_s: float = 0.25
    relax_dark_s: float = 0.6
    gap_relax: float = 0.0  # fractional yield loss over a long dark gap (0 = none) ...
    gap_relax_us: float = 800.0  # ... with this time constant


@dataclass
class Optics:
    laser_fo_cps: float = 100e3  # laser-excited photons at Fo reaching the SPAD (ND OD 2)
    floor_ratio_fm: float = 40.0  # LED-excited / laser-excited at Fm' under actinic 255
    leak_cps: float = 30e3  # stray light, not laser-correlated
    dead_time_s: float = 200e-9
    df_fast: float = 3e-3  # delayed fluorescence after LED off, fraction of the lit rate
    df_fast_us: float = 10.0
    df_slow: float = 3e-4
    df_slow_us: float = 200.0
    irf_sigma_ps: float = 60.0
    irf_tail_ps: float = 50.0
    wander_ps: float = 1250.0  # amplitude of the laser phase wander (sine)
    wander_period_s: float = 3.2


def _states(protocol: list[dict], leaf: Leaf) -> list[float]:
    """Target lifetime of every line."""
    out, lit_before = [], False
    top = max(s["actinic"] for s in protocol)  # the saturating pulse: brightest line
    for s in protocol:
        a = s["actinic"]
        if a <= ag.ACTINIC_ON:
            out.append(leaf.tau_o)
        elif a == top and a >= 150:
            out.append(leaf.tau_mp if lit_before else leaf.tau_m)
        else:
            out.append(leaf.tau_s)
        lit_before = lit_before or a > ag.ACTINIC_ON
    return out


def _schedule(protocol, gap_us, period_us, t_start, clock_ppm, rng, chop=True):
    """Line spans, LED-on spans, gaps (Felix seconds), mimicking run_arr_trigger + chop."""
    k_clk = 1.0 + clock_ppm * 1e-6  # AMBIT seconds -> Felix seconds
    t, lines, gaps = t_start, [], []
    write_max = 0.0
    for i, s in enumerate(protocol):
        f, n, lit = s["freq"], s["pulses"], s["actinic"] > ag.ACTINIC_ON
        P = 1.0 / f
        t_led = t + 0.0010 if lit else None
        e0 = t + ac.LINE_SETUP_S + (ac.FIRST_LINE_WARMUP_S if i == 0 else 0.0)
        g_line = []
        if lit and chop and gap_us:
            for k in range(n - 1):
                edge = e0 + k * P
                ready = edge + ac.ready_after_edge_s(f) + (rng.uniform(-30e-6, 30e-6) if P >= 1.1e-3 else 0.0)
                deadline = edge + P - ac.GAP_EDGE_GUARD_US * 1e-6
                slot = ready + ac.GAP_SETTLE_US * 1e-6
                now = slot
                while True:
                    w_plan = min(max(ac.GAP_WRITE_PLAN_US * 1e-6, write_max), 400e-6)
                    start = max(slot, now)
                    if start + 2 * w_plan + gap_us * 1e-6 > deadline:
                        break
                    w_off, w_on = rng.normal(118e-6, 3e-6), rng.normal(118e-6, 3e-6)
                    write_max = max(write_max, w_off, w_on)
                    d0 = start + w_off
                    d1 = d0 + gap_us * 1e-6 + w_on
                    g_line.append((d0, d1))
                    now = d1
                    slot += period_us * 1e-6
        t_end = e0 + (n - 1) * P + ac.ready_after_edge_s(f) + ac.LINE_END_S
        lines.append(dict(t0=t, t_led=t_led, t1=t_end, lit=lit, actinic=s["actinic"], freq=f, e0=e0))
        gaps.append(g_line)
        t = t_end + 0.0003
    # AMBIT clock -> Felix clock, about the run start
    for ln, g in zip(lines, gaps):
        for key in ("t0", "t1", "e0", "t_led"):
            if ln[key] is not None:
                ln[key] = t_start + (ln[key] - t_start) * k_clk
        g[:] = [(t_start + (a - t_start) * k_clk, t_start + (b - t_start) * k_clk) for a, b in g]
    return lines, gaps, write_max


def _tau_at(t, lines, taus, leaf: Leaf):
    """Lifetime (ns) at times t: exponential approach to each line's target."""
    t = np.asarray(t, float)
    out = np.full(t.size, taus[0])
    prev = taus[0]
    for ln, tau in zip(lines, taus):
        m = (t >= ln["t0"]) & (t < ln["t1"] + 0.0003)
        if tau > prev:
            tc = leaf.rise_s
        else:
            tc = leaf.relax_light_s if ln["lit"] else leaf.relax_dark_s
        if ln["t0"] == lines[0]["t0"]:
            start = tau  # first line: already in its state (pre-adapted)
        else:
            start = prev
        out[m] = tau + (start - tau) * np.exp(-(t[m] - ln["t0"]) / tc)
        prev = tau + (start - tau) * np.exp(-(ln["t1"] - ln["t0"]) / tc)
    out[t >= lines[-1]["t1"]] = prev
    return out


def _deadtime(t: np.ndarray, tau: float) -> np.ndarray:
    """Indices kept by a non-paralyzable dead time (iterated to convergence)."""
    keep = np.ones(t.size, bool)
    for _ in range(8):
        tk = t[keep]
        j = np.searchsorted(tk, t, side="left") - 1
        prev = np.where(j >= 0, tk[np.clip(j, 0, None)], -np.inf)
        new = (t - prev) >= tau
        if np.array_equal(new, keep):
            break
        keep = new
    return keep


def _renewal(a: float, b: float, rate: float, dead: float, rng) -> np.ndarray:
    """Detections of a dead-time-limited SPAD under a constant rate in [a, b)."""
    if rate <= 0 or b <= a:
        return np.zeros(0)
    n = int((b - a) / (dead + 1.0 / rate) * 1.2) + 20
    t = a + np.cumsum(dead + rng.exponential(1.0 / rate, n)) - dead
    return t[t < b]


def simulate_run(
    out_dir,
    protocol: list[dict],
    gap_us: int = ac.GAP_US,
    period_us: int = ac.CHOP_PERIOD_US,
    leaf: Leaf | None = None,
    optics: Optics | None = None,
    pre_s: float = 0.3,
    post_s: float = 0.3,
    clock_ppm: float = -28.0,
    seed: int = 1,
    chop: bool = True,
    log=print,
) -> Path:
    """Write ttm.bin, ambit.csv and run.json for one synthetic chopped run."""
    t_wall = time.monotonic()
    leaf, opt = leaf or Leaf(), optics or Optics()
    rng = np.random.default_rng(seed)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t_start = pre_s + ag.NOMINAL_LATENCY_S
    lines, gaps, write_max = _schedule(protocol, gap_us, period_us, t_start, clock_ppm, rng, chop)
    taus = _states(protocol, leaf)
    t_total = lines[-1]["t1"] + post_s
    k_led = opt.floor_ratio_fm * opt.laser_fo_cps / 255.0  # LED-excited cps per mA at Fo yield

    # LED-on pieces (between gaps) and LED-off pieces
    on, off = [], []
    t_prev_off = 0.0
    for ln, g in zip(lines, gaps):
        if not ln["lit"]:
            continue
        a = ln["t_led"]
        off.append((t_prev_off, a, None))
        for g0, g1 in g:
            on.append((a, g0, ln["actinic"]))
            off.append((g0, g1, (ln["actinic"], g0)))
            a = g1
        on.append((a, ln["t1"], ln["actinic"]))
        t_prev_off = ln["t1"]
        off.append((ln["t1"], ln["t1"], (ln["actinic"], ln["t1"])))  # marks the DF after switch-off
    off.append((t_prev_off, t_total, None))

    # photons while the LED is on: saturated renewal process, laser fraction by rate
    lit_t, lit_is_laser = [], []
    mids = np.array([(a + b) / 2 for a, b, _ in on])
    y_mid = _tau_at(mids, lines, taus, leaf) / leaf.tau_o if mids.size else np.zeros(0)
    for (a, b, act), y in zip(on, y_mid):
        r_led = k_led * act * y
        r_las = opt.laser_fo_cps * y
        tt = _renewal(a, b, r_led + r_las + opt.leak_cps, opt.dead_time_s, rng)
        lit_t.append(tt)
        lit_is_laser.append(rng.random(tt.size) < r_las / (r_led + r_las + opt.leak_cps))
    # photons while the LED is off: Poisson laser + leak + delayed fluorescence, then dead time
    dk_t, dk_is_laser = [], []
    for a, b, src in off:
        if b <= a and src is None:
            continue
        dur = max(b - a, 0.0)
        tt_l = a + dur * rng.random(rng.poisson(opt.laser_fo_cps * 3.0 * dur))  # thinned below
        tt_k = a + dur * rng.random(rng.poisson(opt.leak_cps * dur))
        parts_t, parts_l = [tt_l, tt_k], [np.ones(tt_l.size, bool), np.zeros(tt_k.size, bool)]
        if src is not None:
            act, t_off = src
            y_off = _tau_at(np.array([t_off]), lines, taus, leaf)[0] / leaf.tau_o
            r_lit = k_led * act * y_off
            for frac, tc in ((opt.df_fast, opt.df_fast_us * 1e-6), (opt.df_slow, opt.df_slow_us * 1e-6)):
                span = (b - a) if b > a else 5 * tc
                n_df = rng.poisson(r_lit * frac * tc * (1 - np.exp(-span / tc)))
                u = rng.random(n_df)
                td = -tc * np.log(1 - u * (1 - np.exp(-span / tc)))
                parts_t.append(t_off + td)
                parts_l.append(np.zeros(n_df, bool))
        dk_t += parts_t
        dk_is_laser += parts_l
    dk_t = np.concatenate(dk_t)
    dk_l = np.concatenate(dk_is_laser)
    # yield loss inside a gap (QA- reoxidation in the dark): scales intensity AND lifetime
    g_s = np.array([g0 for g in gaps for g0, _ in g])
    g_e = np.array([g1 for g in gaps for _, g1 in g])

    def relax(t):
        if leaf.gap_relax <= 0 or g_s.size == 0:
            return np.ones(np.size(t))
        j = np.searchsorted(g_s, t, side="right") - 1
        jj = np.clip(j, 0, None)
        inside = (j >= 0) & (t < g_e[jj])
        tig = np.where(inside, t - g_s[jj], 0.0)
        return 1 - leaf.gap_relax * (1 - np.exp(-tig / (leaf.gap_relax_us * 1e-6)))

    # laser photons off: thin the 3x over-generated stream by the yield
    m = dk_l
    y = _tau_at(dk_t[m], lines, taus, leaf) / leaf.tau_o * relax(dk_t[m])
    keep = np.ones(dk_t.size, bool)
    keep[np.flatnonzero(m)] = rng.random(int(m.sum())) < y / 3.0
    dk_t, dk_l = dk_t[keep], dk_l[keep]
    o = np.argsort(dk_t)
    dk_t, dk_l = dk_t[o], dk_l[o]
    kd = _deadtime(dk_t, opt.dead_time_s)
    dk_t, dk_l = dk_t[kd], dk_l[kd]

    t_all = np.concatenate(lit_t + [dk_t])
    is_laser = np.concatenate(lit_is_laser + [dk_l])
    # laser photons: lock to a pulse, add IRF and the decay of the state at that time
    T = T_LASER_S
    phase = lambda t: opt.wander_ps * 1e-12 * np.sin(2 * np.pi * t / opt.wander_period_s + 0.7)
    tl = t_all[is_laser]
    k = np.floor((tl - phase(tl)) / T)
    tau_amp = _tau_at(tl, lines, taus, leaf) * relax(tl)
    f2 = np.clip((tau_amp - TAU1_NS) / (TAU2_NS - TAU1_NS), 0, 1)
    p2_ = f2 * TAU2_NS / ((1 - f2) * TAU1_NS + f2 * TAU2_NS)  # photon fraction of component 2
    tau_draw = np.where(rng.random(tl.size) < p2_, TAU2_NS, TAU1_NS) * 1e-9
    irf = rng.normal(0, opt.irf_sigma_ps * 1e-12, tl.size) + rng.exponential(opt.irf_tail_ps * 1e-12, tl.size)
    t_all[is_laser] = phase(tl) + k * T + 4e-9 + irf + rng.exponential(tau_draw)
    t_all = np.sort(t_all)

    # TTM binary file
    header = np.array([ac.TTM_MAGIC, 10, int(time.time() * 1000), 0, 2_400_000, 0, 16, 3, 1, 0], dtype="<u8")
    rec = np.zeros(t_all.size, ac.REC)
    rec["ce"] = 0x80 | ag.SPAD_CH
    rec["t"] = np.rint((t_all + 1234.5) / LSB_S).astype(np.uint64)  # board up for 20 minutes
    with open(out_dir / "ttm.bin", "wb") as fh:
        fh.write(header.tobytes())
        fh.write(rec.tobytes())

    # AMBIT trace: ratio follows the yield
    s_630, r_630 = [], []
    for ln, seg in zip(lines, protocol):
        ts = ln["e0"] + np.arange(seg["pulses"]) / seg["freq"] * (1 + clock_ppm * 1e-6)
        ratio = 0.0266 * _tau_at(ts, lines, taus, leaf) / leaf.tau_o * (1 + 0.01 * rng.standard_normal(ts.size))
        r_630 += [5300] * ts.size
        s_630 += list(np.rint(ratio * 5300).astype(int))
    data = dict(s_630=np.array(s_630), r_630=np.array(r_630))
    dark = [b - a for g in gaps for a, b in g]
    truth = dict(
        tau_o=leaf.tau_o,
        tau_m=leaf.tau_m,
        tau_s=leaf.tau_s,
        tau_mp=leaf.tau_mp,
        phi2=1 - leaf.tau_s / leaf.tau_mp,
        fvfm=1 - leaf.tau_o / leaf.tau_m,
    )
    ac.save_run_meta(
        out_dir,
        protocol=protocol,
        cmd=ag.build_arrun(protocol),
        data=data,
        chop=dict(gap_us=gap_us if chop else 0, period_us=period_us if chop else 0),
        chopstat=dict(
            result=0,
            gap_us=gap_us,
            period_us=period_us,
            gaps=len(dark),
            dark_min_us=int(min(dark) * 1e6) if dark else 0,
            dark_max_us=int(max(dark) * 1e6) if dark else 0,
            write_max_us=int(write_max * 1e6),
        ),
        clock=dict(f_sync_hz=1.0 / T_LASER_S * (1 + 3e-9), rates={0: 1 / T_LASER_S, 1: 0.0}),
        ttm=dict(captured=int(t_all.size), lost=0),
        ttm_mode=dict(sync=False),
        laser_hz=80e6,
        pre_s=pre_s,
        post_s=post_s,
        persist=0,
        ambit=dict(name="simulated", fw="sim"),
        felix="simulated",
        simulated=dict(**truth, leaf=asdict(leaf), optics=asdict(opt), clock_ppm=clock_ppm, seed=seed),
    )
    log(
        f"simulated {t_all.size:,} SPAD events over {t_total:.2f} s, {len(dark):,} gaps "
        f"({time.monotonic() - t_wall:.1f} s) -> {out_dir}"
    )
    return out_dir
