# ambit-gate

Fluorescence yield from the AMBIT and lifetime windows from the Felix TDC around a
saturating flash, with **software gating**: no wire between the two instruments.

## Why software gating

The AMBIT board has no free output that reaches a connector (GPIO1 drives the flash
MOSFET, the only FFC lines are BOOT, Rxd, Txd, reset_EN, GND, V+), and the Felix
has no gate input. So the Felix runs in period mode with contiguous windows
(default 100 ms) and the `arrun` command is sent while it records. The SPAD
saturates during the actinic segment, which makes the flash edges visible in the
window totals; `ambit_gate.align()` uses those edges to put the AMBIT samples on the
Felix clock (a few ms). Windows that overlap a segment edge are labelled `edge` and
excluded; `flash` windows are kept in the files but not analysed.

## Files

| File | Purpose |
| --- | --- |
| `ambit_gate.py` | Library and command line: detect AMBIT, detect Felix, run, align, analyse, save, plot |
| `1_ambit_gate.ipynb` | Step-by-step notebook calling the library |
| `ambit_chop.py` | Chopped AMBIT light + Felix time tags: run, gap finding, micro-times, fits (see below) |
| `chop_sim.py` | Synthetic chopped runs in the Felix TTM format, to test `ambit_chop` without hardware |
| `2_ambit_chop.ipynb` | Step-by-step notebook for the chopped-light measurement (`SIMULATE = True` runs it offline) |
| `blank_dose.py` | Blank and laser-dose tests: Felix histogram recordings, AMBIT blank correction, laser-on transient |
| `3_blank_laser_dose.ipynb` | Step-by-step notebook for the blank (no leaf) and the laser-dose test |
| `ambit_pause.py` | Actinic pauses: ordinary arrun lines with dark pauses, Felix 5 ms histograms, SPAD recovery and decays per pause |
| `4_ambit_gate_2.ipynb` | Step-by-step notebook for the actinic-pause run |
| `reconv.py` | IRF reconvolution fit (2 exponentials, whole decay with its rise, Poisson likelihood) and its plot |
| `runs/<timestamp>/` | Output of each run (ignored by git) |

## Requirements

Python 3.13 with `pyserial`, `numpy`, `matplotlib`, `scipy` (for the tail fit through
`flimgui/fitting.py`) and `pylibtdc` (cp313 wheel from `docs/pylibtdc-2026-09-01.zip`).
Close TDC-Software first: only one program can open the Felix.

## Wiring

- laser SYNC OUT 1 -> Felix SYNC (channel 0), laser at 80 MHz
- SPAD SMA -> Felix CH1 (channel 1)
- AMBIT on a USB-serial adapter, 115200 8N1 (DTR/RTS are held low so opening the port
  does not reset the board)

## Use

```
cd scripts/ambit-gate
python ambit_gate.py --check            # detect both, print count rates
python ambit_gate.py                    # default protocol, plot, save under runs/
python ambit_gate.py --port COM5 --window 0.05 --threshold 1.0 --out C:\data\flash01
```

Default protocol (edit `DEFAULT_PROTOCOL` or the notebook cell):

| Segment | Line | Meaning |
| --- | --- | --- |
| 1 | `2,0,0,10,0,10,0,0` | 10 samples at 10 Hz, LED off (1 s dark) |
| 2 | `2,0,0,10,0,10,250,0` | 10 samples at 10 Hz, LED at 250 mA (1 s flash) |
| 3 | `2,0,0,30,0,10,0,0` | 30 samples at 10 Hz, LED off (3 s relaxation) |

Type 2 lines do not fire the 730 nm reference LED, so only the 625 nm measuring pulses
(about 70 µs each) reach the SPAD in the dark phases; they add a flat, negligible background.

## Output

- `ambit.csv`: time since sample 0, `s_630`, `r_630`, ratio, actinic current
- `felix_windows.csv`: one row per window with phase, total counts, mean arrival time
  after the peak (`tau_mean_ns`, model free) and a one-exponential tail fit (`tau_fit_ns`)
- `felix_windows.npz`: the full bins matrix for your own analysis
- `decay_pre.csv`, `decay_post.csv`: summed decays before and after the flash
- `run.json`: protocol, alignment mode and latency, count rates, instrument identities
- `summary.png`: the four-panel figure

## Fv/Fm from the lifetime (`fvfm_lifetime.py`)

Runs the AMBIT Fv/Fm protocol (3 s dark, 0.8 s saturating pulse, 0.3 s dark) with 100 ms
Felix windows, then compares three estimates of Fv/Fm: the AMBIT ratio, the laser-correlated
Felix counts per second, and the amplitude-weighted lifetime (`1 - tau_o / tau_m`; the
fluorescence yield is proportional to that lifetime when the excitation is constant). It needs
the SPAD to stay linear during the pulse: with the laser at 100 % an NE20A (OD 2) in front of
the SPAD gives about 3 Mcps in the pulse and 50 kcps in the dark. Every window is fitted after
subtracting the flat floor of pulse-excited fluorescence, estimated from the bins before the
laser peak, with a 2-exponential tail fit weighted by the raw Poisson counts.

```
python fvfm_lifetime.py                     # measure, analyse, save fvfm.json and fvfm.png
python fvfm_lifetime.py --run runs/<dir>    # analyse a saved run
python fvfm_lifetime.py --fm-last 0.2 --fo-last 1.0 --window 0.1 --port COM13
```

Known limits: no IRF, so components near 0.2 ns are only rough; no dead-time correction on the
intensity estimate at 3 Mcps; the laser at 100 % may itself lower the AMBIT's Fv/Fm a little
(compare with the laser off).

## Phi2 in the light (`phi2_lifetime.py`)

Ten Phi2 protocols back to back (5 s at actinic X with 10 Hz sampling, 0.1 s saturating
pulse at 255, 0.5 s at X at 100 Hz). One `arrun` per repeat, sent with the persist flag so the
LED stays on while the AMBIT dumps its data; the firmware limit of 16 segments rules out a
single command. By default X steps 10, 20, ... 100 so Phi2 spans a range; `--actinic 50`
gives ten repeats at one level. The Felix records 20 ms windows per repeat; the pulse windows
are those within 10 % of the maximum total, Fs is the second before them. Lifetimes come from
a 2-exponential tail fit with the two lifetimes shared between Fs and Fm' and free amplitudes,
which keeps the amplitude-weighted lifetime stable on the 80 ms of clean pulse. Output:
`phi2.json`, `phi2.csv`, `phi2.png` with the Phi2 of both methods per repeat, a scatter with
the 1:1 line and a linear fit (slope and offset), and the dark Fv/Fm point from the latest
`*_fvfm` run for comparison.

```
python phi2_lifetime.py                      # ramp 10..100
python phi2_lifetime.py --actinic 50         # ten repeats at X = 50
python phi2_lifetime.py --run runs/<dir>     # re-analyse
```

## Limits

- Alignment needs the flash to be visible on the SPAD (contrast at least 3 over the dark
  windows). If your filters block it completely the script falls back to the nominal
  command latency (`NOMINAL_LATENCY_S`, 40 ms) and says so.
- Lifetimes are not measured during the flash in `ambit_gate.py`. The chopped-light method
  below does that, without a hardware trigger.
- The ADPD gains and currents are whatever the AMBIT has stored; set them beforehand
  with `set_currents` / `set_gains` on the console if needed.

## Chopped light: lifetimes of Fs and Fm' (`ambit_chop.py`, `2_ambit_chop.ipynb`)

During the actinic light and the saturating pulse the LED-excited fluorescence reaches the
SPAD at ~40x the laser-excited signal, and an ND filter cannot change that ratio (runs of
3 Nov 2026). So the AMBIT cuts short dark gaps into every lit line and only the SPAD photons
inside the gaps are used.

**Firmware.** Branch `feat/chop-gaps` of the ambit repo adds two console verbs:
`chop,<gap_us>,<period_us>,` arms the gaps for the next array run only (one-shot), and
`chopstat` reports what the last run did. Each gap is one AS7341 LED-register write off and
one on (~118 us each at 400 kHz), placed after the ADPD read-out and at least 100 us before
the next edge, so the AMBIT's own samples are not disturbed. The leaf is dark for
`gap_us` + one write (265-272 us at `chop,150,1000`). At 100 Hz sampling 9 gaps fit per
ADPD period, at 10 Hz 99; about 25 % of the light is lost at the defaults. On AmbitV003
(3 Nov 2026) chopped and unchopped pulses gave the same Fm (0.1002-0.1015 against
0.1003-0.1024).

**Felix.** Only the SPAD is time-tagged (TTM, binary file, ~30 MB per second of light).
Tagging SYNC is impractical: pylibtdc rebuilds all 80 M laser pulses per second into the file
(0.7 GB per second, written at a tenth of real time), and its collapse switch is read-only.
The laser micro-time is instead rebuilt from the photons: the period is found by epoch
folding and the phase is tracked in 5 ms blocks. The laser wanders ~2.5 ns per second against
the Felix clock (measured from SYNC tags); 5 ms blocks follow that to a few ps. This needs
~50-100 kcps of laser-correlated photons: about OD 2, not OD 4 (`ac.selftest` checks it
against the hardware histogram and photon by photon against SYNC tags).

**Analysis.** `ac.analyse(run_dir)`:
1. places the lines on the Felix clock from the LED-on intervals (the brightest line anchors);
2. finds the gaps of each lit line by folding it at the ADPD period (the AMBIT clock is
   measured, ~-30..-50 ppm), sliding that template per period, and checks the result by
   stacking all gaps (the dark step must sit at 0 and at the measured dark length);
   a line whose gaps fail the check is analysed with all its photons and the floor subtracted;
3. counts photons from 10 us after the LED goes dark to 8 us before the shortest dark time;
4. builds Fo / Fs / Fm' decays and fits them with `fvfm_lifetime.fit_window` and the shared
   2-exponential `phi2_lifetime.global_fit`.

In simulations with the bench numbers (single run: 3 s Fs, 0.8 s pulse, OD 2, true Phi2 0.273)
`intensity` gave 0.270-0.275 (sd 0.015) and `tau_shared` 0.28-0.30 (sd 0.03-0.04); the
independent 2-exponential fits scatter by 0.1. The plots show the yield inside the gap: a fall
of the laser-correlated rate with time since the LED went dark would mean the yield relaxes
in the gaps (`intensity_early` uses only the first 80 us).

```
python ambit_chop.py --check                 # chop verb, Felix, TTM module
python ambit_chop.py --analyse runs/<dir>    # re-analyse a saved run
python ambit_chop.py --simulate runs/sim     # synthetic run through the same analysis
```

## Blank and laser dose (`blank_dose.py`, `3_blank_laser_dose.ipynb`)

Two checks of what the instruments see, with the Felix in histogram mode (micro-times from the
hardware against SYNC: no chopping, gating or reconstruction).

**Blank.** Felix: 10 s of histograms with the leaf, without it, and without it with the laser
blocked; `blank_felix` gives the laser-correlated rate (above the flat floor, from 0.5 ns before
the leaf's peak), the floor, and each recording's own peak, FWHM, mean arrival and tail lifetime,
plus the leaf with the blank's laser-correlated part removed. AMBIT: the leaf runs' chopped
protocols without a leaf; `ambit_blank_correct` subtracts the blank s/r window by window (Fo or
Fs, Fm', after) and the no-leaf Felix rate from the leaf run's gap intensities.

**Laser dose.** A leaf dark-adapted with the laser blocked; the AMBIT runs measuring pulses at
20 Hz, the cell cues the user to unblock the laser after `dark_s`, then a saturating pulse.
`analyse_dose` finds the laser-on moment from the Felix (5-window step search against the
laser-blocked noise, about one 20 ms window), and gives the AMBIT yield before and after laser-on
(true dark Fo, step, rise, Fv/Fm), the Felix rate and mean arrival in log-spaced bins after
laser-on, and the decays of the first ~5000 laser photons and of the last 5 s. Mean arrival times
are floor-subtracted without clipping (clipping the tail noise biases them up at low laser
fractions) and come with their uncertainty.

## Actinic pauses (`ambit_pause.py`, `4_ambit_gate_2.ipynb`)

What the SPAD counts when the AMBIT light stops, without chop firmware or time tags. The AMBIT
runs ordinary `arrun` lines at 50 Hz: per actinic level (0, 100, 250), 1 s of light, a pause of
2 samples, 1 s, 3 samples, 1 s, 5 samples. Every line carries about 9 ms of fixed time, so the
pauses are about 29, 49 and 89 ms dark. Actinic 0 is one dark line, because the AMBIT takes at
most 16 lines. The Felix records 5 ms hardware histograms throughout. The lines are placed on the
Felix clock from the measured LED edges (`led_intervals`, `place`): the nominal `ag.time_axis`
ignores the per-line overheads, which add up to about 0.1 s over 13 lines.

`analyse_pauses` gives the reference (the actinic-0 line: laser plus measuring pulses), the SPAD
rate while lit, and per pause the first fully dark window and the summed decay against the
reference. `plot_pauses` draws the overview, a zoom on every pause against the time since
LED-off (SPAD total and AMBIT s/r), and every pause (Fm n) against the decay before the run (Fo)
with both fits and phi2_lifetime = (tau2,Fm - tau2,Fo) / tau2,Fm.

Many actinic levels: one `arrun` takes at most 16 lines, so `ap.run_levels` runs one `arrun` per
lit level (the actinic-0 reference line, then the level's lit blocks and pauses: 7 lines), each
saved as its own run; `ap.levels_report` puts all levels in one table (reference, lit block, every
pause: tau2, mean arrival, fluorescence relative yield).

Pulse train: `ap.train_protocol(level, dark_points, block_s, pause_samples, freq)` builds dark
(dark_points) | lit block_s, pause | ... | dark (dark_points), one (lit, pause) pair per entry of
pause_samples (2 + 2 x pauses lines, so up to 7 pauses in 16 lines); the first dark line is the
reference. Run it with `run_pauses` like any other protocol. Rates are raw counts per
second.

Lifetimes (`reconv.py`, used by `ap.fit_report` and `ap.plot_fits`): every decay is fitted as a
whole, rise included, with two exponentials convolved with the measured IRF (periodic over the laser
period), the floor and an IRF shift free, from 1.5 ns before the laser peak to 7 ns after it (the
echo ~8 ns after the peak is outside). Reported per component: lifetime, peak (counts per 47 ps bin
at the decay peak, as drawn: peaks + floor = data peak) and share of the photons; the amplitude
before the IRF only for a resolved component (tau >= half the IRF FWHM), since for a shorter one it
is not measurable. t_avg = amplitude-weighted lifetime of the resolved components; t_int =
photon-weighted mean lifetime of all components (the model counterpart of t_mean, which also contains the IRF delay after its peak). Tested on synthetic decays
made with the same IRF: unbiased, errors from the covariance match the scatter.
