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
| `ambit_gate.ipynb` | Step-by-step notebook calling the library |
| `runs/<timestamp>/` | Output of each run (ignored by git if you add it to `.gitignore`) |

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

## Limits

- Alignment needs the flash to be visible on the SPAD (contrast at least 3 over the dark
  windows). If your filters block it completely the script falls back to the nominal
  command latency (`NOMINAL_LATENCY_S`, 40 ms) and says so.
- Lifetimes are not measured during the flash. That needs a chopped actinic light and a
  hardware trigger into Felix CH2 (see the AMBIT gating notes in the main readme).
- The ADPD gains and currents are whatever the AMBIT has stored; set them beforehand
  with `set_currents` / `set_gains` on the console if needed.
