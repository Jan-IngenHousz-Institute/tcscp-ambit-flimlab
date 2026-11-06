# PSII quantum yield from fluorescence intensity and lifetime

Photosystem II quantum yield, phi2, comes from chlorophyll fluorescence measured before and during a
flash of light that saturates photosynthesis. A fluorimeter usually measures the fluorescence
intensity during a short pulse of light. Fluorescence lifetime gives the same quantity, and it does
not change when chloroplasts move or the leaf's absorbance shifts.

This repository runs an AMBIT fluorimeter and a FLIM LABS TCSPC set-up together, so one
measurement gives phi2 from both the intensity and the lifetime.

![AMBIT–TCSPC adapter](docs/ambit-tcscp-adapter.png)

## Wiring

- The Felix runs on a 5 V supply. The SPAD and the laser run on 9 V or on USB-C.
- Laser SYNC OUT 1 goes to Felix SYNC, channel 0. The laser runs at 80 MHz and 50 % intensity.
- The SPAD's SMA output goes to Felix CH1, channel 1.
- The AMBIT sits on a USB-serial adapter at 115200 8N1. The code holds DTR and RTS low, so opening
  the port does not reset the board.

## Measurement approach

The SPAD saturates during the actinic light, above about 3 Mcps. So during the saturating flash the
AMBIT stops its light for about 10 ms, and the Felix records the decay histogram in that gap.

## Get the phi2 GUI

The phi2 GUI runs the actinic-pauses measurement of notebook 4, shows its results and uploads each
run to openJII.

Download the zip for your computer from the
[latest release](https://github.com/Jan-IngenHousz-Institute/tcscp-ambit-flimlab/releases/latest),
extract all of it and read `README.txt` inside. Each zip holds everything the app needs, Python and
TEDIEL's `pylibtdc` included, so you install nothing else.

| Computer | Zip | Drives the Felix |
| --- | --- | --- |
| Windows 10 or 11, 64-bit | `phi2gui-v<version>-win64.zip` | Yes. On first use the app offers to install the Felix's WinUSB driver. |
| Linux x64 with glibc 2.35 or newer | `phi2gui-v<version>-linux-x64.zip` | Yes, after a udev rule that `README.txt` gives. Nobody has measured with it yet. |
| macOS on Apple silicon | `phi2gui-v<version>-macos-arm64.zip` | No. TEDIEL ships no macOS library, so a Mac only opens and uploads saved runs. |

The Felix itself still needs a valid TEDIEL licence. The app ships the IRF from `irf-library/` as
its default, and "Add IRF file..." takes another TDC-Software export. Every run gets a "Raw data"
tab and a `raw.png`, also when the fits fail. The app refuses a run whose Felix window is longer
than half a pause.

## Run from source

[uv](https://docs.astral.sh/uv/) manages Python 3.13 and every package. On Windows, install it with
`winget install astral-sh.uv`. Then, in the repository folder:

```
uv sync                       # Python 3.13 and every package, pylibtdc included, into .venv/
uv run nbstripout --install   # once per clone, so notebook outputs stay out of commits
```

`pyproject.toml` lists the packages and `uv.lock` pins them. `pylibtdc`, the Felix bindings, is not
on PyPI, so uv takes TEDIEL's wheels from this repository's `deps-pylibtdc-*` release. Run anything
with `uv run`, for example `uv run pytest` or `uv run python scripts/ambit-gate/ambit_gate.py`, and
pick `.venv` as the notebook kernel. `uv add <package>` adds a package and updates `uv.lock`. Commit
both files.

On Windows, double-clicking `phi2gui.bat` sets up `.venv` on first use and starts the phi2 GUI.
Without uv it falls back to an installed Python 3.13. From a terminal, run
`uv run python -m phi2gui` in `scripts/ambit-gate`.

Close TDC-Software before measuring. Only one program can open the Felix at a time.

## Upload to openJII

Uploading needs a personal API key, made in the openJII web app under Platform > Account > API
keys, and an experiment you are a member of.

Each run becomes one row of an upload table in the experiment, named `ambit_tcspc_runs` by default.
The row holds the metadata, the phi2 results and the fits as nested structures, and the AMBIT trace
and the Felix windows as nested arrays, such as `ambit_points.ratio` and `felix_windows.counts`.

Appending to an existing table needs contribute access. The first upload creates the table, and that
needs manage access, which experiment admins have. A member who is not an admin gets "there is no
upload table ... creating one needs manage access". An admin then uploads the first run under that
table name, and everyone's uploads append to it from then on. Tables made by GUI versions up to
0.1.2 held one row per AMBIT sample and per Felix window and cannot take the new rows, so the GUI
switches to the new default name once.

"Felix histograms" sets how much of each Felix window goes into the row:

| Choice | Bins | Row size for a 5 s run |
| --- | --- | --- |
| 47 ps bins, the default | 266 per window, each the sum of 20 raw bins | about 0.8 MB |
| Full resolution | every raw 2.34 ps bin, about 5300 per window | about 13 MB |
| None | results and the AMBIT trace only | about 20 kB |

`felix_windows.rebin` and `felix_windows.bin_ps` say which one a row holds. openJII stores a
full-resolution row whole. A 13 MB run uploaded to openJII dev came back from SQL with every count.
The experiment's Data page cannot list a table that holds such a row yet. It answers "This page
holds more data than one read can return". So send full-resolution runs to a table of their own,
such as `ambit_tcspc_runs_full`, and read them with SQL, a notebook or an export. A full-resolution
row above 16 MiB, from a run much longer than 5 s, goes up with 47 ps bins instead, and the log
says so.

## Files

All in `scripts/ambit-gate/`.

| File | Purpose |
| --- | --- |
| `ambit_gate.py` | Library and command line that detect the AMBIT and the Felix, run, align, analyse, save and plot |
| `1_ambit_gate.ipynb` | Step-by-step notebook that calls the library |
| `ambit_chop.py` | Chopped AMBIT light with Felix time tags, gap finding, micro-times and fits |
| `chop_sim.py` | Synthetic chopped runs in the Felix TTM format, to test `ambit_chop` without hardware |
| `2_ambit_chop.ipynb` | Step-by-step notebook for the chopped-light measurement. `SIMULATE = True` runs it offline. |
| `blank_dose.py` | Blank and laser-dose tests, with Felix histogram recordings, the AMBIT blank correction and the laser-on transient |
| `3_blank_laser_dose.ipynb` | Step-by-step notebook for the blank test without a leaf and the laser-dose test |
| `ambit_pause.py` | Actinic pauses, with ordinary arrun lines and dark pauses, 5 ms Felix histograms, SPAD recovery and decays per pause |
| `4_ambit_gate_2.ipynb` | Step-by-step notebook for the actinic-pause run |
| `reconv.py` | IRF reconvolution fit with 2 exponentials over the whole decay and its rise, Poisson likelihood, and its plot |
| `phi2gui/` | The phi2 GUI |
| `runs/<timestamp>/` | Output of each run, not in git |

## Command line

```
cd scripts/ambit-gate
uv run python ambit_gate.py --check     # detect both instruments, print count rates
uv run python ambit_gate.py             # default protocol, plot, save under runs/
uv run python ambit_gate.py --port COM5 --window 0.05 --threshold 1.0 --out C:\data\flash01
```

The default protocol lives in `DEFAULT_PROTOCOL` and in the notebook cell:

| Segment | Line | Meaning |
| --- | --- | --- |
| 1 | `2,0,0,10,0,10,0,0` | 10 samples at 10 Hz with the LED off, 1 s dark |
| 2 | `2,0,0,10,0,10,250,0` | 10 samples at 10 Hz with the LED at 250 mA, a 1 s flash |
| 3 | `2,0,0,30,0,10,0,0` | 30 samples at 10 Hz with the LED off, 3 s of relaxation |

Type 2 lines do not fire the 730 nm reference LED. In the dark phases only the 625 nm measuring
pulses, about 70 µs each, reach the SPAD, and they add a flat background too small to matter.

### Output

- `ambit.csv` holds the time since sample 0, `s_630`, `r_630`, their ratio and the actinic current.
- `felix_windows.csv` has one row per window with its phase, total counts, the model-free mean
  arrival time after the peak in `tau_mean_ns`, and a one-exponential tail fit in `tau_fit_ns`.
- `felix_windows.npz` holds the full bins matrix for your own analysis.
- `decay_pre.csv` and `decay_post.csv` hold the summed decays before and after the flash.
- `run.json` records the protocol, the alignment mode and latency, count rates and the instruments.
- `summary.png` is the four-panel figure.

## Developing

[docs/](docs/README.md) explains how the pieces fit, the openJII upload and the release pipeline.
Agents start at [AGENTS.md](AGENTS.md).
