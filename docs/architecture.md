# Architecture

How the code in this repository fits together. The README covers using the set-up; this page covers
changing the code.

## Two instruments, one clock problem

The **AMBIT** is a fluorimeter on a USB-serial port (115200 8N1). It runs a protocol of up to 16
`arrun` lines (`ambit_gate.MAX_SEGMENTS`) and reports one row per measuring pulse: the 630 nm signal,
its reference, their ratio and the actinic LED current.

The **Felix** is TEDIEL's time-to-digital converter, driven through `pylibtdc`. With the laser SYNC
on channel 0 and the SPAD on channel 1 it builds hardware histograms of photon arrival times
(TCSPC): 2.34375 ps bins, one histogram per time window (5 ms in the pauses measurement).

Every analysis has to put the AMBIT's samples on the Felix's clock. The modules do it from the
measured LED edges in the Felix count rate, not from the AMBIT's nominal timing, which drifts by
about 9 ms of fixed overhead per `arrun` line (see the `ambit_pause` docstring).

## Layout

| Path | What it is |
| --- | --- |
| `scripts/ambit-gate/*.py` | The measurement and analysis library, one module per experiment |
| `scripts/ambit-gate/*.ipynb` | One notebook per experiment, numbered in the order they were developed |
| `scripts/ambit-gate/phi2gui/` | The phi2 GUI: notebook 4 (actinic pauses) as a desktop app with openJII upload |
| `flimgui/` | A separate live FLIM GUI for the Felix alone: counters, histograms, IRFs, fits |
| `irf-library/` | Measured instrument response functions (TDC-Software CSV exports) |
| `packaging/` | The PyInstaller spec, version stamping and the pylibtdc fetcher for the Windows build |
| `CAD/` | The light-tight leaf jig (FreeCAD); its own `CAD/README.md` |
| `docs/` | These pages and the vendor manuals (`docs/README.md`) |

### The library in `scripts/ambit-gate`

| Module | Role |
| --- | --- |
| `ambit_gate.py` | Base layer: finds and drives the AMBIT and the Felix, runs a protocol, aligns, saves, plots. Has a CLI (`--check`, `--help`). |
| `ambit_chop.py` | Chopped AMBIT light with Felix time tags (TTM): gap finding and lifetimes in the dark gaps |
| `chop_sim.py` | Synthetic chopped runs in the TTM binary format, to test `ambit_chop` without hardware |
| `blank_dose.py` | Blank (no leaf) and laser-dose tests; the tail fits the GUI shows |
| `ambit_pause.py` | Actinic pauses: the measurement the phi2 GUI runs (`run_pauses`, `analyse_pauses`, `run_summary`) |
| `reconv.py` | IRF reconvolution fit of a whole decay, rise included |
| `fvfm_lifetime.py`, `phi2_lifetime.py` | Fv/Fm and phi2 from the AMBIT and from the lifetimes |

The modules import each other as siblings (`import ambit_gate as ag`), with the folder on
`sys.path`: run them from `scripts/ambit-gate`, or let `phi2gui/__init__.py` insert the folder. Nothing
is installed as a package (`[tool.uv] package = false`). `ambit_gate` also puts the repository root on
`sys.path` for the optional `flimgui.fitting` import.

## The phi2 GUI

| Module | Role |
| --- | --- |
| `app.py` | `MainWindow`: instruments, openJII sign-in, the measurement form, results and figures |
| `upload.py` | One run as one openJII row (`run_row`, `write_run`), the access check and the upload loop |
| `openjii.py` | The openJII REST client (API key, experiments, uploads) and the GUI's settings file |
| `driver.py` | Finds the bundled Felix WinUSB driver installer and runs it elevated |
| `__main__.py` | Entry point; `--open RUN_FOLDER` and the hidden `--smoke-test FILE` |

### A run, end to end

1. `_run_measurement` builds the protocol (`ambit_pause.train_protocol`) and starts a `Task`.
2. `ambit_pause.run_pauses` runs the AMBIT and the Felix together and writes the run folder:
   `felix.npz` (the histogram matrix), `ambit.csv` and `run.json` (protocol, instruments, note).
3. `analyse_pauses` places the lines on the Felix clock and fits the decays with the IRF;
   `run_summary` turns that into phi2 from the AMBIT and from the lifetimes, with warnings.
4. `_show_run` writes `results.json`, draws `raw.png` and the result tabs. A failed analysis still
   shows the raw data.
5. Upload (button, or every run when "upload every run" is on): `upload.write_run` writes
   `openjii_upload.ndjson` and `upload_and_wait` sends it and polls until the platform answers;
   the outcome goes to `openjii_upload.json` in that run's folder. See [openjii-upload.md](openjii-upload.md).

Run folders live in `runs/<YYYYmmdd_HHMMSS>_gui/pauses_<label>_<HHMMSS>/`, next to the source
(`scripts/ambit-gate/runs/`, gitignored) or next to `phi2gui.exe` in the standalone build.
`phi2gui --open <run folder>` re-analyses one without instruments.

### Threads

Anything that waits on hardware or the network runs in a `Task` (a `QThread`) that calls
`fn(log)` and emits `done(result)` or `failed(text)`. Connect those signals to bound methods of
`MainWindow`, so Qt queues the call onto the GUI thread; a lambda or `functools.partial` has no
receiver object and can run on the worker thread, where touching a widget crashes. When a slot
needs more than the result, return it from `fn` as a tuple (as `_upload_current` returns the run
folder with the upload record).

`app.py` selects matplotlib's Agg backend before importing the analysis modules, which draw with
pyplot; the window embeds their figures by hand.

### Settings and user files

| What | Where |
| --- | --- |
| phi2 GUI settings, API keys per environment | `%APPDATA%\phi2gui\phi2gui.json` (`openjii.config_dir()`; `~/.config/phi2gui` on Linux) |
| IRFs added in the phi2 GUI | `<config dir>\irf-library\` |
| flimgui settings | `flimgui_settings.json` in the repository root (gitignored) |
| flimgui IRFs | `irf-library/` in the repository |

## flimgui

A live acquisition and fitting app for the Felix alone, independent of the phi2 GUI: `uv run python -m flimgui`, or `--sim` for a simulated device. Its `fitting.py` is also the
optional tail fit `ambit_gate` imports, and `tests/test_fitting.py` covers it.
