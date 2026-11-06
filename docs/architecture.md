# Architecture

How the code fits together. The README covers using the set-up, and this page covers changing the
code.

## Two instruments, one clock

The AMBIT is a fluorimeter on a USB-serial port at 115200 8N1. It runs a protocol of up to 16
`arrun` lines, set by `ambit_gate.MAX_SEGMENTS`, and reports one row per measuring pulse with the
630 nm signal, its reference, their ratio and the actinic LED current.

The Felix is TEDIEL's time-to-digital converter, driven through `pylibtdc`. With the laser SYNC on
channel 0 and the SPAD on channel 1, it builds hardware histograms of photon arrival times. Each
bin is 2.34375 ps, and the Felix makes one histogram per time window, 5 ms in the pauses
measurement.

Every analysis has to put the AMBIT's samples on the Felix's clock. The modules place them from the
LED edges they measure in the Felix count rate. The AMBIT's nominal timing would be off by about
9 ms of fixed overhead per `arrun` line, as the `ambit_pause` docstring explains.

## Layout

| Path | What it is |
| --- | --- |
| `scripts/ambit-gate/*.py` | The measurement and analysis library, one module per experiment |
| `scripts/ambit-gate/*.ipynb` | One notebook per experiment, numbered in the order they were written |
| `scripts/ambit-gate/phi2gui/` | The phi2 GUI, which runs notebook 4's actinic pauses as a desktop app with openJII upload |
| `flimgui/` | A separate live FLIM GUI for the Felix alone, with counters, histograms, IRFs and fits |
| `irf-library/` | Measured instrument response functions, exported from TDC-Software as CSV |
| `packaging/` | The PyInstaller spec, the bundle check and the version stamp for the standalone builds |
| `CAD/` | The light-tight leaf jig in FreeCAD, described in `CAD/README.md` |
| `docs/` | These pages and the vendor manuals, listed in `docs/README.md` |

### The library in `scripts/ambit-gate`

| Module | Role |
| --- | --- |
| `ambit_gate.py` | The base layer. It finds and drives the AMBIT and the Felix, runs a protocol, aligns, saves and plots, and has a CLI. |
| `ambit_chop.py` | Chopped AMBIT light with Felix time tags, gap finding and lifetimes in the dark gaps |
| `chop_sim.py` | Synthetic chopped runs in the TTM binary format, to test `ambit_chop` without hardware |
| `blank_dose.py` | Blank tests without a leaf, laser-dose tests, and the tail fits the GUI shows |
| `ambit_pause.py` | Actinic pauses, the measurement the phi2 GUI runs, in `run_pauses`, `analyse_pauses` and `run_summary` |
| `reconv.py` | IRF reconvolution fit of a whole decay, its rise included |
| `fvfm_lifetime.py`, `phi2_lifetime.py` | Fv/Fm and phi2 from the AMBIT and from the lifetimes |

The modules import each other as siblings, as in `import ambit_gate as ag`, with the folder on
`sys.path`. Run them from `scripts/ambit-gate`, or let `phi2gui/__init__.py` add the folder to the
path. uv does not install the repository as a package, because `pyproject.toml` sets `package = false`.
`ambit_gate` also puts the repository root on `sys.path` for its optional `flimgui.fitting` import.

## The phi2 GUI

| Module | Role |
| --- | --- |
| `app.py` | `MainWindow`, with the instruments, the openJII sign-in, the measurement form, the results and the figures |
| `upload.py` | Turns a run into one openJII row in `run_row` and `write_run`, checks access, and runs the upload |
| `openjii.py` | The openJII REST client for the API key, experiments and uploads, and the GUI's settings file |
| `driver.py` | Finds the bundled Felix WinUSB driver installer and runs it with administrator rights, on Windows |
| `__main__.py` | The entry point, with `--open RUN_FOLDER` and the hidden `--smoke-test FILE` |

### A run, end to end

1. `_run_measurement` builds the protocol with `ambit_pause.train_protocol` and starts a `Task`.
2. `ambit_pause.run_pauses` runs the AMBIT and the Felix together and writes the run folder. It
   saves the histogram matrix in `felix.npz`, the AMBIT rows in `ambit.csv`, and the protocol,
   instruments and note in `run.json`.
3. `analyse_pauses` places the lines on the Felix clock and fits the decays with the IRF.
   `run_summary` turns that into phi2 from the AMBIT and from the lifetimes, with warnings.
4. `_show_run` writes `results.json` and draws `raw.png` and the result tabs. A failed analysis
   still shows the raw data.
5. The upload starts from the button, or after every run when "upload every run" is on.
   `upload.write_run` writes `openjii_upload.ndjson`, and `upload_and_wait` sends it and polls until
   the platform answers. The outcome goes to `openjii_upload.json` in that run's folder, as
   [openjii-upload.md](openjii-upload.md) describes.

Run folders live in `runs/<YYYYmmdd_HHMMSS>_gui/pauses_<label>_<HHMMSS>/`. From source that is
`scripts/ambit-gate/runs/`, which git ignores. A standalone build on Windows or Linux writes `runs/`
next to the executable. On macOS it writes `~/Documents/phi2gui/runs/`, because the `.app` bundle
must stay unchanged. "Open saved run..." and `--open` analyse a run folder again without
instruments.

### Threads

Anything that waits on hardware or the network runs in a `Task`, a `QThread` that calls `fn(log)`
and emits `done(result)` or `failed(text)`. Connect those signals to bound methods of `MainWindow`,
so Qt queues each call onto the GUI thread. A lambda or a `functools.partial` has no receiver
object and can run on the worker thread, and touching a widget there crashes the app. When a slot
needs more than the result, return a tuple from `fn`, the way `_upload_current` returns the run
folder with the upload record.

`app.py` selects matplotlib's Agg backend before it imports the analysis modules, which draw with
pyplot. The window embeds their figures itself.

### Settings and user files

| What | Where |
| --- | --- |
| phi2 GUI settings and API keys per environment | `phi2gui/phi2gui.json` in the user config folder from `openjii.config_dir()`: `%APPDATA%` on Windows, `~/.config` on Linux, `~/Library/Application Support` on macOS |
| IRFs added in the phi2 GUI | `irf-library/` in the same config folder |
| The standalone build's log | `phi2gui.log` next to the executable, or `~/Library/Logs/phi2gui.log` on macOS |
| flimgui settings | `flimgui_settings.json` in the repository root, which git ignores |
| flimgui IRFs | `irf-library/` in the repository |

## flimgui

A live acquisition and fitting app for the Felix alone, separate from the phi2 GUI. Start it with
`uv run python -m flimgui`, or add `--sim` for a simulated device. Its `fitting.py` is also the
optional tail fit that `ambit_gate` imports, and `tests/test_fitting.py` covers it.
