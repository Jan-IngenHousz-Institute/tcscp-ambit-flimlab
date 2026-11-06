# Agent context

A lab repository: an AMBIT fluorimeter and a TEDIEL Felix TDC measure PSII quantum yield (phi2)
from fluorescence intensity and lifetime. The code drives real instruments on a Windows PC, and an
agent almost never has them, so know which checks you can run (below) and say plainly in your
handoff what only the set-up can test.

## Where things live

| Path | What it is | Read first |
| --- | --- | --- |
| `scripts/ambit-gate/` | Measurement and analysis library, one module per experiment, plus numbered notebooks | `docs/architecture.md` |
| `scripts/ambit-gate/phi2gui/` | The phi2 GUI (PySide6): the actinic-pauses measurement with upload to openJII | `docs/architecture.md`, `docs/openjii-upload.md` |
| `flimgui/` | Separate live FLIM GUI for the Felix alone; its `fitting.py` is shared | `docs/architecture.md` |
| `packaging/`, `build_phi2gui.bat`, `.github/workflows/gui-release.yml` | Windows standalone build and automatic releases | `docs/releasing.md` |
| `CAD/` | The light-tight leaf jig (FreeCAD) | `CAD/README.md` |
| `docs/` | Developer pages and vendor manuals | `docs/README.md` |

`README.md` is for the people using the set-up; keep it about using, and put developer detail in
`docs/`.

## Setup and checks

uv manages Python 3.13 and every package (`pyproject.toml`, `uv.lock`). From the repository root:

```bash
uv run --no-project packaging/fetch_pylibtdc.py   # once per clone, see below
uv sync
uv run nbstripout --install                        # once per clone: notebook outputs stay local
uv run pytest -q
```

pylibtdc, the Felix bindings, is not on PyPI. The fetch script downloads TEDIEL's zip from
tediel.com into `docs/` (or uses the copy already there), checks its SHA-256, and extracts the
Python 3.13 wheels into `wheelhouse/`; both are gitignored.

The checks that work without instruments, on Linux too:

- `uv run pytest -q`: fitting, the upload row and access rules, the release version choice.
- The phi2 GUI's startup, offscreen: from `scripts/ambit-gate`,
  `QT_QPA_PLATFORM=offscreen uv run python -m phi2gui --smoke-test /tmp/phi2gui-ok.txt` writes
  `ok <version>` when Qt, pylibtdc and the main window load.
- `QT_QPA_PLATFORM=offscreen uv run python -m flimgui --sim` runs flimgui on a simulated Felix (it
  runs until killed; a `timeout` exit means it started).
- `chop_sim.simulate_run(out_dir, ambit_chop.phi2_protocol(100))` writes a synthetic chopped run
  with its known truth for `ambit_chop.analyse` (a library, no CLI; its docstring has the example).
- `python -m phi2gui --open <run folder>` re-analyses a saved run. Run folders are gitignored
  (`scripts/ambit-gate/runs/`), so ask the person for one.

A measurement, the driver install and anything that talks to the AMBIT or the Felix need the
set-up. Hand those to the person as steps to try.

## Conventions

- **Sibling imports.** The modules in `scripts/ambit-gate` import each other by name with that
  folder on `sys.path`; nothing is installed as a package. Run scripts from that folder.
- **Docstrings carry the physics.** Module and function docstrings explain the measurement, the
  timing and the reason for each choice in full sentences; inline comments are short and state a
  constraint. Match that, and keep units in names (`t_s`, `bin_ps`, `window_s`, `laser_hz`).
- **Notebooks** go through nbstripout (`.gitattributes`, installed in the setup above), so their
  outputs never reach a commit.
- **Dependencies** change with `uv add` / `uv remove`; commit `pyproject.toml` and `uv.lock` together.
- **Commits** use conventional subjects: `feat(upload): ...`, `fix(gui): ...`, `ci: ...`,
  `build: ...`, `docs: ...`.
- **Pull requests** go to `main` with green CI. A merge that touches the GUI publishes a release by
  itself, so `main` is always releasable. Leave `__version__` alone unless the change is a new minor
  or major version (`docs/releasing.md`).
- `.bat` files are Windows batch with CRLF (`.gitattributes`); keep them cmd-compatible.

## Gotchas

- `uv sync --locked` reads every path source, so it fails without `wheelhouse/` even when the
  `felix` group is not wanted: fetch pylibtdc first, or use `uv sync --frozen --no-group felix`.
- Worker results reach the window through Qt signals. Connect `Task.done` and `Task.failed` to bound
  `MainWindow` methods; a lambda runs on the worker thread and crashes when it touches a widget.
  Return extra data from the task as a tuple.
- `phi2gui/app.py` selects matplotlib's Agg backend before importing the analysis modules; keep
  that order.
- Renaming or retyping a key of the uploaded row breaks existing openJII tables: raise
  `upload.ROW_FORMAT` and change `DEFAULT_TABLE` (`docs/openjii-upload.md`).
- Test uploads on openJII **dev** with a dev API key and a `QA-<date>-` experiment; a production key
  is for the people measuring.
- A Felix opens in one program at a time; TDC-Software must be closed. The AMBIT takes at most 16
  `arrun` lines (`ambit_gate.MAX_SEGMENTS`).
- Conda's Python makes PySide6 fail with "DLL load failed while importing QtWidgets". Use uv's
  `.venv`, and build releases only through the workflow or `build_phi2gui.bat`.
- `flimgui` writes `flimgui_settings.json` into the repository root (gitignored); the phi2 GUI keeps
  its settings and API keys in `%APPDATA%\phi2gui\phi2gui.json`.
