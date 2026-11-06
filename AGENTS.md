# Agent context

A lab repository. An AMBIT fluorimeter and a TEDIEL Felix TDC measure the PSII quantum yield of a
leaf, phi2, from fluorescence intensity and lifetime. The code drives real instruments, and an agent
almost never has them. Know which checks you can run, listed below, and say plainly in your handoff
what only the set-up can test.

## Where things live

| Path | What it is | Read first |
| --- | --- | --- |
| `scripts/ambit-gate/` | Measurement and analysis library, one module per experiment, and the numbered notebooks | `docs/architecture.md` |
| `scripts/ambit-gate/phi2gui/` | The phi2 GUI in PySide6, the actinic-pauses measurement with upload to openJII | `docs/architecture.md`, `docs/openjii-upload.md` |
| `flimgui/` | A separate live FLIM GUI for the Felix alone. `ambit_gate` also imports its `fitting.py`. | `docs/architecture.md` |
| `packaging/`, `build_phi2gui.bat`, `.github/workflows/gui-release.yml` | The standalone builds for Windows, Linux and macOS, and the automatic releases | `docs/releasing.md` |
| `CAD/` | The light-tight leaf jig in FreeCAD | `CAD/README.md` |
| `docs/` | Developer pages and the vendor manuals | `docs/README.md` |

`README.md` is for the people who use the set-up. Developer detail goes in `docs/`.

## Setup and checks

uv manages Python 3.13 and every package, pinned in `uv.lock`. From the repository root:

```bash
uv sync                       # every package into .venv, TEDIEL's pylibtdc included
uv run nbstripout --install   # once per clone, so notebook outputs stay out of commits
uv run pytest -q
```

pylibtdc is not on PyPI. `pyproject.toml` points uv at TEDIEL's wheels on this repository's
`deps-pylibtdc-*` release. TEDIEL makes them for Windows and Linux x64 only, so on a Mac the `felix`
group installs nothing.

These checks work without instruments, on Linux and macOS too:

- `uv run pytest -q` covers fitting, the upload row and access rules, and the release version choice.
- From `scripts/ambit-gate`, `QT_QPA_PLATFORM=offscreen uv run python -m phi2gui --smoke-test
  /tmp/phi2gui-ok.txt` writes `ok <version>` once Qt, pylibtdc and the main window load.
- `QT_QPA_PLATFORM=offscreen uv run python -m flimgui --sim` runs flimgui on a simulated Felix. It
  runs until killed, so a `timeout` exit means it started.
- `chop_sim.simulate_run(out_dir, ambit_chop.phi2_protocol(100))` writes a synthetic chopped run
  with its known truth for `ambit_chop.analyse`. It is a library without a CLI, and its docstring
  has the example.
- `uv run python -m phi2gui --open <run folder>` re-analyses a saved run. Run folders are not in git,
  so ask the person for one.

A measurement, the driver install and anything that talks to the AMBIT or the Felix need the set-up.
Hand those to the person as steps to try.

## Conventions

- **Sibling imports.** The modules in `scripts/ambit-gate` import each other by name, with that
  folder on `sys.path`. uv does not install the repository as a package. Run scripts from that folder.
- **Docstrings carry the physics.** Module and function docstrings explain the measurement, the
  timing and the reason for each choice in full sentences. Inline comments are short and state a
  constraint. Match that, and keep units in names, such as `t_s`, `bin_ps`, `window_s`, `laser_hz`.
- **Notebooks** go through nbstripout, set in `.gitattributes`, so their outputs never reach a commit.
- **Dependencies** change with `uv add` and `uv remove`. Commit `pyproject.toml` and `uv.lock` together.
- **Commits** use conventional subjects, such as `feat(upload): ...`, `fix(gui): ...`, `ci: ...`,
  `build: ...` and `docs: ...`.
- **Pull requests** go to `main` with green CI. A merge that touches the GUI publishes a release for
  all three platforms by itself, so `main` must always be releasable. Leave `__version__` alone
  unless the change starts a new minor or major version, as `docs/releasing.md` explains.
- **Self-contained builds.** A standalone zip carries everything except the operating system's own
  libraries. `packaging/check_bundle.py` fails a Linux or macOS build that misses one, and the
  smoke test starts every build on a fresh runner.
- **Plain prose.** Docs follow openJII's unslop rules: no em dashes and no parentheses in prose,
  active voice, sentence-case headings.
- `.bat` files are Windows batch files with CRLF line endings, set in `.gitattributes`. Keep them
  compatible with cmd.

## Gotchas

- Worker results reach the window through Qt signals. Connect `Task.done` and `Task.failed` to
  bound `MainWindow` methods. A lambda runs on the worker thread and crashes when it touches a
  widget. Return extra data from the task as a tuple.
- `phi2gui/app.py` selects matplotlib's Agg backend before it imports the analysis modules. Keep
  that order.
- Renaming or retyping a key of the uploaded row breaks existing openJII tables. Raise
  `upload.ROW_FORMAT` and change `DEFAULT_TABLE`, as `docs/openjii-upload.md` describes.
- Test uploads on openJII dev with a dev API key and a `QA-<date>-` experiment. Production keys
  belong to the people measuring.
- A Felix opens in one program at a time, so close TDC-Software first. The AMBIT takes at most 16
  `arrun` lines, `ambit_gate.MAX_SEGMENTS`.
- Conda's Python makes PySide6 fail with "DLL load failed while importing QtWidgets". Use uv's
  `.venv`, and build releases only through the workflow or `build_phi2gui.bat`.
- `flimgui` writes `flimgui_settings.json` into the repository root, which git ignores. The phi2 GUI
  keeps its settings and API keys in `%APPDATA%\phi2gui\phi2gui.json`, or the platform's equivalent.
