# Releasing the phi2 GUI

The standalone Windows build (`phi2gui-v<version>-win64.zip` on the Releases page) is made by
`.github/workflows/gui-release.yml`. Nobody builds or tags a release by hand.

## What happens on a merge

| Event | Result |
| --- | --- |
| Pull request touching the GUI | Tests, Windows build, smoke test; the zip stays on the run as an artifact for 30 days |
| Merge to `main` touching the GUI | The same, then a GitHub release `v<version>` with the zip and `SHA256SUMS`, marked Latest |
| Change to CAD, docs or notebooks only | No build, no release |
| Manual run (`workflow_dispatch`) | Build and artifact, never a release |

"Touching the GUI" is the `paths` list at the top of the workflow: the library and `phi2gui`,
`flimgui`, `irf-library`, `packaging`, `pyproject.toml`, `uv.lock`, `build_phi2gui.bat` and the
workflow file itself.

## The version

`packaging/stamp_version.py` takes the latest `v*` tag and adds one to the patch number, unless
`__version__` in `scripts/ambit-gate/phi2gui/__init__.py` is already higher. It writes the chosen
number into the build, so the exe and the release always agree.

- A normal change: leave `__version__` alone. It lags behind the tags on `main`, which is expected.
- A new minor or major version: set `__version__` to it (say `0.2.0`) in the pull request.

A run that finds its tag already taken by another commit fails rather than overwrite it; run the
workflow again to build the next number. A rerun of a commit that is already released does nothing.

## The smoke test

The workflow starts the frozen exe with `--smoke-test FILE`. It loads Qt, `pylibtdc` with its DLLs
and the main window offscreen, writes `ok <version>` to `FILE` and exits; anything else fails the
build, and the job prints `phi2gui.log`. This catches a bundle that cannot start, such as v0.1.3's
"DLL load failed while importing QtWidgets". It cannot catch a fault that needs the instruments: try
a new release on the set-up before relying on it.

## Building on Windows

PyInstaller resolves the DLLs a binary needs through `PATH`, and a DLL from elsewhere on the PC can
end up in the bundle. v0.1.3, built from an Anaconda prompt, shipped conda's ICU and did not start
on any PC. Three defences, each documented where it lives:

- The build runs with `PATH` cut down to Windows, the build Python and `.venv` (the workflow's build
  step and `build_phi2gui.bat`).
- `packaging/phi2gui.spec` never bundles ICU, keeps one Visual C++ runtime (PySide6's), and drops
  or replaces any binary from outside this Python, counting both `sys.prefix` and
  `sys.base_prefix` as this Python.
- The smoke test above.

`build_phi2gui.bat` makes the same zip on your own PC with uv, smoke test included, to try a build
before merging.

The exe is not code-signed, so Windows SmartScreen may warn on first start ("More info" > "Run
anyway"); `packaging/README-dist.txt`, shipped as `README.txt` in the zip, says so.

## A new pylibtdc

TEDIEL publishes pylibtdc as a zip of wheels at https://tediel.com/resources. To move to a new one:

1. In `packaging/fetch_pylibtdc.py`, set `ZIP_NAME`, `URL`, `SHA256` and the two `WHEELS` names
   (the Python 3.13 Windows and Linux wheels inside the zip).
2. Point the two `[tool.uv.sources]` paths in `pyproject.toml` at the new wheel names.
3. Delete `wheelhouse/`, run `uv run --no-project packaging/fetch_pylibtdc.py`, then `uv lock` and
   `uv sync`, and run the tests and `python -m phi2gui --smoke-test <file>` from `scripts/ambit-gate`.
4. Update the version in `packaging/README-dist.txt` ("Python API 2026.09.01").
