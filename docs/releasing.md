# Releasing the phi2 GUI

`.github/workflows/gui-release.yml` builds the standalone phi2 GUI for Windows, Linux and macOS and
publishes it on the Releases page. Nobody builds or tags a release by hand.

## What happens on a merge

| Event | Result |
| --- | --- |
| Pull request that touches the GUI | Tests, three builds, the bundle check and the smoke test. The zips stay on the run as artifacts for 30 days. |
| Merge to `main` that touches the GUI | The same, then a GitHub release `v<version>` with three zips and `SHA256SUMS`, marked Latest |
| Change to CAD, docs or notebooks only | No build and no release |
| Manual run, `workflow_dispatch` | Builds and artifacts, never a release |

The `paths` list at the top of the workflow defines "touches the GUI". It covers the library and
`phi2gui`, `flimgui`, `irf-library`, `packaging`, `pyproject.toml`, `uv.lock`,
`build_phi2gui.bat` and the workflow file itself.

| Zip | Built on | Runs on |
| --- | --- | --- |
| `phi2gui-v<version>-win64.zip` | `windows-latest` | Windows 10 and 11, 64-bit |
| `phi2gui-v<version>-linux-x64.zip` | `ubuntu-22.04` | Linux x64 with glibc 2.35 or newer |
| `phi2gui-v<version>-macos-arm64.zip` | `macos-15` | macOS on Apple silicon |

The Linux build runs on the oldest supported Ubuntu, because a binary built against an older glibc
runs on newer ones and not the other way round.

## The version

`packaging/stamp_version.py` takes the latest `v*` tag and adds one to the patch number, unless
`__version__` in `scripts/ambit-gate/phi2gui/__init__.py` is already higher. One job picks the
number, and every build writes it into its executable, so the three apps and the release agree.

- For a normal change, leave `__version__` alone. It lags behind the tags on `main`, and that is
  expected.
- For a new minor or major version, set `__version__` to it, say `0.2.0`, in the pull request.

If another commit already holds the tag, the run fails instead of overwriting it, and running the
workflow again builds the next number. A rerun for a commit that is already released does nothing.

## Self-contained builds

A user downloads a zip, extracts it and starts the app without installing anything. Each bundle
carries Python, Qt, numpy, scipy, matplotlib, `pylibtdc` and their native libraries. It leaves out
only what has to match the user's machine:

- On Windows, the system DLLs, and the ICU that Windows 10 1903 and later provides for Qt.
- On Linux, glibc, the GPU driver stack of `libGL`, `libEGL` and `libdrm`, and the display server's
  client libraries, `libxcb` and `libwayland-*`. Every desktop Linux has these.
- On macOS, the system frameworks and `/usr/lib`.

PyInstaller bundles a library only when the build machine has it, so the Linux job first installs
the X11 libraries Qt's window plugin uses, `libxcb-cursor0` among them. `packaging/check_bundle.py`
then fails a Linux or macOS build that needs any library outside those lists and lacks it. On
Windows the smoke test on a fresh runner covers the same ground.

## The smoke test

The workflow starts each frozen app with `--smoke-test FILE`. The app loads Qt, `pylibtdc` with its
libraries on Windows and Linux, and the main window offscreen, writes `ok <version>` to `FILE` and
exits. Anything else fails the build and prints `phi2gui.log`. This catches a bundle that cannot
start, such as v0.1.3, which stopped with "DLL load failed while importing QtWidgets". It cannot
catch a fault that needs the instruments, so try a new release on the set-up before relying on it.

## Windows specifics

PyInstaller resolves the DLLs a binary needs through `PATH`, and a DLL from elsewhere on the PC can
end up in the bundle. v0.1.3, built from an Anaconda prompt, shipped conda's ICU and did not start
on any PC. Three defences guard against that, each documented where it lives:

- The Windows build runs with `PATH` cut down to Windows, the build Python and `.venv`. Both the
  workflow and `build_phi2gui.bat` do this.
- On Windows, `packaging/phi2gui.spec` never bundles ICU. It keeps one Visual C++ runtime,
  PySide6's. It also drops or replaces any binary from outside this Python, and it counts both
  `sys.prefix` and `sys.base_prefix` as this Python.
- The smoke test.

`build_phi2gui.bat` makes the same Windows zip on your own PC with uv, smoke test included, so you
can try a build before merging. On Linux or macOS, run
`uv run --group build python -m PyInstaller packaging/phi2gui.spec` and then
`uv run python packaging/check_bundle.py dist/phi2gui`, or `dist/phi2gui.app` on macOS.

Nothing is code-signed. Windows SmartScreen may warn on first start, and macOS needs a right-click
and Open, because the app is not notarized. `packaging/README-dist.txt` ships as `README.txt` in
every zip and tells users what to do.

## A new pylibtdc

TEDIEL publishes pylibtdc as a zip of wheels at https://tediel.com/resources, and only for
Windows and Linux. uv cannot read a zip of wheels, so this repository republishes the two
Python 3.13 wheels as assets of a `deps-pylibtdc-<date>` release. That release is a pre-release and
never Latest. To move to a new pylibtdc:

1. Download TEDIEL's new zip and note its SHA-256.
2. Create the release with the two wheels, unchanged:
   `gh release create deps-pylibtdc-<date> <win_amd64 wheel> <manylinux x86_64 wheel> --prerelease --latest=false`.
   Put the zip's URL and hash and each wheel's hash in the notes, as `deps-pylibtdc-2026-09-01` does.
3. Point the two `[tool.uv.sources]` URLs in `pyproject.toml` at the new assets.
4. Run `uv lock` and `uv sync`, then the tests and
   `uv run python -m phi2gui --smoke-test <file>` from `scripts/ambit-gate`.
5. Update the version in `packaging/README-dist.txt`, where it says "Python API 2026.09.01".
