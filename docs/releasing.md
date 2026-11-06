# Releasing the phi2 GUI

`.github/workflows/gui-release.yml` builds the standalone phi2 GUI for Windows, Linux and macOS and
publishes it on the Releases page. Nobody builds or tags a release by hand.

## Pull requests and checks

Every change reaches `main` through a pull request, and GitHub merges it as one squash commit whose
subject is the PR title. A ruleset on `main` enforces that and lets nobody bypass it. Before the
merge button works, these checks must pass:

| Check | Workflow | What it runs |
| --- | --- | --- |
| `PR title` | `pr-title.yml` | The title is a Conventional Commit, such as `fix(upload): keep each record with its run` |
| `lint` | `ci.yml` | `ruff check` and `ruff format --check` |
| `test` | `ci.yml` | `pytest` |
| `Build (win64)`, `Build (linux-x64)`, `Build (macos-arm64)` | `gui-release.yml` | The tests on each platform, the build, the bundle check and the smoke test |

The three builds keep their zips on the run as artifacts for 30 days, so a reviewer can try the app
before merging. `uv run ruff check --fix . && uv run ruff format .` fixes most lint findings
locally.

## What happens on a merge

semantic-release reads the subjects merged since the last tag and decides, following openJII's rules
in `.releaserc.json`:

| PR title type | Release |
| --- | --- |
| `feat` | A new minor version, 0.1.7 to 0.2.0 |
| `fix`, `perf`, `revert` | A new patch version, 0.1.7 to 0.1.8 |
| Any type with `!`, such as `feat!:`, or a `BREAKING CHANGE:` footer | A new major version |
| `docs`, `style`, `chore`, `refactor`, `test`, `build`, `ci` | None |

When a release is due, the workflow builds the three zips under that version, and semantic-release
tags the commit `v<version>` and publishes the release with the zips, `SHA256SUMS` and notes grouped
by type. When none is due, nothing is built or published. A manual run, `workflow_dispatch`, builds
the zips and never releases.

| Zip | Built on | Runs on |
| --- | --- | --- |
| `phi2gui-v<version>-win64.zip` | `windows-latest` | Windows 10 and 11, 64-bit |
| `phi2gui-v<version>-linux-x64.zip` | `ubuntu-22.04` | Linux x64 with glibc 2.35 or newer |
| `phi2gui-v<version>-macos-arm64.zip` | `macos-15` | macOS on Apple silicon |

The Linux build runs on the oldest supported Ubuntu, because a binary built against an older glibc
runs on newer ones and not the other way round.

## The version

The version comes from the PR titles, so nobody edits it. `__version__` in
`scripts/ambit-gate/phi2gui/__init__.py` stays `0.0.0` in git. One job picks the number, and every
build writes it into its executable with `packaging/stamp_version.py`, so the three apps and the
release agree. A pull request builds as the last tag plus `+pr.<number>`, such as `0.1.7+pr.8`.

## When a release fails halfway

semantic-release pushes the tag `v<version>` before it creates the GitHub release and uploads the
zips. If the run fails in between, the tag exists and the release is missing, still a draft, or
short of assets. Re-run all jobs of that workflow run. The release-version job then finds the tag
at the merged commit with an unfinished release, builds that same version again, and the release
job uploads the zips and publishes it. The last step passes only when the release is published with
exactly the three zips and `SHA256SUMS`.

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

Nothing is code-signed. Windows SmartScreen may warn on first start. macOS blocks the first start
because the app is not notarized, and since macOS 15 the only way past that is System Settings >
Privacy & Security > Open Anyway. `packaging/README-dist.txt` ships as `README.txt` in
every zip and tells users what to do.

## Supply chain

Everything the release job runs or ships is pinned, so a change upstream cannot reach a release
without a pull request here:

- `uv.lock` pins every Python package by hash, the pylibtdc wheels included.
- `packaging/fetch_driver.py` takes TEDIEL's Felix driver installer only at its pinned SHA-256,
  whether it downloads it or finds it in `docs/`. The phi2 GUI runs that installer with
  administrator rights, so the spec refuses to bundle any other file. When TEDIEL publishes a new
  installer, check it and change `SHA256` in a pull request.
- semantic-release and its preset come from `packaging/release-tools/package-lock.json` through
  `npm ci`. Update them with `npm install --save-exact` in that folder and commit the lockfile.
- Every action in `.github/workflows/` is pinned to a commit SHA, with the version in a comment.

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
