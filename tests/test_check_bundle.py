"""The bundle check's macOS path resolution (packaging/check_bundle.py), with mocked otool output.

Run with:  uv run pytest
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_bundle", ROOT / "packaging" / "check_bundle.py")
cb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cb)

MACHO = b"\xcf\xfa\xed\xfe" + b"\0" * 28


def otool_l(libs: list[tuple[str, str]], rpaths: list[str]) -> str:
    out = ["binary:"]
    for i, (cmd, name) in enumerate(libs):
        out += [f"Load command {i}", f"          cmd {cmd}", "      cmdsize 56", f"         name {name} (offset 24)"]
    for r in rpaths:
        out += ["Load command 99", "          cmd LC_RPATH", "      cmdsize 32", f"         path {r} (offset 12)"]
    return "\n".join(out) + "\n"


@pytest.fixture
def app(tmp_path):
    """A bundle with the executable in Contents/MacOS and one library in Contents/Frameworks."""
    root = tmp_path / "phi2gui.app"
    (root / "Contents" / "MacOS").mkdir(parents=True)
    (root / "Contents" / "Frameworks").mkdir()
    (root / "Contents" / "MacOS" / "phi2gui").write_bytes(MACHO)
    (root / "Contents" / "Frameworks" / "QtCore").write_bytes(MACHO)
    return root


def check(app, monkeypatch, libs, rpaths=("@loader_path/../Frameworks",)):
    outputs = {"phi2gui": otool_l(libs, list(rpaths)), "QtCore": otool_l([], [])}
    monkeypatch.setattr(cb, "_run", lambda *cmd: outputs[Path(cmd[-1]).name])
    return cb.macos(app)


def test_rpath_library_in_the_bundle_passes(app, monkeypatch):
    assert check(app, monkeypatch, [("LC_LOAD_DYLIB", "@rpath/QtCore")]) == []


def test_missing_rpath_library_fails(app, monkeypatch):
    problems = check(app, monkeypatch, [("LC_LOAD_DYLIB", "@rpath/Missing.dylib")])
    assert len(problems) == 1 and problems[0].startswith("@rpath/Missing.dylib")


def test_loader_and_executable_paths_resolve(app, monkeypatch):
    libs = [("LC_LOAD_DYLIB", "@loader_path/../Frameworks/QtCore"), ("LC_LOAD_DYLIB", "@executable_path/x.dylib")]
    problems = check(app, monkeypatch, libs)
    assert [p.split()[0] for p in problems] == ["@executable_path/x.dylib"]


def test_system_and_weak_libraries_pass(app, monkeypatch):
    libs = [("LC_LOAD_DYLIB", "/usr/lib/libSystem.B.dylib"), ("LC_LOAD_WEAK_DYLIB", "@rpath/Optional.dylib")]
    assert check(app, monkeypatch, libs) == []


def test_library_outside_the_bundle_fails(app, monkeypatch, tmp_path):
    (tmp_path / "elsewhere.dylib").write_bytes(MACHO)
    problems = check(app, monkeypatch, [("LC_LOAD_DYLIB", str(tmp_path / "elsewhere.dylib"))])
    assert len(problems) == 1


def test_a_failing_tool_fails_the_check(monkeypatch):
    with pytest.raises(SystemExit, match="failed"):
        cb._run("false")
