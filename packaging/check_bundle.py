"""Fail when a built bundle needs a library it does not carry, beyond the operating system's own.

    python packaging/check_bundle.py dist/phi2gui            # Linux
    python packaging/check_bundle.py dist/phi2gui.app        # macOS

The standalone build must run on a PC with nothing installed. A few libraries cannot travel with
it, because they have to match the machine: on Linux the C library, the GPU driver stack (libGL,
libEGL, libdrm) and the display server's client libraries (libxcb, libwayland-*); on macOS the
system frameworks and /usr/lib. Everything else the binaries link against must be in the bundle.
PyInstaller only bundles what is installed on the build machine, so a library missing on the CI
runner would otherwise leave the bundle depending on the user's system without anyone noticing.
Windows has no equivalent check here: the smoke test runs the exe on a fresh runner instead.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

LINUX_SYSTEM = re.compile(
    r"^(ld-linux.*|linux-vdso.*|lib(c|m|dl|pthread|rt|util|resolv|anl|nsl)\.so\..*"
    r"|lib(GL|EGL|GLX|GLdispatch|OpenGL|gbm|drm)\.so\..*"
    r"|libxcb\.so\..*|libwayland-(client|cursor|egl|server)\.so\..*)$"
)
MAC_SYSTEM = ("/System/Library/", "/usr/lib/")
MAC_LOCAL = ("@rpath/", "@loader_path/", "@executable_path/")


def _run(*cmd: str) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=False).stdout


def linux(root: Path) -> list[str]:
    bundled = {p.name for p in root.rglob("*") if p.is_file()}
    missing: dict[str, set[str]] = {}
    for f in (p for p in root.rglob("*") if p.is_file()):
        if f.read_bytes()[:4] != b"\x7fELF":
            continue
        for so in re.findall(r"NEEDED\s+(\S+)", _run("objdump", "-p", str(f))):
            if so not in bundled and not LINUX_SYSTEM.match(so):
                missing.setdefault(so, set()).add(str(f.relative_to(root)))
    return [f"{so}  (needed by {', '.join(sorted(users)[:3])})" for so, users in sorted(missing.items())]


def macos(root: Path) -> list[str]:
    missing: dict[str, set[str]] = {}
    for f in (p for p in root.rglob("*") if p.is_file() and not p.is_symlink()):
        head = f.read_bytes()[:4]
        if head not in (b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe", b"\xce\xfa\xed\xfe"):
            continue
        for line in _run("otool", "-L", str(f)).splitlines()[1:]:
            lib = line.strip().split(" (")[0]
            if lib and not lib.startswith(MAC_SYSTEM + MAC_LOCAL):
                missing.setdefault(lib, set()).add(str(f.relative_to(root)))
    return [f"{lib}  (needed by {', '.join(sorted(users)[:3])})" for lib, users in sorted(missing.items())]


def main(argv: list[str]) -> int:
    root = Path(argv[1])
    problems = macos(root) if sys.platform == "darwin" else linux(root)
    if problems:
        print(f"{root} needs libraries it does not carry:")
        print("\n".join(f"  {p}" for p in problems))
        return 1
    print(f"{root}: every library outside the operating system is bundled")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
