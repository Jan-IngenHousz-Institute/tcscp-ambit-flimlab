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
MACHO_MAGIC = (b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe", b"\xce\xfa\xed\xfe")
# Load commands whose library must be there. LC_LOAD_WEAK_DYLIB is optional by definition.
MAC_REQUIRED = {"LC_LOAD_DYLIB", "LC_REEXPORT_DYLIB", "LC_LOAD_UPWARD_DYLIB"}


def _run(*cmd: str) -> str:
    """A tool's output. A tool that fails fails the check, so it never passes by seeing nothing."""
    done = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise SystemExit(f"{' '.join(cmd)} failed ({done.returncode}): {done.stderr.strip()[:300]}")
    return done.stdout


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


def load_commands(otool_l: str) -> tuple[list[str], list[str]]:
    """The required libraries and the LC_RPATH entries in `otool -l` output."""
    libs: list[str] = []
    rpaths: list[str] = []
    cmd = ""
    for line in otool_l.splitlines():
        words = line.split()
        if len(words) >= 2 and words[0] == "cmd":
            cmd = words[1]
        elif len(words) >= 2 and words[0] == "name" and cmd in MAC_REQUIRED:
            libs.append(line.split("name", 1)[1].rsplit("(offset", 1)[0].strip())
        elif len(words) >= 2 and words[0] == "path" and cmd == "LC_RPATH":
            rpaths.append(line.split("path", 1)[1].rsplit("(offset", 1)[0].strip())
    return libs, rpaths


def _expand(path: str, loader: Path, executable_dir: Path) -> Path:
    if path.startswith("@loader_path"):
        return loader.parent / path[len("@loader_path") :].lstrip("/")
    if path.startswith("@executable_path"):
        return executable_dir / path[len("@executable_path") :].lstrip("/")
    return Path(path)


def macos(root: Path) -> list[str]:
    """Every required library resolves to a file inside the bundle, or to the operating system.
    @rpath, @loader_path and @executable_path are resolved the way dyld does, through the binary's
    own LC_RPATH entries, so a reference to a library nobody bundled fails here."""
    inside = root.resolve()
    executable_dir = root / "Contents" / "MacOS" if root.suffix == ".app" else root
    missing: dict[str, set[str]] = {}
    for f in (p for p in root.rglob("*") if p.is_file() and not p.is_symlink()):
        if f.read_bytes()[:4] not in MACHO_MAGIC:
            continue
        libs, rpaths = load_commands(_run("otool", "-l", str(f)))
        for lib in libs:
            if lib.startswith(MAC_SYSTEM):
                continue
            if lib.startswith("@rpath/"):
                found = [_expand(r, f, executable_dir) / lib[len("@rpath/") :] for r in rpaths]
            else:
                found = [_expand(lib, f, executable_dir)]
            if not any(c.exists() and (c.resolve() == inside or inside in c.resolve().parents) for c in found):
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
