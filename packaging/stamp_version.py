"""Pick the release version and write it into phi2gui/__init__.py, for the release workflow.

usage: python packaging/stamp_version.py <latest tag, e.g. v0.1.4, or empty>

The next version is the latest tag's patch + 1, unless ``__version__`` in the source is already
higher (a person raised the minor or major version on purpose), in which case that one is used.
Prints the version without the ``v``.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

INIT = Path(__file__).resolve().parents[1] / "scripts" / "ambit-gate" / "phi2gui" / "__init__.py"
VERSION_LINE = re.compile(rb'^__version__ = "(\d+)\.(\d+)\.(\d+)"', re.M)


def next_version(latest_tag: str, source: tuple[int, int, int]) -> tuple[int, int, int]:
    m = re.fullmatch(r"v(\d+)\.(\d+)\.(\d+)", latest_tag.strip())
    if not m:
        return source
    major, minor, patch = map(int, m.groups())
    bumped = (major, minor, patch + 1)
    return max(bumped, source)


def main(argv: list[str]) -> int:
    raw = INIT.read_bytes()                       # bytes: the file starts with a UTF-8 BOM
    m = VERSION_LINE.search(raw)
    if not m:
        raise SystemExit(f"no __version__ line in {INIT}")
    source = tuple(int(g) for g in m.groups())
    version = ".".join(map(str, next_version(argv[1] if len(argv) > 1 else "", source)))
    INIT.write_bytes(VERSION_LINE.sub(f'__version__ = "{version}"'.encode(), raw, count=1))
    print(version)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
