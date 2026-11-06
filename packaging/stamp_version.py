"""Write a version into phi2gui/__init__.py, for the release workflow and local builds.

    python packaging/stamp_version.py 0.2.0

semantic-release picks the number from the Conventional Commit titles merged since the last tag
(see .releaserc.json and docs/releasing.md), and every build writes it in, so the executables and
the release agree. A pull request's build uses the last tag plus "+pr.<number>".
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

INIT = Path(__file__).resolve().parents[1] / "scripts" / "ambit-gate" / "phi2gui" / "__init__.py"
VERSION_LINE = re.compile(rb'^__version__ = "[^"]*"', re.M)
VERSION = re.compile(r"\d+\.\d+\.\d+([-+][0-9A-Za-z.+-]+)?")


def main(argv: list[str]) -> int:
    if len(argv) != 2 or not VERSION.fullmatch(argv[1]):
        raise SystemExit(f"usage: stamp_version.py <version like 1.2.3 or 1.2.3+pr.4>, got {argv[1:]}")
    raw = INIT.read_bytes()  # bytes, so a UTF-8 BOM survives if the file has one
    if not VERSION_LINE.search(raw):
        raise SystemExit(f"no __version__ line in {INIT}")
    INIT.write_bytes(VERSION_LINE.sub(f'__version__ = "{argv[1]}"'.encode(), raw, count=1))
    print(argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
