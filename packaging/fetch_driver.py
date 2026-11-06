"""Put TEDIEL's Felix WinUSB driver installer into docs/, checked against a pinned SHA-256.

    python packaging/fetch_driver.py

The Windows build bundles this installer, and the phi2 GUI runs it with administrator rights on a
measurement PC that has no Felix driver yet. So the build takes it only when its SHA-256 matches
SHA256 below, whether it was just downloaded or already in docs/. A mismatch deletes the file and
fails. The pinned hash is the installer TEDIEL served on 2026-11-06, identical to the copy in the
v0.1.4 release zip that Ludovico built by hand. When TEDIEL publishes a new installer, check it and
change SHA256 in a pull request.
"""

from __future__ import annotations

import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

URL = "https://tediel.com/wp-content/Resources/Felix/Felix_WinUSB_driver.exe"
SHA256 = "421a1f7015868a581d345ec53ad63074b3f5811cea2acaa633cc25e8f933845f"
DRIVER = Path(__file__).resolve().parents[1] / "docs" / "Felix_WinUSB_driver.exe"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verified(path: Path, expected: str = SHA256) -> bool:
    """True when ``path`` holds the pinned installer. A file with another hash is deleted."""
    if not path.is_file():
        return False
    if sha256(path) == expected:
        return True
    path.unlink()
    return False


def main() -> int:
    if verified(DRIVER):
        print(f"{DRIVER.name}: SHA-256 matches")
        return 0
    DRIVER.parent.mkdir(parents=True, exist_ok=True)
    tmp = DRIVER.with_suffix(".part")
    print(f"downloading {URL}")
    with urllib.request.urlopen(URL, timeout=120) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh)
    tmp.replace(DRIVER)
    if not verified(DRIVER):
        print(f"{URL} served a file whose SHA-256 is not {SHA256}; deleted it", file=sys.stderr)
        return 1
    print(f"{DRIVER.name}: downloaded, SHA-256 matches")
    return 0


if __name__ == "__main__":
    sys.exit(main())
