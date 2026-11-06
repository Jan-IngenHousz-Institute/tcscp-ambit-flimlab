"""Put TEDIEL's pylibtdc wheels into wheelhouse/, where pyproject.toml's `felix` group finds them.

    uv run --no-project packaging/fetch_pylibtdc.py

pylibtdc (the Felix TDC bindings) is not on PyPI. TEDIEL offers it as a free zip at
https://tediel.com/resources; this takes that zip from docs/ when it is there, downloads it
otherwise, checks its SHA-256 and extracts the Windows and Linux wheels for Python 3.13.
Standard library only, so it runs before the project environment exists.
"""
from __future__ import annotations

import hashlib
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ZIP_NAME = "pylibtdc-2026-09-01.zip"
URL = f"https://tediel.com/wp-content/Resources/Common/2026.09.01/{ZIP_NAME}"
SHA256 = "1856e096897ac5eadf2f6a15a3a39af86ffa7dbfcf550cecafb59bbab655e73a"
WHEELS = (
    "pylibtdc-1.18.2.1-cp313-cp313-win_amd64.whl",
    "pylibtdc-1.18.2.1-cp313-cp313-manylinux2014_x86_64.manylinux_2_17_x86_64.whl",
)
WHEELHOUSE = ROOT / "wheelhouse"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if all((WHEELHOUSE / w).is_file() for w in WHEELS):
        print(f"pylibtdc wheels already in {WHEELHOUSE}")
        return 0
    archive = ROOT / "docs" / ZIP_NAME
    if not archive.is_file():
        archive.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {URL}")
        tmp = archive.with_suffix(".part")
        with urllib.request.urlopen(URL, timeout=120) as resp, open(tmp, "wb") as fh:
            shutil.copyfileobj(resp, fh)
        tmp.replace(archive)
    digest = _sha256(archive)
    if digest != SHA256:
        raise SystemExit(f"{archive}: SHA-256 {digest}, expected {SHA256}; delete it and run again")
    WHEELHOUSE.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive) as zf:
        for name in WHEELS:
            with zf.open(f"wheelhouse/{name}") as src, open(WHEELHOUSE / name, "wb") as dst:
                shutil.copyfileobj(src, dst)
            print(f"extracted {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
