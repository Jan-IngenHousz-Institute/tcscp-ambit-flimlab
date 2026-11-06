"""The pinned-hash check on the Felix driver installer (packaging/fetch_driver.py), without network.

Run with:  uv run pytest
"""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fetch_driver", ROOT / "packaging" / "fetch_driver.py")
fd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fd)


def test_matching_file_is_kept(tmp_path):
    f = tmp_path / "driver.exe"
    f.write_bytes(b"installer")
    assert fd.verified(f, hashlib.sha256(b"installer").hexdigest())
    assert f.exists()


def test_other_file_is_deleted(tmp_path):
    f = tmp_path / "driver.exe"
    f.write_bytes(b"something else")
    assert not fd.verified(f, hashlib.sha256(b"installer").hexdigest())
    assert not f.exists()


def test_missing_file_is_not_verified(tmp_path):
    assert not fd.verified(tmp_path / "driver.exe")


def test_pinned_hash_is_a_sha256():
    assert len(fd.SHA256) == 64 and int(fd.SHA256, 16) >= 0
