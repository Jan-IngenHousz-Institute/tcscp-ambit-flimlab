"""The release workflow's version choice (packaging/stamp_version.py).

Run with:  python -m pytest tests
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("stamp_version", ROOT / "packaging" / "stamp_version.py")
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)


def test_next_patch_after_latest_tag():
    assert sv.next_version("v0.1.4", (0, 1, 4)) == (0, 1, 5)


def test_source_version_wins_when_raised_by_hand():
    assert sv.next_version("v0.1.4", (0, 2, 0)) == (0, 2, 0)


def test_no_tag_takes_source_version():
    assert sv.next_version("", (0, 1, 4)) == (0, 1, 4)


def test_stamp_keeps_bom_and_rest_of_file(tmp_path, monkeypatch, capsys):
    init = tmp_path / "__init__.py"
    init.write_bytes(b'\xef\xbb\xbf"""doc"""\n\n__version__ = "0.1.4"\n# tail\n')
    monkeypatch.setattr(sv, "INIT", init)
    assert sv.main(["stamp_version.py", "v0.1.4"]) == 0
    assert capsys.readouterr().out.strip() == "0.1.5"
    assert init.read_bytes() == b'\xef\xbb\xbf"""doc"""\n\n__version__ = "0.1.5"\n# tail\n'


if __name__ == "__main__":
    import pytest

    sys.exit(pytest.main([__file__, "-q"]))
