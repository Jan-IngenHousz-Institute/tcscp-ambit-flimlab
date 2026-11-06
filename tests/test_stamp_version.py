"""Writing the release version into the build (packaging/stamp_version.py).

Run with:  uv run pytest
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("stamp_version", ROOT / "packaging" / "stamp_version.py")
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)


@pytest.fixture
def init(tmp_path, monkeypatch):
    path = tmp_path / "__init__.py"
    path.write_bytes(b'\xef\xbb\xbf"""doc"""\n\n__version__ = "0.0.0"\n# tail\n')
    monkeypatch.setattr(sv, "INIT", path)
    return path


@pytest.mark.parametrize("version", ["0.2.0", "1.0.0", "0.1.7+pr.6"])
def test_stamp_keeps_bom_and_rest_of_file(init, capsys, version):
    assert sv.main(["stamp_version.py", version]) == 0
    assert capsys.readouterr().out.strip() == version
    assert init.read_bytes() == b'\xef\xbb\xbf"""doc"""\n\n__version__ = "' + version.encode() + b'"\n# tail\n'


@pytest.mark.parametrize("bad", [[], ["v0.2.0"], ["0.2"], ["0.2.0", "extra"]])
def test_refuses_anything_but_a_version(init, bad):
    with pytest.raises(SystemExit):
        sv.main(["stamp_version.py", *bad])
    assert b'__version__ = "0.0.0"' in init.read_bytes()
