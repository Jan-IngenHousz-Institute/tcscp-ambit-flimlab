"""The upload pre-check of the phi2 GUI: who may append to or create the upload table.

Run with:  python -m pytest tests  (or simply  python tests/test_upload_access.py)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "ambit-gate"))

from phi2gui.openjii import OpenJIIError  # noqa: E402
from phi2gui.upload import check_access  # noqa: E402


class FakeClient:
    """Answers the two calls check_access makes, and refuses to upload."""

    def __init__(self, caps: dict, table_id: str | None) -> None:
        self.caps, self.table_id = caps, table_id

    def capabilities(self, experiment_id: str) -> dict:
        return self.caps

    def find_table(self, experiment_id: str, table_name: str):
        return self.table_id

    def upload_files(self, *a, **k):
        raise AssertionError("upload_files must not be called by check_access")


def test_member_appends_to_existing_table():
    client = FakeClient({"canContribute": True, "canManage": False}, "11111111-2222-3333-4444-555555555555")
    assert check_access(client, "exp", "ambit_tcspc_pauses") == "11111111-2222-3333-4444-555555555555"


def test_admin_creates_missing_table():
    client = FakeClient({"canContribute": True, "canManage": True}, None)
    assert check_access(client, "exp", "ambit_tcspc_pauses") is None


def test_member_cannot_create_table():
    client = FakeClient({"canContribute": True, "canManage": False}, None)
    with pytest.raises(OpenJIIError) as info:
        check_access(client, "exp", "ambit_tcspc_pauses", experiment_name="Leaf study")
    text = str(info.value)
    assert "ambit_tcspc_pauses" in text and "'Leaf study'" in text and "manage access" in text


def test_non_member_cannot_contribute():
    client = FakeClient({"canContribute": False, "canManage": False}, None)
    with pytest.raises(OpenJIIError, match="cannot contribute"):
        check_access(client, "exp", "ambit_tcspc_pauses")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
