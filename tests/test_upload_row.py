"""The one-row-per-run layout of the phi2 GUI's openJII upload.

Run with:  python -m pytest tests  (or simply  python tests/test_upload_row.py)
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "ambit-gate"))

from phi2gui import upload as up  # noqa: E402
from phi2gui.openjii import OpenJIIError  # noqa: E402

N_PTS, N_WIN, N_RAW = 30, 7, 266 * up.REBIN + 5      # 5 raw bins too many: dropped by the rebinning


def fake_run(tmp_path: Path) -> tuple[dict, dict]:
    run_dir = tmp_path / "pauses_leaf_101010"
    run_dir.mkdir()
    (run_dir / "run.json").write_text(json.dumps({"saved": "2026-11-05 10:10:10", "note": "n", "cmd": "arrun2,1",
                                                  "laser_hz": 8e7, "ambit": {"name": "A", "fw": "1"},
                                                  "felix": "F", "protocol": [{"type": 2}]}), encoding="utf-8")
    rng = np.random.default_rng(0)
    bins = rng.poisson(3.0, size=(N_WIN, N_RAW)).astype(np.int64)
    res = dict(run_dir=run_dir, W=0.005, bin_s=2.34375e-12, bins=bins, ts=np.linspace(-0.1, 0.4, N_WIN),
               lines=[{"line": 1, "role": "dark"}, {"line": 2, "role": "lit"}],
               ambit=dict(t=np.linspace(0, 3, N_PTS), s_630=rng.normal(size=N_PTS), r_630=rng.normal(size=N_PTS),
                          ratio=rng.normal(size=N_PTS), actinic=np.zeros(N_PTS), line=np.array([1, 2] * (N_PTS // 2))))
    summary = dict(phi2_ambit=0.5, phi2_tau_mean=float("nan"), t_amp_fm=[0.7, np.float64(0.8)], n_fo=np.int64(50),
                   fits={"reference": {"tau1": 0.15, "empty1": np.bool_(False)}, "pauses": []},
                   warnings=["w1"])
    return res, summary


def test_one_row_with_nested_traces(tmp_path):
    res, summary = fake_run(tmp_path)
    row = up.run_row(res, summary, include_felix=True, label="leaf")
    path, n = up.write_ndjson(tmp_path, [row])
    assert n == 1 and path.read_text(encoding="utf-8").count("\n") == 1
    back = json.loads(path.read_text(encoding="utf-8"))          # strict JSON: no NaN, no numpy types
    assert back["row_format"] == up.ROW_FORMAT and back["run_id"] == "pauses_leaf_101010"
    assert back["results"]["phi2_tau_mean"] is None and back["results"]["n_fo"] == 50
    assert "fits" not in back["results"] and back["fits"]["reference"]["empty1"] is False
    assert back["warnings"] == ["w1"] and back["label"] == "leaf" and back["timestamp_local"].startswith("2026-11-05T10:10:10")
    pts = back["ambit_points"]
    assert pts["n"] == N_PTS and len(pts["ratio"]) == N_PTS and pts["role"][:2] == ["dark", "lit"] and pts["line"][1] == 2
    fw = back["felix_windows"]
    assert fw["included"] is True and fw["n"] == N_WIN and fw["n_bins"] == 266 and fw["rebin"] == up.REBIN
    assert len(fw["t_start_s"]) == N_WIN and len(fw["counts"]) == N_WIN and len(fw["counts"][0]) == 266
    assert math.isclose(fw["bin_ps"], 46.875)
    raw = np.asarray(res["bins"])[:, :266 * up.REBIN]
    assert np.array_equal(np.asarray(fw["counts"]), raw.reshape(N_WIN, 266, up.REBIN).sum(2))


def test_without_felix_windows(tmp_path):
    res, summary = fake_run(tmp_path)
    row = up.run_row(res, summary, include_felix=False)
    fw = row["felix_windows"]
    assert fw["included"] is False and fw["counts"] is None and fw["n"] is None
    assert row["ambit_points"]["n"] == N_PTS


def test_row_too_large_is_refused(tmp_path, monkeypatch):
    res, summary = fake_run(tmp_path)
    monkeypatch.setattr(up, "MAX_ROW_BYTES", 1000)
    with pytest.raises(OpenJIIError, match="without the Felix windows"):
        up.write_ndjson(tmp_path, [up.run_row(res, summary)])


def _felix(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))["felix_windows"]


def test_write_run_full_resolution_keeps_every_raw_bin(tmp_path):
    res, summary = fake_run(tmp_path)
    fw = _felix(up.write_run(res["run_dir"], res, summary, felix="full", log=lambda _: None))
    assert fw["rebin"] == 1 and fw["n_bins"] == N_RAW and math.isclose(fw["bin_ps"], 2.34375)
    assert np.array_equal(np.asarray(fw["counts"]), res["bins"])


def test_write_run_rebinned_and_none(tmp_path):
    res, summary = fake_run(tmp_path)
    fw = _felix(up.write_run(res["run_dir"], res, summary, felix="rebinned"))
    assert fw["rebin"] == up.REBIN and fw["n_bins"] == 266
    fw = _felix(up.write_run(res["run_dir"], res, summary, felix="none"))
    assert fw["included"] is False and fw["counts"] is None


def test_write_run_falls_back_to_47ps_when_full_resolution_is_too_large(tmp_path, monkeypatch):
    res, summary = fake_run(tmp_path)
    rebinned = len(json.dumps(up.run_row(res, summary, rebin=up.REBIN), separators=(",", ":")))
    monkeypatch.setattr(up, "MAX_ROW_BYTES", rebinned + 10)
    logged = []
    fw = _felix(up.write_run(res["run_dir"], res, summary, felix="full", log=logged.append))
    assert fw["rebin"] == up.REBIN and fw["n_bins"] == 266
    assert len(logged) == 1 and "47 ps bins instead" in logged[0]


def test_write_run_rebinned_too_large_is_refused(tmp_path, monkeypatch):
    res, summary = fake_run(tmp_path)
    monkeypatch.setattr(up, "MAX_ROW_BYTES", 1000)
    with pytest.raises(up.RowTooLarge):
        up.write_run(res["run_dir"], res, summary, felix="rebinned")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
