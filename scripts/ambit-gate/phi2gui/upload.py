"""One pauses run as rows of an openJII upload table, the NDJSON file, and the upload itself.

Rows, told apart by ``kind``:
- ``run``: one row with the run's metadata (protocol, instruments, note), the phi2 results
  (``ambit_pause.run_summary``) and the fit numbers of the reference decay and the pauses;
- ``ambit_point``: every AMBIT sample (time, s_630, r_630, ratio, actinic, role);
- ``felix_window``: every Felix histogram window, rebinned to 47 ps (``REBIN`` x 2.34 ps; 266
  bins per 5 ms window, ~0.8 MB per run), optional.

Every row carries ``run_id`` (the run folder's name) and an index, so duplicates can be told.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime
from pathlib import Path
from typing import Iterator, Optional

import numpy as np

from .openjii import OpenJIIClient, OpenJIIError, UploadRecord

REBIN = 20                       # 20 x 2.34 ps = 47 ps bins, as the fits use
UPLOAD_FILE = "openjii_upload.ndjson"
MAX_FILE_BYTES = 100_000_000     # platform limit per NDJSON file


def _plain(value):
    """JSON-safe: numpy scalars and arrays to Python, nan/inf to null (strict JSON has no NaN)."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        v = float(value)
        return v if math.isfinite(v) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def _timestamp(meta: dict, run_dir: Path) -> str:
    """The run's save time as local ISO 8601 (from run.json, else the folder's mtime)."""
    try:
        t = datetime.strptime(meta["saved"], "%Y-%m-%d %H:%M:%S")
    except (KeyError, ValueError):
        t = datetime.fromtimestamp(run_dir.stat().st_mtime)
    return t.astimezone().isoformat(timespec="seconds")


def run_rows(res: dict, summary: dict, include_felix: bool = True, rebin: int = REBIN,
             label: str = "") -> Iterator[dict]:
    """The rows of one analysed run: ``res`` from ``analyse_pauses``, ``summary`` from ``run_summary``."""
    run_dir = Path(res["run_dir"])
    run_id = run_dir.name
    meta = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    ts = _timestamp(meta, run_dir)
    roles = {int(r["line"]): r["role"] for r in res["lines"]}
    results = {k: v for k, v in summary.items() if k not in ("fits", "warnings")}
    yield _plain(dict(kind="run", run_id=run_id, timestamp_local=ts, label=label, note=meta.get("note", ""),
                      protocol=meta.get("protocol"), ambit_cmd=meta.get("cmd"), laser_hz=meta.get("laser_hz"),
                      window_s=res["W"], bin_s=res["bin_s"], ambit=meta.get("ambit"), felix=meta.get("felix"),
                      pre_illumination=meta.get("pre_illumination"), alignment=meta.get("align"),
                      results=results, fits=summary["fits"], warnings=summary["warnings"],
                      felix_windows_included=bool(include_felix), felix_rebin=int(rebin) if include_felix else None))
    a = res["ambit"]
    for i in range(a["t"].size):
        yield _plain(dict(kind="ambit_point", run_id=run_id, timestamp_local=ts, index=i, t_s=a["t"][i],
                          s_630=a["s_630"][i], r_630=a["r_630"][i], ratio=a["ratio"][i], actinic=a["actinic"][i],
                          line=a["line"][i], role=roles.get(int(a["line"][i]), "")))
    if include_felix:
        bins = np.asarray(res["bins"])
        n = bins.shape[1] // rebin
        reb = bins[:, :n * rebin].reshape(bins.shape[0], n, rebin).sum(1 + 1).astype(np.int64)
        bin_ps = res["bin_s"] * rebin * 1e12
        for i in range(bins.shape[0]):
            yield _plain(dict(kind="felix_window", run_id=run_id, timestamp_local=ts, index=i,
                              t_start_s=round(float(res["ts"][i]), 6), window_s=res["W"], bin_ps=bin_ps,
                              n_bins=n, counts=reb[i].tolist()))


def write_ndjson(run_dir: Path, rows, name: str = UPLOAD_FILE) -> tuple[Path, int]:
    """Write the rows as NDJSON into the run folder; returns (path, number of rows). Raises
    OpenJIIError when the file would exceed the platform's size limit."""
    path = Path(run_dir) / name
    n = 0
    with open(path, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, separators=(",", ":"), allow_nan=False) + "\n")
            n += 1
    size = path.stat().st_size
    if size > MAX_FILE_BYTES:
        raise OpenJIIError(f"{path.name} is {size / 1e6:.0f} MB, above the platform's {MAX_FILE_BYTES / 1e6:.0f} MB "
                           "per file: upload without the Felix windows")
    return path, n


def upload_and_wait(client: OpenJIIClient, experiment_id: str, table_name: str, path: Path, log=print,
                    poll_s: float = 10.0, timeout_s: float = 900.0) -> UploadRecord:
    """Upload one NDJSON file into the experiment's table (created on the first upload), then poll
    until the platform reports the rows queryable or the upload failed. Returns the final record;
    a record with status 'uploaded' means the platform had not finished within ``timeout_s``."""
    table_id = client.find_table(experiment_id, table_name)
    log(f"upload: {'appending to' if table_id else 'creating'} table {table_name} "
        f"({path.stat().st_size / 1e6:.1f} MB)")
    result = client.upload_files(experiment_id, table_name, table_id, [path])
    upload_id = result["uploadId"]
    accepted = {f.get("fileName") for f in result.get("files") or []}
    if path.name not in accepted:
        return UploadRecord(upload_id, "failed", result.get("uploadTableId"), table_name, None,
                            "file_not_accepted", [path.name])
    log(f"upload {upload_id[:8]} sent; waiting for the platform (job on Databricks, polled every {poll_s:.0f} s)")
    t_end = time.monotonic() + timeout_s
    last = ""
    while time.monotonic() < t_end:
        time.sleep(poll_s)
        rec = client.upload_status(experiment_id, table_name, {upload_id}).get(upload_id)
        if rec is None:
            continue
        if rec.status != last:
            log(f"upload {upload_id[:8]}: {rec.status}")
            last = rec.status
        if rec.finished:
            return rec
    return UploadRecord(upload_id, "uploaded", result.get("uploadTableId"), table_name, None, None)


def describe(rec: UploadRecord) -> str:
    if rec.status == "completed":
        return f"completed: {rec.row_count} rows queryable in table {rec.table_name or ''}"
    if rec.status == "failed":
        return f"FAILED on the platform: {rec.error_code or '?'} {rec.rejected_files or ''}"
    return f"{rec.status}: not confirmed yet; check the experiment's Data page later"
