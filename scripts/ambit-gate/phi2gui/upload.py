"""One pauses run as ONE row of an openJII upload table, the NDJSON file, and the upload itself.

The row (``run_row``) holds the run's metadata (protocol, instruments, note), the phi2 results
(``ambit_pause.run_summary``), the fit numbers of the reference decay and the pauses, and the raw
traces as nested arrays: ``ambit_points`` (one array per AMBIT column: t_s, s_630, r_630, ratio,
actinic, line, role) and ``felix_windows`` (t_start_s and the counts of every Felix histogram
window). The Felix windows go up rebinned to 47 ps (``REBIN`` x 2.34 ps, 266 bins per 5 ms
window, ~0.8 MB for a 5 s run), at full resolution (2.34 ps, ~5300 bins per window, ~13 MB), or
not at all (``FELIX_UPLOAD``); ``felix_windows.rebin`` and ``bin_ps`` say which. A
full-resolution row above ``MAX_ROW_BYTES`` falls back to 47 ps (``write_run``). openJII stores a
13 MB row whole (checked on dev: every count back from SQL), but its Data page cannot list a table
holding one yet, so full-resolution runs belong in a table of their own. The platform keeps nested values as STRUCT/ARRAY columns, so a table row is a
measurement and ``ambit_points.ratio[i]`` or ``felix_windows.counts[j]`` is reachable in SQL.
``run_id`` (the run folder's name) tells duplicates.

Rows before ``ROW_FORMAT`` 2 were one per AMBIT sample and per Felix window (``kind`` column);
those tables cannot take the new rows, hence a new default table name in the GUI.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from .openjii import OpenJIIClient, OpenJIIError, UploadRecord

ROW_FORMAT = 2                   # one row per run (1: one row per sample and per window)
REBIN = 20                       # 20 x 2.34 ps = 47 ps bins, as the fits use
FELIX_UPLOAD = {"full": 1, "rebinned": REBIN, "none": None}   # choice -> rebin of the Felix windows
UPLOAD_FILE = "openjii_upload.ndjson"
MAX_FILE_BYTES = 100_000_000     # platform limit per NDJSON file
MAX_ROW_BYTES = 16 * 1024 * 1024  # one row stays below this (a 5 s run at full resolution is ~13 MB)


class RowTooLarge(OpenJIIError):
    """The run's row is above ``MAX_ROW_BYTES``."""


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


def run_row(res: dict, summary: dict, include_felix: bool = True, rebin: int = REBIN, label: str = "") -> dict:
    """The one row of an analysed run: ``res`` from ``analyse_pauses``, ``summary`` from ``run_summary``."""
    run_dir = Path(res["run_dir"])
    meta = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    roles = {int(r["line"]): r["role"] for r in res["lines"]}
    results = {k: v for k, v in summary.items() if k not in ("fits", "warnings")}
    a = res["ambit"]
    lines = np.asarray(a["line"]).astype(int)
    ambit_points = dict(n=int(lines.size), t_s=a["t"], s_630=a["s_630"], r_630=a["r_630"], ratio=a["ratio"],
                        actinic=a["actinic"], line=lines, role=[roles.get(int(ln), "") for ln in lines])
    felix_windows = dict(included=bool(include_felix), n=None, rebin=None, bin_ps=None, n_bins=None,
                         t_start_s=None, counts=None)
    if include_felix:
        bins = np.asarray(res["bins"])
        n = bins.shape[1] // rebin
        reb = bins[:, :n * rebin].reshape(bins.shape[0], n, rebin).sum(2).astype(np.int64)
        felix_windows.update(n=int(bins.shape[0]), rebin=int(rebin), bin_ps=res["bin_s"] * rebin * 1e12, n_bins=n,
                             t_start_s=[round(float(t), 6) for t in res["ts"]], counts=reb)
    return _plain(dict(row_format=ROW_FORMAT, run_id=run_dir.name, timestamp_local=_timestamp(meta, run_dir),
                       label=label, note=meta.get("note", ""), protocol=meta.get("protocol"), ambit_cmd=meta.get("cmd"),
                       laser_hz=meta.get("laser_hz"), window_s=res["W"], bin_s=res["bin_s"], ambit=meta.get("ambit"),
                       felix=meta.get("felix"), pre_illumination=meta.get("pre_illumination"),
                       alignment=meta.get("align"), results=results, fits=summary["fits"],
                       warnings=summary["warnings"], ambit_points=ambit_points, felix_windows=felix_windows))


def write_ndjson(run_dir: Path, rows, name: str = UPLOAD_FILE) -> tuple[Path, int]:
    """Write the rows as NDJSON into the run folder; returns (path, number of rows). Raises
    OpenJIIError when a row or the file would exceed the platform's size limits."""
    path = Path(run_dir) / name
    n = 0
    with open(path, "w", encoding="utf-8") as fh:
        for row in rows:
            line = json.dumps(row, separators=(",", ":"), allow_nan=False)
            if len(line) > MAX_ROW_BYTES:
                raise RowTooLarge(f"the run's row is {len(line) / 1e6:.0f} MB of JSON, above the "
                                  f"{MAX_ROW_BYTES / 2 ** 20:.0f} MiB per row: upload the Felix windows in "
                                  "47 ps bins or without the Felix windows")
            fh.write(line + "\n")
            n += 1
    size = path.stat().st_size
    if size > MAX_FILE_BYTES:
        raise OpenJIIError(f"{path.name} is {size / 1e6:.0f} MB, above the platform's {MAX_FILE_BYTES / 1e6:.0f} MB "
                           "per file: upload without the Felix windows")
    return path, n


def write_run(run_dir: Path, res: dict, summary: dict, felix: str = "rebinned", label: str = "",
              log=print) -> Path:
    """The run's one-row NDJSON file, with the Felix windows as ``felix`` (a ``FELIX_UPLOAD`` key)
    says. A full-resolution row above ``MAX_ROW_BYTES`` (a run much longer than 5 s) is written
    with 47 ps bins instead, and the log says so."""
    rebin = FELIX_UPLOAD[felix]
    try:
        path, _ = write_ndjson(run_dir, [run_row(res, summary, include_felix=rebin is not None,
                                                 rebin=rebin or REBIN, label=label)])
    except RowTooLarge as exc:
        if rebin != 1:
            raise
        log(f"upload: {str(exc).split(':')[0]}; uploading the Felix windows in 47 ps bins instead")
        path, _ = write_ndjson(run_dir, [run_row(res, summary, rebin=REBIN, label=label)])
    return path


NO_TABLE_TEXT = ("there is no upload table '{table}' in {exp} yet, and creating one needs manage access "
                 "(experiment admin), which your account does not have. Ask an admin of the experiment "
                 "to upload the first run from this GUI with the same table name (later uploads only "
                 "need contribute access), or to make you an admin; or type the name of an upload table "
                 "that already exists (Data page of the experiment)")


def check_access(client: OpenJIIClient, experiment_id: str, table_name: str,
                 experiment_name: str = "") -> Optional[str]:
    """The id of the existing upload table ``table_name`` (None when the first upload will create
    it), or OpenJIIError when this account cannot upload: appending needs *contribute* access on
    the experiment, creating a table *manage* access (experiment admin). The platform refuses a
    creating upload by a contributor with 403 "Creating a new data table requires experiment
    management access"; this says the same before the file is sent."""
    exp = f"'{experiment_name}'" if experiment_name else "this experiment"
    caps = client.capabilities(experiment_id)
    if not caps.get("canContribute"):
        raise OpenJIIError(f"your account cannot contribute data to {exp}: ask its admin to add you "
                           "as a member")
    table_id = client.find_table(experiment_id, table_name)
    if table_id is None and not caps.get("canManage"):
        raise OpenJIIError(NO_TABLE_TEXT.format(table=table_name, exp=exp))
    return table_id


def upload_and_wait(client: OpenJIIClient, experiment_id: str, table_name: str, path: Path, log=print,
                    poll_s: float = 10.0, timeout_s: float = 900.0, experiment_name: str = "") -> UploadRecord:
    """Upload one NDJSON file into the experiment's table (created on the first upload), then poll
    until the platform reports the rows queryable or the upload failed. Returns the final record;
    a record with status 'uploaded' means the platform had not finished within ``timeout_s``.
    OpenJIIError before anything is sent when the account cannot upload (``check_access``)."""
    table_id = check_access(client, experiment_id, table_name, experiment_name)
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
