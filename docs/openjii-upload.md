# openJII upload

How the phi2 GUI puts a run on [openJII](https://openjii.org), what the platform expects, and how to
change the upload without breaking someone's table. The user-facing summary is in the README
("Uploading to openJII"); the code is `scripts/ambit-gate/phi2gui/upload.py` and `openjii.py`.

## One run, one row

`upload.run_row` builds one JSON object per run, and `write_run` writes it as a one-line NDJSON file
(`openjii_upload.ndjson`) into the run folder. The platform keeps nested values, so the row is a
structure, queryable in SQL as `ambit_points.ratio[i]` or `felix_windows.counts[j]`:

| Key | Holds |
| --- | --- |
| `row_format` | `upload.ROW_FORMAT`, currently 2 |
| `run_id`, `timestamp_local`, `label`, `note` | Identity: `run_id` is the run folder's name and tells duplicates apart |
| `protocol`, `ambit_cmd`, `laser_hz`, `window_s`, `bin_s`, `ambit`, `felix`, `pre_illumination`, `alignment` | From `run.json` and the analysis |
| `results`, `fits`, `warnings` | `ambit_pause.run_summary`: phi2 from the AMBIT and from the lifetimes |
| `ambit_points` | One array per AMBIT column: `t_s`, `s_630`, `r_630`, `ratio`, `actinic`, `line`, `role` |
| `felix_windows` | `t_start_s` and `counts` (one array per window), with `rebin`, `bin_ps`, `n_bins`, `included` |

NaN and infinity become `null` (`upload._plain`), since strict JSON has neither.

**Changing the row's shape** (renaming, retyping or removing a key) makes it incompatible with rows
already in a table. Raise `ROW_FORMAT` and change `DEFAULT_TABLE` in `app.py`; the GUI then switches
users to the new default table name once (`_load_settings`). Adding a key is fine within one format.

### Felix windows

| Choice (`upload.FELIX_UPLOAD`) | `rebin` | Bins per 5 ms window | 5 s run |
| --- | --- | --- | --- |
| `rebinned` (default) | 20 | 266 (46.875 ps) | ~0.8 MB |
| `full` | 1 | ~5333 (2.34375 ps) | ~13 MB |
| `none` | - | - | ~20 kB |

A full-resolution row over `MAX_ROW_BYTES` (16 MiB) is written with 47 ps bins instead
(`write_run`), and the log says so.

## The platform's side

Endpoints, all under `/api/v1` with the user's personal API key in `x-api-key`:

| Call | Used for |
| --- | --- |
| `GET /auth/get-session` | Checking the key (200 with a null body means rejected) |
| `GET /experiments?filter=member&status=active` | The experiment list |
| `GET /experiments/{id}/access` | `canContribute`, `canManage` |
| `GET /experiments/{id}/data/uploads?uploadTableName=` | Finding the table id, and polling an upload's status |
| `POST /experiments/{id}/data/uploads` (multipart) | The upload: `sourceKind=ndjson`, `uploadId`, then `targetKind=new` + `targetName`, or `existing` + `uploadTableId` |

Rules the client already handles, and the reason for each:

- Creating a table needs **manage** access (experiment admin); appending needs **contribute**.
  `check_access` refuses before sending, with the same message as the platform's 403.
- A new table appears in the platform's table list a few minutes after its first upload finished.
  An append before that gets 404 "No upload table with id"; the GUI says to retry.
- Each upload starts a Databricks job. `upload_and_wait` polls every 10 s for up to 15 minutes; a
  full-resolution row took about 3.5 minutes on dev.
- Limits: 100 MB per NDJSON file, 100 files per upload.

Checked on openJII dev on 2026-11-06 with mock runs built from a real prod run:

- A 12.9 MB full-resolution row uploads, completes, and reads back from SQL with every count
  identical to the file.
- openJII's Data page cannot list a table holding such a row: the API answers `PAGE_TOO_LARGE`
  even for a page of one row, because the warehouse sizes each count as an 8-byte integer (about
  46 MB) against a 25 MiB limit per read. Tracked as openJII ticket OJD-2658. Hence the 47 ps
  default and the advice to keep full-resolution runs in a table of their own.
- A table of 11 runs at 47 ps opens at the Data page's default 10 rows per page.

## Testing an upload change

Test against **dev** (`https://dev.openjii.org`, API `https://api.dev.openjii.org`), never with a
production key. The GUI's environment selector switches between the two; `openjii.ENVIRONMENTS`
holds the URLs.

1. Create a dev API key in the dev web app (Platform > Account > API keys) and a private experiment
   named with a `QA-<date>-` prefix, so test data is recognisable.
2. Without hardware, build `res` and `summary` from an existing run folder
   (`ambit_pause.analyse_pauses(run_dir, irf)` and `run_summary(res)`), or from a row already on the
   platform, then call `upload.write_run` and `upload.upload_and_wait` with an `OpenJIIClient` for
   `ENVIRONMENTS["dev"]`.
3. Check the outcome three ways: the returned record says `completed`, SQL reads the row back, and
   the experiment's Data page opens the table.

`tests/test_upload_row.py` and `tests/test_upload_access.py` cover the row layout, the size
fallback and the access check without the network; extend them with any change to `upload.py`.
