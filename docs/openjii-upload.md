# openJII upload

How the phi2 GUI puts a run on [openJII](https://openjii.org), what the platform expects, and how to
change the upload without breaking someone's table. The README's "Upload to openJII" section is the
summary for users. The code lives in `scripts/ambit-gate/phi2gui/upload.py` and `openjii.py`.

## One run, one row

`upload.run_row` builds one JSON object per run, and `write_run` writes it as a one-line NDJSON file,
`openjii_upload.ndjson`, in the run folder. The platform keeps nested values, so SQL can reach into
the row with `ambit_points.ratio[i]` or `felix_windows.counts[j]`.

| Key | Holds |
| --- | --- |
| `row_format` | `upload.ROW_FORMAT`, now 2 |
| `run_id`, `timestamp_local`, `label`, `note` | The run's identity. `run_id` is the run folder's name and tells duplicates apart. |
| `protocol`, `ambit_cmd`, `laser_hz`, `window_s`, `bin_s`, `ambit`, `felix`, `pre_illumination`, `alignment` | Values from `run.json` and the analysis |
| `results`, `fits`, `warnings` | `ambit_pause.run_summary`, with phi2 from the AMBIT and from the lifetimes |
| `ambit_points` | One array per AMBIT column: `t_s`, `s_630`, `r_630`, `ratio`, `actinic`, `line`, `role` |
| `felix_windows` | `t_start_s` and `counts`, one array per window, with `rebin`, `bin_ps`, `n_bins` and `included` |

`upload._plain` turns NaN and infinity into `null`, since strict JSON has neither.

Changing the row's shape breaks existing tables. That means renaming a key, changing its type or
removing it. In that case raise `ROW_FORMAT` and change `DEFAULT_TABLE` in `app.py`, and
`_load_settings` moves users to the new default table name once. Adding a key is safe within a
format.

### Felix windows

| Choice in `upload.FELIX_UPLOAD` | `rebin` | Bins per 5 ms window | 5 s run |
| --- | --- | --- | --- |
| `rebinned`, the default | 20 | 266 of 46.875 ps | about 0.8 MB |
| `full` | 1 | about 5333 of 2.34375 ps | about 13 MB |
| `none` | none | none | about 20 kB |

`write_run` writes a full-resolution row larger than `MAX_ROW_BYTES`, 16 MiB, with 47 ps bins
instead, and the log says so.

## The platform's side

Every call goes to `/api/v1` with the user's personal API key in the `x-api-key` header.

| Call | Used for |
| --- | --- |
| `GET /auth/get-session` | Checking the key. A 200 with a null body means openJII rejected it. |
| `GET /experiments?filter=member&status=active` | The experiment list |
| `GET /experiments/{id}/access` | Reading `canContribute` and `canManage` |
| `GET /experiments/{id}/data/uploads?uploadTableName=` | Finding the table's id and polling an upload's status |
| `POST /experiments/{id}/data/uploads`, multipart | The upload itself, with `sourceKind=ndjson` and `uploadId`, then either `targetKind=new` and `targetName`, or `existing` and `uploadTableId` |

The client already follows these rules:

- Creating a table needs manage access, which experiment admins have, and appending needs
  contribute access. `check_access` refuses before sending, with the same message as the
  platform's 403.
- openJII lists a new table a few minutes after its first upload finishes. An append before then
  gets 404 "No upload table with id", and the GUI tells the user to retry.
- Each upload starts a Databricks job. `upload_and_wait` polls every 10 s for up to 15 minutes. A
  full-resolution row took about 3.5 minutes on dev.
- openJII takes at most 100 MB per NDJSON file and 100 files per upload.

On 6 November 2026 I tested on openJII dev with mock runs built from a real production run:

- A 12.9 MB full-resolution row uploaded, completed, and read back from SQL with every count equal
  to the file.
- openJII's Data page cannot list a table that holds such a row. The API answers `PAGE_TOO_LARGE`
  even for a page of one row. The warehouse counts every number as an 8-byte integer, about 46 MB,
  against a 25 MiB limit per read. openJII tracks this as OJD-2658. That is why 47 ps is the
  default and why full-resolution runs belong in a table of their own.
- A table of 11 runs at 47 ps opens at the Data page's default of 10 rows per page.

## Testing an upload change

Test on openJII dev, at `https://dev.openjii.org` with the API at `https://api.dev.openjii.org`,
and never with a production key. The GUI's environment selector switches between the two, and
`openjii.ENVIRONMENTS` holds the URLs.

1. Create a dev API key in the dev web app under Platform > Account > API keys, and a private
   experiment whose name starts with `QA-<date>-`, so test data is easy to spot.
2. Without hardware, build `res` and `summary` from a saved run folder with
   `ambit_pause.analyse_pauses(run_dir, irf)` and `run_summary(res)`, or from a row already on the
   platform. Then call `upload.write_run` and `upload.upload_and_wait` with an `OpenJIIClient` for
   `ENVIRONMENTS["dev"]`.
3. Check the outcome three ways. The returned record says `completed`, SQL reads the row back, and
   the experiment's Data page opens the table.

`tests/test_upload_row.py` and `tests/test_upload_access.py` cover the row layout, the size fallback
and the access check without the network. Extend them with any change to `upload.py`.
