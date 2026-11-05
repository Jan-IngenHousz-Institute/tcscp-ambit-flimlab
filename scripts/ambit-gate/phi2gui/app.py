"""The window: instruments, openJII, the measurement, and the results of the current run.

Everything that waits on hardware or the network runs in a ``Task`` thread and reports through
signals; the window stays responsive and the buttons that would collide are disabled meanwhile.
The analysis modules draw with pyplot on the Agg backend; their figures are embedded by hand.
"""
from __future__ import annotations

import json
import shutil
import sys
import time
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                                   # before the analysis modules import pyplot
import matplotlib.pyplot as plt                         # noqa: E402
import numpy as np                                      # noqa: E402
import serial.tools.list_ports                          # noqa: E402
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT   # noqa: E402
from PySide6 import QtCore, QtGui, QtWidgets            # noqa: E402

import ambit_gate as ag                                 # noqa: E402
import ambit_pause as ap                                # noqa: E402
import blank_dose as bd                                 # noqa: E402

from . import driver as drv                             # noqa: E402
from . import openjii as oj                             # noqa: E402
from . import upload as up                              # noqa: E402

HERE = Path(__file__).resolve().parent
if getattr(sys, "frozen", False):                       # PyInstaller build: runs/ next to phi2gui.exe
    RUNS_DIR = Path(sys.executable).resolve().parent / "runs"
else:
    RUNS_DIR = HERE.parent / "runs"
DEFAULT_TABLE = "ambit_tcspc_pauses"
HIGHLIGHT = QtGui.QColor("#fff0a8")
RED = "#b00000"


class Task(QtCore.QThread):
    """Runs ``fn(log)`` in a thread; ``log`` posts text to the window, ``done`` carries the result."""
    message = QtCore.Signal(str)
    done = QtCore.Signal(object)
    failed = QtCore.Signal(str)

    def __init__(self, fn, parent=None) -> None:
        super().__init__(parent)
        self._fn = fn

    def run(self) -> None:
        try:
            self.done.emit(self._fn(self.message.emit))
        except Exception as exc:                        # reported to the user, never lost
            self.failed.emit(f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=4)}")


def _spin(lo, hi, value, decimals: int = 0, step=1, suffix: str = ""):
    if decimals:
        s = QtWidgets.QDoubleSpinBox()
        s.setDecimals(decimals)
        s.setSingleStep(step)
    else:
        s = QtWidgets.QSpinBox()
        s.setSingleStep(int(step))
    s.setRange(lo, hi)
    s.setValue(value)
    if suffix:
        s.setSuffix(suffix)
    return s


def _num(v, p: int = 3) -> str:
    return "-" if v is None or not np.isfinite(v) else f"{v:.{p}f}"


def _pm(v, e, p: int = 3) -> str:
    return _num(v, p) if e is None or not np.isfinite(e) else f"{_num(v, p)} ± {_num(e, p)}"


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, settings: oj.Settings) -> None:
        super().__init__()
        self.setWindowTitle("phi2: AMBIT + Felix")
        self.settings = settings
        self.ambit = None
        self.felix = None
        self.client: oj.OpenJIIClient | None = None
        self.experiments: list[oj.Experiment] = []
        self.irf: dict | None = None
        self.irf_dir = oj.config_dir() / "irf-library"        # the user's own IRF files (Add IRF file...)
        self.session = RUNS_DIR / (time.strftime("%Y%m%d_%H%M%S") + "_gui")
        self.run_dir: Path | None = None
        self.res: dict | None = None
        self.summary: dict | None = None
        self._tasks: list[Task] = []
        self._build_ui()
        self._load_settings()
        self._rescan_ports()
        self._refresh_irfs(select=str(self.settings.get("irf", "")))
        QtCore.QTimer.singleShot(800, self._check_driver_at_start)     # once the window is up

    # ------------------------------------------------------------------ ui
    def _build_ui(self) -> None:
        left = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(left)
        lay.setContentsMargins(6, 6, 6, 6)

        # instruments
        g = QtWidgets.QGroupBox("Instruments")
        f = QtWidgets.QFormLayout(g)
        self.port_combo = QtWidgets.QComboBox()
        self.btn_rescan = QtWidgets.QPushButton("Rescan")
        self.btn_rescan.clicked.connect(self._rescan_ports)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.port_combo, 1)
        row.addWidget(self.btn_rescan)
        f.addRow("AMBIT port", row)
        self.btn_ambit = QtWidgets.QPushButton("Connect AMBIT")
        self.btn_ambit.clicked.connect(self._toggle_ambit)
        self.lbl_ambit = QtWidgets.QLabel("not connected")
        self.lbl_ambit.setWordWrap(True)
        f.addRow(self.btn_ambit, self.lbl_ambit)
        self.btn_felix = QtWidgets.QPushButton("Connect Felix")
        self.btn_felix.clicked.connect(self._toggle_felix)
        self.lbl_felix = QtWidgets.QLabel("not connected")
        self.lbl_felix.setWordWrap(True)
        f.addRow(self.btn_felix, self.lbl_felix)
        self.spin_threshold = _spin(0.0, 2.5, ag.THRESHOLD_V, 2, 0.1, " V")
        f.addRow("Felix thresholds", self.spin_threshold)
        self.btn_rates = QtWidgets.QPushButton("Read count rates")
        self.btn_rates.clicked.connect(self._read_rates)
        self.lbl_rates = QtWidgets.QLabel("-")
        f.addRow(self.btn_rates, self.lbl_rates)
        self.irf_combo = QtWidgets.QComboBox()
        self.irf_combo.setToolTip("Instrument response for the reconvolution fits: the shipped default, "
                                  "or a TDC-Software histogram export you added")
        self.irf_combo.currentIndexChanged.connect(self._irf_changed)
        self.btn_irf = QtWidgets.QPushButton("Add IRF file...")
        self.btn_irf.clicked.connect(self._add_irf)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.irf_combo, 1)
        row.addWidget(self.btn_irf)
        f.addRow("IRF", row)
        self.lbl_irf = QtWidgets.QLabel("")
        self.lbl_irf.setWordWrap(True)
        f.addRow("", self.lbl_irf)
        lay.addWidget(g)

        # openJII
        g = QtWidgets.QGroupBox("openJII")
        f = QtWidgets.QFormLayout(g)
        self.env_combo = QtWidgets.QComboBox()
        for key in oj.ENVIRONMENTS:
            self.env_combo.addItem(key, key)
        self.env_combo.currentIndexChanged.connect(self._env_changed)
        f.addRow("Environment", self.env_combo)
        self.key_edit = QtWidgets.QLineEdit()
        self.key_edit.setEchoMode(QtWidgets.QLineEdit.EchoMode.Password)
        self.key_edit.setPlaceholderText("jii_... (Platform > Account > API keys)")
        f.addRow("API key", self.key_edit)
        self.btn_signin = QtWidgets.QPushButton("Sign in")
        self.btn_signin.clicked.connect(self._sign_in)
        self.lbl_user = QtWidgets.QLabel("not signed in")
        self.lbl_user.setWordWrap(True)
        f.addRow(self.btn_signin, self.lbl_user)
        self.exp_combo = QtWidgets.QComboBox()
        self.btn_exps = QtWidgets.QPushButton("Refresh")
        self.btn_exps.clicked.connect(self._refresh_experiments)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.exp_combo, 1)
        row.addWidget(self.btn_exps)
        f.addRow("Experiment", row)
        self.table_edit = QtWidgets.QLineEdit(DEFAULT_TABLE)
        f.addRow("Upload table", self.table_edit)
        self.chk_upload = QtWidgets.QCheckBox("upload every run when it is analysed")
        f.addRow("", self.chk_upload)
        self.chk_felix = QtWidgets.QCheckBox("include the Felix histograms (47 ps bins, ~1 MB per run)")
        self.chk_felix.setChecked(True)
        f.addRow("", self.chk_felix)
        self.btn_upload = QtWidgets.QPushButton("Upload this run")
        self.btn_upload.clicked.connect(self._upload_current)
        self.lbl_upload = QtWidgets.QLabel("-")
        self.lbl_upload.setWordWrap(True)
        f.addRow(self.btn_upload, self.lbl_upload)
        lay.addWidget(g)

        # measurement
        g = QtWidgets.QGroupBox("Measurement (section 1 of notebook 4)")
        f = QtWidgets.QFormLayout(g)
        self.spin_level = _spin(ag.ACTINIC_ON + 1, 255, 255)
        f.addRow("Actinic byte of the lit blocks", self.spin_level)
        self.spin_dark = _spin(1, ag.MAX_POINTS, 100)
        f.addRow("Dark points before and after", self.spin_dark)
        self.spin_block = _spin(0.02, 5.0, 0.2, 2, 0.05, " s")
        f.addRow("Lit block", self.spin_block)
        self.pauses_edit = QtWidgets.QLineEdit("2,2,2,2")
        self.pauses_edit.setToolTip("AMBIT points per pause, one entry per pause (1 point = ~10 ms dark at 100 Hz; "
                                    "the Felix window must be at most half of a pause)")
        f.addRow("Pauses (points each)", self.pauses_edit)
        self.spin_freq = _spin(1, 2000, 100, suffix=" Hz")
        f.addRow("AMBIT sampling", self.spin_freq)
        self.spin_window = _spin(1.0, 100.0, 5.0, 1, 1.0, " ms")
        f.addRow("Felix window", self.spin_window)
        self.spin_fo = _spin(1, ag.MAX_POINTS, ap.FO_POINTS)
        f.addRow("Fo: last dark points", self.spin_fo)
        self.spin_fm = _spin(1, ag.MAX_POINTS, ap.FM_POINTS)
        f.addRow("Fm: last lit points", self.spin_fm)
        self.label_edit = QtWidgets.QLineEdit("leaf")
        f.addRow("Label", self.label_edit)
        self.note_edit = QtWidgets.QLineEdit()
        self.note_edit.setPlaceholderText("laser power, ND, leaf ...")
        f.addRow("Note", self.note_edit)
        self.lbl_protocol = QtWidgets.QLabel("")
        self.lbl_protocol.setWordWrap(True)
        f.addRow("", self.lbl_protocol)
        for w in (self.spin_level, self.spin_dark, self.spin_block, self.spin_freq):
            w.valueChanged.connect(self._update_protocol_label)
        self.pauses_edit.textChanged.connect(self._update_protocol_label)
        self.btn_run = QtWidgets.QPushButton("Run measurement")
        self.btn_run.setStyleSheet("font-weight: bold")
        self.btn_run.clicked.connect(self._run_measurement)
        self.btn_open = QtWidgets.QPushButton("Open saved run...")
        self.btn_open.clicked.connect(self._open_run_dialog)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.btn_run, 1)
        row.addWidget(self.btn_open)
        f.addRow(row)
        lay.addWidget(g)

        self.log_box = QtWidgets.QPlainTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumBlockCount(3000)
        self.log_box.setFont(QtGui.QFont("Consolas", 8))
        lay.addWidget(self.log_box, 1)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(left)
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(430)

        # results
        self.tabs = QtWidgets.QTabWidget()
        res_w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(res_w)
        self.lbl_title = QtWidgets.QLabel("no run yet: connect the instruments and run a measurement, or open a saved run")
        self.lbl_title.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.lbl_title.setWordWrap(True)
        v.addWidget(self.lbl_title)
        self.lbl_summary = QtWidgets.QLabel("")
        self.lbl_summary.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.lbl_summary.setWordWrap(True)
        self.lbl_summary.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.lbl_summary)
        self.lbl_warnings = QtWidgets.QLabel("")
        self.lbl_warnings.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.lbl_warnings.setWordWrap(True)
        v.addWidget(self.lbl_warnings)
        self.table = QtWidgets.QTableWidget(0, 0)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        v.addWidget(self.table, 1)
        self.lbl_table_note = QtWidgets.QLabel("")
        self.lbl_table_note.setWordWrap(True)
        self.lbl_table_note.setStyleSheet("color: gray")
        v.addWidget(self.lbl_table_note)
        self.tabs.addTab(res_w, "Results")

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        split.addWidget(scroll)
        split.addWidget(self.tabs)
        split.setStretchFactor(1, 1)
        self.setCentralWidget(split)
        self.statusBar()
        self._update_protocol_label()

    # --------------------------------------------------------------- helpers
    def _log(self, text: str) -> None:
        for line in str(text).rstrip().splitlines():
            self.log_box.appendPlainText(time.strftime("%H:%M:%S ") + line)
        self.statusBar().showMessage(str(text).splitlines()[0] if text else "", 10000)

    def _start_task(self, fn, on_done, on_failed=None) -> Task:
        self._tasks = [t for t in self._tasks if t.isRunning()]
        task = Task(fn, self)
        task.message.connect(self._log)
        task.done.connect(on_done)
        task.failed.connect(on_failed or self._task_failed)
        self._tasks.append(task)
        task.start()
        return task

    def _task_failed(self, text: str) -> None:
        self._log("ERROR: " + text)
        self._busy(False)

    def _busy(self, busy: bool) -> None:
        for w in (self.btn_run, self.btn_open, self.btn_upload, self.btn_ambit, self.btn_felix, self.btn_rates):
            w.setEnabled(not busy)

    def _env(self) -> oj.Environment:
        return oj.ENVIRONMENTS[self.env_combo.currentData()]

    # ------------------------------------------------------------------- IRF
    def _refresh_irfs(self, select: str = "") -> None:
        """List the IRF files (the user's own folder first, then the shipped ones) and load
        ``select`` (a file name) or the shipped default."""
        files = bd.list_irfs([self.irf_dir])
        self.irf_combo.blockSignals(True)
        self.irf_combo.clear()
        for p in files:
            self.irf_combo.addItem(p.name, str(p))
        i = self.irf_combo.findText(select) if select else -1
        if i < 0:
            i = self.irf_combo.findText(bd.DEFAULT_IRF_NAME)
        self.irf_combo.setCurrentIndex(max(i, 0))
        self.irf_combo.blockSignals(False)
        self._irf_changed()

    def _irf_changed(self) -> None:
        path = self.irf_combo.currentData()
        try:
            self.irf = bd.load_irf(path) if path else None
        except Exception as exc:
            self.irf = None
            self._log(f"ERROR: could not read the IRF {path}: {exc}")
        if self.irf is None:
            self.lbl_irf.setText("NO IRF: the decays will not be fitted, no lifetime phi2. "
                                 "Add one with 'Add IRF file...' (a TDC-Software histogram export).")
            self.lbl_irf.setStyleSheet(f"color: {RED}; font-weight: bold")
            self._log("WARNING: no IRF: the decays will not be fitted, no lifetime phi2")
            return
        self.lbl_irf.setText(f"{self.irf['name']}: {self.irf['bin_s'] * 1e12:.2f} ps bins, "
                             f"{self.irf['counts'].sum():,.0f} counts")
        self.lbl_irf.setStyleSheet("color: gray")
        self.settings.set("irf", self.irf["name"])
        self._log(f"IRF: {self.irf['name']}")

    def _add_irf(self) -> None:
        """Copy a TDC-Software histogram export into the user's IRF folder and select it."""
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "IRF file (TDC-Software histogram export)", "",
                                                        "CSV files (*.csv);;All files (*)")
        if not path:
            return
        src = Path(path)
        try:
            bd.load_irf(src)                                 # readable as an IRF?
            self.irf_dir.mkdir(parents=True, exist_ok=True)
            dst = self.irf_dir / src.name
            if src.resolve() != dst.resolve():
                shutil.copy2(src, dst)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "IRF", f"{src.name} could not be used as an IRF:\n{exc}")
            return
        self._log(f"IRF added: {dst}")
        self._refresh_irfs(select=src.name)

    # ----------------------------------------------------------- instruments
    def _rescan_ports(self) -> None:
        current = self.port_combo.currentData()
        self.port_combo.clear()
        self.port_combo.addItem("scan all ports", None)
        for p in sorted(serial.tools.list_ports.comports(), key=lambda p: p.device):
            self.port_combo.addItem(f"{p.device}: {p.description}", p.device)
        remembered = current or self.settings.get("port")
        if remembered:
            i = self.port_combo.findData(remembered)
            if i >= 0:
                self.port_combo.setCurrentIndex(i)

    def _toggle_ambit(self) -> None:
        if self.ambit is not None:
            try:
                self.ambit.close()
            finally:
                self.ambit = None
            self.lbl_ambit.setText("not connected")
            self.btn_ambit.setText("Connect AMBIT")
            return
        port = self.port_combo.currentData()
        self.btn_ambit.setEnabled(False)
        self._start_task(lambda log: ag.find_ambit(port, log=log), self._ambit_connected,
                         lambda e: (self._log("ERROR: " + e), self.btn_ambit.setEnabled(True)))

    def _ambit_connected(self, link) -> None:
        self.ambit = link
        self.btn_ambit.setEnabled(True)
        self.btn_ambit.setText("Disconnect AMBIT")
        self.lbl_ambit.setText(f"{link.name} on {link.ser.port}, firmware {link.fw or '?'}")
        self.settings.set("port", link.ser.port)

    def _toggle_felix(self) -> None:
        if self.felix is not None:
            try:
                self.felix.close()
            finally:
                self.felix = None
            self.lbl_felix.setText("not connected")
            self.btn_felix.setText("Connect Felix")
            self.lbl_rates.setText("-")
            return
        thr = self.spin_threshold.value()

        def fn(log):
            f = ag.FelixLink(log=log).open()
            f.set_thresholds(thr)
            rates = f.count_rates(0.5)
            rates["tdc_rejection"] = f.tdc_rejection(0.5)
            return f, rates

        self.btn_felix.setEnabled(False)
        self._start_task(fn, self._felix_connected, self._felix_failed)

    def _felix_connected(self, result) -> None:
        self.felix, rates = result
        self.btn_felix.setEnabled(True)
        self.btn_felix.setText("Disconnect Felix")
        self.lbl_felix.setText(self.felix.description)
        self._show_rates(rates)

    def _felix_failed(self, text: str) -> None:
        self._log("ERROR: " + text)
        self.btn_felix.setEnabled(True)
        if ag.felix_usb_state().get("problem") == ag.CM_PROB_FAILED_INSTALL:
            self._offer_driver_install()

    # Felix driver: Windows shows the Felix with problem code 28 when it has no driver yet.
    def _check_driver_at_start(self) -> None:
        if ag.felix_usb_state().get("problem") == ag.CM_PROB_FAILED_INSTALL:
            self._log("the Felix is plugged in but Windows has no driver for it")
            self._offer_driver_install()

    def _offer_driver_install(self) -> None:
        path = drv.installer_path()
        if path is None:
            self._log(f"the driver installer is not bundled with this build: get it from {drv.DOWNLOAD_URL}")
            return
        yes = QtWidgets.QMessageBox.StandardButton.Yes
        answer = QtWidgets.QMessageBox.question(
            self, "Felix driver",
            "The Felix is plugged in, but Windows has no driver for it yet.\n\n"
            "Install TEDIEL's Felix WinUSB driver now? Windows asks for administrator permission, "
            "then a short wizard runs (Next, Finish). The Felix is connected afterwards.",
            yes | QtWidgets.QMessageBox.StandardButton.No, yes)
        if answer != yes:
            self._log("driver not installed; 'Connect Felix' asks again")
            return

        def fn(log):
            log(f"running {path.name} (administrator permission needed) ...")
            code = drv.run_installer(path)
            state = ag.felix_usb_state()
            for _ in range(20):                                  # Windows binds the driver after dpinst exits
                if state.get("present") and not state.get("problem"):
                    break
                time.sleep(0.5)
                state = ag.felix_usb_state()
            return code, state

        self.btn_felix.setEnabled(False)
        self._start_task(fn, self._driver_installed,
                         lambda e: (self._log("ERROR: " + e), self.btn_felix.setEnabled(True)))

    def _driver_installed(self, result) -> None:
        code, state = result
        self.btn_felix.setEnabled(True)
        if state.get("present") and not state.get("problem"):
            self._log(f"Felix driver installed (installer exit code {code}); connecting ...")
            if self.felix is None:
                self._toggle_felix()
        else:
            self._log(f"the installer finished (exit code {code}) but the Felix still reports problem "
                      f"{state.get('problem')}: unplug and replug it, then 'Connect Felix'")

    def _show_rates(self, rates: dict) -> None:
        sync, spad = rates.get(ag.SYNC_CH, 0.0), rates.get(ag.SPAD_CH, 0.0)
        text = f"SYNC {sync / 1e6:.2f} MHz, SPAD {spad / 1e3:,.0f} kcps"
        if sync < 0.9 * ag.LASER_REP_RATE_HZ:
            text += "  (SYNC below the laser rate: laser off or cable?)"
        if spad > bd.SPAD_MAX_CPS:
            text += "  (SPAD saturated)"
        warn = ag.rejection_warning(rates.get("tdc_rejection", float("nan")), spad)
        if warn:
            text += "  (TDC DROPS THE SPAD EVENTS, see the log)"
            self._log("ERROR: " + warn)
        self.lbl_rates.setText(text)

    def _read_rates(self) -> None:
        if self.felix is None:
            self._log("connect the Felix first")
            return
        felix, thr = self.felix, self.spin_threshold.value()

        def fn(log):
            felix.set_thresholds(thr)
            rates = felix.count_rates(0.5)
            rates["tdc_rejection"] = felix.tdc_rejection(0.5)
            return rates

        self._start_task(fn, self._show_rates)

    # --------------------------------------------------------------- openJII
    def _env_changed(self) -> None:
        self.key_edit.setText(self.settings.api_key(self.env_combo.currentData()))
        self.client = None
        self.lbl_user.setText("not signed in")
        self.exp_combo.clear()

    def _sign_in(self) -> None:
        env = self._env()
        key = self.key_edit.text().strip()
        if not key:
            self._log(f"paste an API key: {env.api_keys_url}")
            return
        client = oj.OpenJIIClient(env, key)

        def fn(log):
            user = client.validate_key()
            log(f"openJII {env.key}: signed in as {user.get('name') or user.get('email') or user.get('id')}")
            return client, user, client.list_experiments()

        self.btn_signin.setEnabled(False)
        self._start_task(fn, self._signed_in, lambda e: (self._log("ERROR: " + e), self.btn_signin.setEnabled(True)))

    def _signed_in(self, result) -> None:
        self.client, user, exps = result
        self.btn_signin.setEnabled(True)
        self.lbl_user.setText(f"{user.get('name') or ''} <{user.get('email') or ''}>".strip())
        self.settings.set_api_key(self._env().key, self.client.api_key)
        self.settings.set("env", self._env().key)
        self.settings.save()
        self._fill_experiments(exps)

    def _fill_experiments(self, exps: list[oj.Experiment]) -> None:
        self.experiments = exps
        self.exp_combo.clear()
        for e in exps:
            self.exp_combo.addItem(f"{e.name}  ({e.id[:8]})", e.id)
        i = self.exp_combo.findData(self.settings.get("experiment_id"))
        if i >= 0:
            self.exp_combo.setCurrentIndex(i)
        self._log(f"{len(exps)} active experiment(s) you are a member of")

    def _refresh_experiments(self) -> None:
        if self.client is None:
            self._log("sign in first")
            return
        client = self.client
        self._start_task(lambda log: client.list_experiments(), self._fill_experiments)

    def _upload_current(self) -> None:
        if self.res is None or self.summary is None:
            self._log("nothing to upload: run a measurement or open a saved run first")
            return
        if self.client is None or not self.exp_combo.currentData():
            self._log("upload: sign in to openJII and choose an experiment first")
            return
        table = self.table_edit.text().strip()
        if not oj.TABLE_NAME_RE.match(table):
            self._log("upload: the table name needs letters, digits and underscores, starting with a letter")
            return
        if not np.isfinite(self.summary["phi2_ambit"]) and not np.isfinite(self.summary["phi2_tau_mean"]):
            self._log("upload: this run has no phi2 at all (see the warnings); uploading its raw data anyway")
        client, exp_id = self.client, self.exp_combo.currentData()
        res, summ, run_dir = self.res, self.summary, self.run_dir
        include = self.chk_felix.isChecked()
        label = self.label_edit.text().strip()
        self.settings.set("experiment_id", exp_id)
        self.settings.set("table_name", table)

        def fn(log):
            path, n = up.write_ndjson(run_dir, up.run_rows(res, summ, include_felix=include, label=label))
            log(f"upload: {n} rows written to {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
            return up.upload_and_wait(client, exp_id, table, path, log=log)

        self.btn_upload.setEnabled(False)
        self.lbl_upload.setText("uploading ...")
        self._start_task(fn, self._upload_done, self._upload_failed)

    def _upload_done(self, rec: oj.UploadRecord) -> None:
        self.btn_upload.setEnabled(True)
        text = up.describe(rec)
        self.lbl_upload.setText(text)
        self._log("upload " + text)
        if self.run_dir is not None:
            try:
                (self.run_dir / "openjii_upload.json").write_text(json.dumps(rec.__dict__, indent=1), encoding="utf-8")
            except OSError:
                pass

    def _upload_failed(self, text: str) -> None:
        self.btn_upload.setEnabled(True)
        if "No upload table with id" in text:
            text = ("the platform does not list the new table yet (it appears a few minutes after its "
                    "first upload finished): try again in a minute")
        self.lbl_upload.setText("failed: " + text.splitlines()[0])
        self._log("ERROR upload: " + text)

    # ----------------------------------------------------------- measurement
    def _pauses(self) -> tuple[int, ...]:
        text = self.pauses_edit.text().replace(";", ",").replace(" ", ",")
        vals = tuple(int(v) for v in text.split(",") if v.strip())
        if not vals:
            raise ValueError("at least one pause")
        return vals

    def _protocol(self) -> list[dict]:
        return ap.train_protocol(self.spin_level.value(), self.spin_dark.value(), self.spin_block.value(),
                                 self._pauses(), self.spin_freq.value())

    def _update_protocol_label(self, *_args) -> None:
        try:
            prot = self._protocol()
            pred = ap.predict(prot)
            n_pause = sum(1 for r in pred if r["role"] == "pause")
            dark = [r["dark1"] - r["dark0"] for r in pred if r["role"] == "pause"]
            self.lbl_protocol.setText(f"{len(prot)} AMBIT lines, {sum(s['pulses'] for s in prot)} points, "
                                      f"{pred[-1]['t1']:.2f} s; {n_pause} pauses of "
                                      f"{', '.join(f'{d * 1e3:.0f}' for d in dark)} ms dark (predicted)")
            self.lbl_protocol.setStyleSheet("color: gray")
        except ValueError as exc:
            self.lbl_protocol.setText(str(exc))
            self.lbl_protocol.setStyleSheet(f"color: {RED}")

    def _run_measurement(self) -> None:
        if self.ambit is None or self.felix is None:
            self._log("connect the AMBIT and the Felix first")
            return
        try:
            protocol = self._protocol()
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "Protocol", str(exc))
            return
        warn = ap.window_check(protocol, self.spin_window.value() / 1e3)
        if warn:                                             # the pauses would yield no decay
            QtWidgets.QMessageBox.warning(self, "Felix window vs pauses", warn)
            return
        label = self.label_edit.text().strip() or "run"
        out = self.session / f"pauses_{label}_{time.strftime('%H%M%S')}"
        params = dict(window_s=self.spin_window.value() / 1e3, note=self.note_edit.text().strip(),
                      threshold_v=self.spin_threshold.value(), laser_hz=ag.LASER_REP_RATE_HZ)
        fo, fm = self.spin_fo.value(), self.spin_fm.value()
        ambit, felix, irf = self.ambit, self.felix, self.irf
        self._save_settings()

        def fn(log):
            d = ap.run_pauses(ambit, felix, out, protocol, log=log, **params)
            log("analysing (reconvolution fits) ...")
            try:
                res = ap.analyse_pauses(d, irf=irf)
                return d, res, ap.run_summary(res, fo, fm)
            except Exception as exc:                         # the raw data are shown anyway
                log(f"ERROR: analysis failed: {exc}")
                return d, None, None

        self._busy(True)
        self._log(f"run: {out.name}")
        self._start_task(fn, self._run_done)

    def _run_done(self, result) -> None:
        self._busy(False)
        self._show_run(*result)
        if self.chk_upload.isChecked() and self.summary is not None:
            self._upload_current()

    def _open_run_dialog(self) -> None:
        start = str(self.run_dir.parent if self.run_dir else RUNS_DIR)
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Run folder (with felix.npz)", start)
        if path:
            self.open_run(Path(path))

    def open_run(self, path: Path, threaded: bool = True) -> None:
        """Analyse a saved run folder and show it (no instruments needed)."""
        path = Path(path)
        if not (path / "felix.npz").exists():
            self._log(f"{path}: no felix.npz, not a run folder")
            return
        fo, fm, irf = self.spin_fo.value(), self.spin_fm.value(), self.irf

        def fn(log):
            log(f"analysing {path.name} ...")
            try:
                res = ap.analyse_pauses(path, irf=irf)
                return path, res, ap.run_summary(res, fo, fm)
            except Exception as exc:                         # the raw data are shown anyway
                log(f"ERROR: analysis failed: {exc}")
                return path, None, None

        if threaded:
            self._busy(True)
            self._start_task(fn, self._run_done)
        else:
            self._show_run(*fn(self._log))

    def _show_run(self, run_dir: Path, res: dict | None, summ: dict | None) -> None:
        self.run_dir, self.res, self.summary = Path(run_dir), res, summ
        self._set_figures(res)                               # the raw-data tab first, always
        if summ is None:
            self.lbl_title.setText(f"<b>{self.run_dir.name}</b> &nbsp; <span style='color: {RED}'>the analysis "
                                   "failed (see the log): the raw Felix and AMBIT data are in the 'Raw data' tab</span>")
            self.lbl_summary.setText("")
            self.lbl_warnings.setText("")
            self.table.clear()
            self.table.setRowCount(0)
            self.lbl_table_note.setText("")
            self.tabs.setCurrentIndex(1 if self.tabs.count() > 1 else 0)
            self.lbl_upload.setText("-")
            return
        try:
            (self.run_dir / "results.json").write_text(json.dumps(summ, indent=1, default=float), encoding="utf-8")
        except OSError as exc:
            self._log(f"could not save results.json: {exc}")
        self._fill_results(summ)
        self.tabs.setCurrentIndex(0)
        self._log(f"{summ['run']}: phi2 AMBIT {_pm(summ['phi2_ambit'], summ['phi2_ambit_err'])}, "
                  f"phi2 lifetime {_pm(summ['phi2_tau_mean'], summ['phi2_tau_sd'])}"
                  + "".join(f"\n  WARNING: {w}" for w in summ["warnings"]))
        self.lbl_upload.setText("-")

    # --------------------------------------------------------------- results
    def _fill_results(self, s: dict) -> None:
        self.lbl_title.setText(f"<b>{s['run']}</b>" + (f" &nbsp; {s['note']}" if s.get("note") else ""))
        n = len(s["phi2_tau"])
        big = "font-size: 14pt; font-weight: bold"
        html = (f"<p>phi2 from the AMBIT intensity, (Fm - Fo) / Fm = <span style='{big}'>"
                f"{_pm(s['phi2_ambit'], s['phi2_ambit_err'])}</span><br>"
                f"<span style='color: gray'>Fo = {_pm(s['fo_ambit'], s['fo_ambit_se'], 5)} ({s['n_fo']} points), "
                f"Fm = {_pm(s['fm_ambit'], s['fm_ambit_se'], 5)} ({s['n_fm']} points); raw s_630 median "
                f"{_num(s['s_630_fo_median'], 0)} / {_num(s['s_630_fm_median'], 0)} counts</span></p>"
                f"<p>phi2 from the lifetimes, (t_amp,Fm - t_amp,Fo) / t_amp,Fm = <span style='{big}'>"
                f"{_pm(s['phi2_tau_mean'], s['phi2_tau_sd'])}</span> &nbsp;"
                f"<span style='color: gray'>mean ± SD of {n} pauses: "
                f"{', '.join(_num(p) for p in s['phi2_tau'])}</span><br>"
                f"<b>t_amp</b> Fo = {_pm(s['t_amp_fo'], s['t_amp_fo_err'])} ns, "
                f"<span style='background: {HIGHLIGHT.name()}'><b>mean t_amp of the {n} Fm pauses = "
                f"{_pm(s['t_amp_fm_mean'], s['t_amp_fm_sd'])} ns</b></span><br>"
                f"<span style='color: gray'>t_int Fo = {_pm(s['t_int_fo'], s['t_int_fo_err'])} ns, "
                f"mean t_int of the pauses = {_pm(s['t_int_fm_mean'], s.get('t_int_fm_sd'))} ns; "
                f"tau2 Fo {_pm(s['tau2_fo'], s['tau2_fo_err'])} ns, mean tau2 of the pauses "
                f"{_pm(s['tau2_fm_mean'], s['tau2_fm_sd'])} ns</span></p>")
        self.lbl_summary.setText(html)
        if s["warnings"]:
            self.lbl_warnings.setText(f"<ul style='color: {RED}; margin-top: 0'>"
                                      + "".join(f"<li>{w}</li>" for w in s["warnings"]) + "</ul>")
        else:
            self.lbl_warnings.setText("<span style='color: gray'>no warnings</span>")
        self._fill_table(s)

    def _fill_table(self, s: dict) -> None:
        cols = ["decay", "laser counts", "seconds", "tau1 [ns]", "share 1", "tau2 [ns]", "share 2", "t_amp [ns]",
                "t_int [ns]", "t_mean [ns]", "IRF shift [ps]", "chi2r", "phi2 lifetime"]
        fits = [("reference (Fo): " + s["ref_name"], s["fits"]["reference"], None, None)]
        for k, f in enumerate(s["fits"]["pauses"]):
            fits.append((f"Fm pause {k + 1} ({s['dark_ms'][k]:.0f} ms dark)", f, s["phi2_tau"][k], s["phi2_tau_err"][k]))
        n = len(s["fits"]["pauses"])
        self.table.clear()
        self.table.setColumnCount(len(cols))
        self.table.setHorizontalHeaderLabels(cols)
        self.table.setRowCount(len(fits) + 1)

        def put(r, c, text, bold=False, highlight=False, color=None):
            it = QtWidgets.QTableWidgetItem(text)
            it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter if c else QtCore.Qt.AlignmentFlag.AlignVCenter)
            if bold:
                font = it.font()
                font.setBold(True)
                it.setFont(font)
            if highlight:
                it.setBackground(HIGHLIGHT)
            if color:
                it.setForeground(QtGui.QColor(color))
            self.table.setItem(r, c, it)

        for r, (name, f, phi, perr) in enumerate(fits):
            put(r, 0, name)
            if f is None:
                put(r, 1, "no fit", color=RED)
                for c in range(2, len(cols)):
                    put(r, c, "-")
                continue
            put(r, 1, f"{f['counts']:,.0f}")
            put(r, 2, f"{f['seconds']:.3f}")
            for c, i in ((3, 1), (5, 2)):
                if f.get(f"empty{i}"):
                    put(r, c, "none", color="gray")
                    put(r, c + 1, "0 %", color="gray")
                else:
                    put(r, c, _pm(f[f"tau{i}"], f[f"tau{i}_err"]) + ("" if f[f"resolved{i}"] else " *"))
                    put(r, c + 1, f"{100 * f[f'share{i}']:.0f} %")
            put(r, 7, _pm(f["t_amp"], f["t_amp_err"]))
            put(r, 8, _pm(f["t_int"], f["t_int_err"]))
            put(r, 9, _pm(f["t_mean"], f["t_mean_err"]))
            put(r, 10, f"{f['shift_ps']:+.0f}")
            put(r, 11, f"{f['chi2r']:.2f}" + ("  (bound)" if f["at_bound"] else ""))
            put(r, 12, _pm(phi, perr) if phi is not None else "")
        r = len(fits)
        put(r, 0, f"mean of the {n} Fm pauses", bold=True, highlight=True)
        for c in range(1, len(cols)):
            put(r, c, "", highlight=True)
        put(r, 5, _pm(s["tau2_fm_mean"], s["tau2_fm_sd"]) + "  (± SD)", highlight=True)
        put(r, 7, _pm(s["t_amp_fm_mean"], s["t_amp_fm_sd"]) + "  (± SD)", bold=True, highlight=True)
        put(r, 8, _pm(s["t_int_fm_mean"], s.get("t_int_fm_sd")) + "  (± SD)", highlight=True)
        put(r, 12, _pm(s["phi2_tau_mean"], s["phi2_tau_sd"]) + "  (± SD)", bold=True, highlight=True)
        self.table.resizeColumnsToContents()
        self.table.horizontalHeader().setStretchLastSection(True)
        self.lbl_table_note.setText("2-exponential IRF reconvolution fit of every decay (reconv.py): tau1/tau2 with "
                                    "their fit errors (*: shorter than half the IRF FWHM, unresolved; none: no photons, "
                                    "the decay is a single exponential), share = photons of the component. t_amp = "
                                    "amplitude-weighted lifetime, proportional to the yield (the amplitude of an "
                                    "unresolved component is poorly determined and carries most of its error); t_int = "
                                    "photon-weighted mean lifetime; t_mean = model-free mean arrival time. "
                                    "phi2 lifetime = (t_amp,Fm - t_amp,Fo) / t_amp,Fm per pause.")

    def _set_figures(self, res: dict | None) -> None:
        while self.tabs.count() > 1:
            w = self.tabs.widget(1)
            self.tabs.removeTab(1)
            w.deleteLater()
        try:
            self._add_figure_tab("Raw data", ap.plot_raw(self.run_dir, irf=self.irf, save=True, show=False))
        except Exception as exc:
            self._log(f"ERROR: raw-data plot failed: {exc}")
        if res is None:
            return
        figs = ap.plot_pauses(res, save=True, show=False, ref="reference")
        for name, fig in zip(("Overview", "Pause zoom", "Fm vs Fo decays"), figs):
            if fig is not None:
                self._add_figure_tab(name, fig)
        keys = [k for k, _ in ap.fit_groups(res)]
        for key, fig in zip(keys, ap.plot_fits(res, save=True, show=False)):
            if fig is not None:
                self._add_figure_tab("Decay fits: reference" if key == "laser" else f"Decay fits: Fm pauses ({key})", fig)

    def _add_figure_tab(self, title: str, fig) -> None:
        canvas = FigureCanvasQTAgg(fig)
        w, h = (fig.get_size_inches() * fig.dpi).astype(int)
        canvas.setMinimumSize(int(w), int(h))
        inner = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(inner)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(NavigationToolbar2QT(canvas, inner))
        v.addWidget(canvas, 1)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(inner)
        scroll.setWidgetResizable(True)
        self.tabs.addTab(scroll, title)
        plt.close(fig)                                   # pyplot's bookkeeping; the figure lives in the canvas

    # -------------------------------------------------------------- settings
    def _load_settings(self) -> None:
        s = self.settings
        i = self.env_combo.findData(s.get("env", "prod"))
        self.env_combo.setCurrentIndex(max(i, 0))
        self.key_edit.setText(s.api_key(self.env_combo.currentData()))
        self.table_edit.setText(s.get("table_name", DEFAULT_TABLE))
        self.chk_upload.setChecked(bool(s.get("upload_after_run", False)))
        self.chk_felix.setChecked(bool(s.get("include_felix", True)))
        self.spin_threshold.setValue(float(s.get("threshold_v", ag.THRESHOLD_V)))
        self.spin_level.setValue(int(s.get("level", 255)))
        self.spin_dark.setValue(int(s.get("dark_points", 100)))
        self.spin_block.setValue(float(s.get("block_s", 0.2)))
        self.pauses_edit.setText(str(s.get("pauses", "2,2,2,2")))
        self.spin_freq.setValue(int(s.get("freq", 100)))
        self.spin_window.setValue(float(s.get("window_ms", 5.0)))
        self.spin_fo.setValue(int(s.get("fo_points", ap.FO_POINTS)))
        self.spin_fm.setValue(int(s.get("fm_points", ap.FM_POINTS)))
        self.label_edit.setText(str(s.get("label", "leaf")))
        self.note_edit.setText(str(s.get("note", "")))
        geo = s.get("geometry")
        if geo:
            self.restoreGeometry(QtCore.QByteArray.fromBase64(geo.encode()))
        else:
            self.resize(1500, 900)

    def _save_settings(self) -> None:
        s = self.settings
        s.set("env", self.env_combo.currentData())
        s.set_api_key(self.env_combo.currentData(), self.key_edit.text().strip())
        if self.exp_combo.currentData():
            s.set("experiment_id", self.exp_combo.currentData())
        s.set("table_name", self.table_edit.text().strip())
        s.set("upload_after_run", self.chk_upload.isChecked())
        s.set("include_felix", self.chk_felix.isChecked())
        s.set("threshold_v", self.spin_threshold.value())
        s.set("level", self.spin_level.value())
        s.set("dark_points", self.spin_dark.value())
        s.set("block_s", self.spin_block.value())
        s.set("pauses", self.pauses_edit.text())
        s.set("freq", self.spin_freq.value())
        s.set("window_ms", self.spin_window.value())
        s.set("fo_points", self.spin_fo.value())
        s.set("fm_points", self.spin_fm.value())
        s.set("label", self.label_edit.text())
        s.set("note", self.note_edit.text())
        s.set("geometry", bytes(self.saveGeometry().toBase64()).decode())
        try:
            s.save()
        except OSError as exc:
            self._log(f"could not save the settings: {exc}")

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self._save_settings()
        for name in ("ambit", "felix"):
            dev = getattr(self, name)
            if dev is not None:
                try:
                    dev.close()
                except Exception:
                    pass
        event.accept()
