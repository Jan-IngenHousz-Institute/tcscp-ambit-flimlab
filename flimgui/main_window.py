"""Main window: wires device, settings, IRF store, fit worker and widgets together."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field

import numpy as np
from PySide6.QtCore import QObject, QThread, QTimer, Signal, Slot
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QMainWindow,
    QMessageBox,
    QScrollArea,
    QStackedWidget,
    QStatusBar,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .device import SYNC_CH, DeviceBase, HistStatus
from .fitting import FitResult, fit_decay
from .irf import IRF, IRFStore
from .settings import Settings
from .widgets.channels import ChannelPanel
from .widgets.counters import CounterPanel
from .widgets.fit_panel import FitPanel
from .widgets.histogram import HistogramConfigPanel, HistogramView
from .widgets.irf_panel import IRFPanel

log = logging.getLogger(__name__)

COUNTER_INTEGRATION_S = 0.5
FIT_MIN_INTERVAL_S = 0.4


class FitWorker(QObject):
    done = Signal(int, object)

    def __init__(self):
        super().__init__()
        self._p0 = {}

    @Slot(int, object, float, int, object, object)
    def fit(self, idx, counts, bin_ns, n_exp, irf, period_ns):
        key = (idx, n_exp, irf is not None)
        try:
            r = fit_decay(counts, bin_ns, n_exp, irf, period_ns, p0=self._p0.get(key))
            if r.success:
                self._p0[key] = r.params
        except Exception as exc:  # never let the worker die
            log.exception("fit failed")
            r = fit_decay(np.zeros(8), bin_ns, n_exp)  # yields an empty result
            r.message = f"fit error: {exc}"
        self.done.emit(idx, r)


@dataclass
class HistTab:
    idx: int
    view: HistogramView
    config: HistogramConfigPanel
    left_ns: float = 0.0
    bin_ns: float = 0.0
    n: int = 0
    t_ns: np.ndarray | None = None
    bins: np.ndarray | None = None
    status: HistStatus = field(default_factory=HistStatus)
    started: float = 0.0
    fit_busy: bool = False
    last_fit_time: float = 0.0
    last_fit: FitResult | None = None
    irf_axis: np.ndarray | None = None  # active IRF on this tab's axis
    irf_resampled: bool = False


class MainWindow(QMainWindow):
    fit_requested = Signal(int, object, float, int, object, object)

    def __init__(self, device: DeviceBase, settings: Settings, irf_store: IRFStore):
        super().__init__()
        self.device = device
        self.settings = settings
        self.irf_store = irf_store
        self.active_irf: IRF | None = None
        self.tabs: list[HistTab] = []
        self._rec_timer: QTimer | None = None
        self._rec_tab: HistTab | None = None
        self._rec_t0 = 0.0
        self._rec_duration = 0.0
        self._rec_name = ""
        self._closing = False

        self.setWindowTitle(f"FLIM GUI - {device.name}")
        self.names = device.channel_names()
        self._build_ui()
        self._start_fit_thread()
        self._connect_device()
        self._restore_settings()

    # ------------------------------------------------------------------ UI
    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        # left column: counters + channels
        left = QVBoxLayout()
        self.counters = CounterPanel(self.names, SYNC_CH)
        left.addWidget(self.counters)
        thr_ranges = {ch: self.device.threshold_range(ch) for ch in range(len(self.names))}
        dly_ranges = {ch: tuple(x * 1e9 for x in self.device.delay_range(ch)) for ch in range(len(self.names))}
        self.channels = ChannelPanel(self.names, thr_ranges, dly_ranges)
        left.addWidget(self.channels)
        if self.device.is_simulator:
            left.addWidget(self._build_sim_controls())
        left.addStretch(1)
        left_w = QWidget()
        left_w.setLayout(left)
        left_w.setFixedWidth(360)
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setWidget(left_w)
        left_scroll.setFixedWidth(380)
        root.addWidget(left_scroll)

        # centre: histogram tabs
        self.tabw = QTabWidget()
        root.addWidget(self.tabw, 1)

        # right column: per-tab config (stacked) + IRF + fit
        right = QVBoxLayout()
        self.config_stack = QStackedWidget()
        right.addWidget(self.config_stack)
        self.irf_panel = IRFPanel()
        right.addWidget(self.irf_panel)
        self.fit_panel = FitPanel()
        right.addWidget(self.fit_panel)
        right.addStretch(1)
        right_w = QWidget()
        right_w.setLayout(right)
        right_w.setFixedWidth(380)
        right_scroll = QScrollArea()
        right_scroll.setWidgetResizable(True)
        right_scroll.setWidget(right_w)
        right_scroll.setFixedWidth(400)
        root.addWidget(right_scroll)

        for idx in range(max(self.device.n_histograms(), 1)):
            view = HistogramView()
            cfg = HistogramConfigPanel(self.names, f"Histogram {idx + 1}")
            tab = HistTab(idx=idx, view=view, config=cfg)
            self.tabs.append(tab)
            self.tabw.addTab(view, f"Histogram {idx + 1}")
            self.config_stack.addWidget(cfg)
            cfg.apply_requested.connect(lambda _i=idx: self._apply_config(_i))
            cfg.enable_toggled.connect(lambda on, _i=idx: self._toggle_hist(_i, on))
            cfg.reset_requested.connect(lambda _i=idx: self._reset_hist(_i))
            cfg.log_toggled.connect(lambda on, _i=idx: self._set_log(_i, on))
        self.tabw.currentChanged.connect(self._on_tab_changed)

        self.setStatusBar(QStatusBar())
        self.statusBar().showMessage("ready")

        self.channels.threshold_changed.connect(self._on_threshold)
        self.channels.delay_changed.connect(self._on_delay)
        self.channels.auto_requested.connect(self._on_auto_threshold)
        self.irf_panel.record_requested.connect(self._start_irf_record)
        self.irf_panel.cancel_requested.connect(self._cancel_irf_record)
        self.irf_panel.active_changed.connect(self._set_active_irf)
        self.irf_panel.delete_requested.connect(self._delete_irf)
        self.fit_panel.n_exp_changed.connect(self._on_fit_params)
        self.fit_panel.enabled_changed.connect(self._on_fit_params)

    def _build_sim_controls(self) -> QGroupBox:
        box = QGroupBox("Simulator")
        form = QFormLayout(box)
        self.sim_sample = QComboBox()
        self.sim_sample.addItems(["fluorescence", "scatter"])
        self.sim_sample.currentTextChanged.connect(lambda s: setattr(self.device, "sample", s))
        form.addRow("Sample", self.sim_sample)
        rate = QDoubleSpinBox()
        rate.setRange(0, 1e7)
        rate.setDecimals(0)
        rate.setSingleStep(50000)
        rate.setValue(self.device.signal_rate.get(1, 0))
        rate.setSuffix(" cps")
        rate.valueChanged.connect(lambda v: self.device.signal_rate.__setitem__(1, float(v)))
        form.addRow("CH1 rate", rate)
        laser = QDoubleSpinBox()
        laser.setRange(1, 100)
        laser.setValue(self.device.laser_hz / 1e6)
        laser.setSuffix(" MHz")
        laser.valueChanged.connect(lambda v: setattr(self.device, "laser_hz", float(v) * 1e6))
        form.addRow("Laser", laser)
        return box

    # ------------------------------------------------------------- threads
    def _start_fit_thread(self) -> None:
        self._fit_thread = QThread(self)
        self._fit_worker = FitWorker()
        self._fit_worker.moveToThread(self._fit_thread)
        self.fit_requested.connect(self._fit_worker.fit)
        self._fit_worker.done.connect(self._on_fit_done)
        self._fit_thread.start()

    # ------------------------------------------------------------- device
    def _connect_device(self) -> None:
        self.device.counts_updated.connect(self._on_counts)
        self.device.histogram_updated.connect(self._on_histogram)
        self.device.connection_lost.connect(self._on_connection_lost)
        self.device.message.connect(lambda m: self.statusBar().showMessage(m, 15000))

    def _restore_settings(self) -> None:
        s = self.settings
        # thresholds and delays: saved value if present, otherwise what the device has
        for ch in range(len(self.names)):
            row = self.channels.rows[ch]
            key = str(ch)
            try:
                if key in s.thresholds_v:
                    v = self.device.set_threshold(ch, s.thresholds_v[key])
                else:
                    v = self.device.get_threshold(ch)
                row.set_threshold(v)
            except Exception as exc:
                self._warn(f"threshold restore failed on {self.names[ch]}: {exc}")
            try:
                if key in s.delays_ns:
                    d = self.device.set_delay(ch, s.delays_ns[key] * 1e-9)
                else:
                    d = self.device.get_delay(ch)
                row.set_delay_ns(d * 1e9)
            except Exception as exc:
                self._warn(f"delay restore failed on {self.names[ch]}: {exc}")
        try:
            self.device.start_counters(COUNTER_INTEGRATION_S)
        except Exception as exc:
            self._warn(f"could not start counters: {exc}")

        self.fit_panel.n_spin.setValue(max(1, min(4, s.n_exp)))
        self.fit_panel.enable_chk.setChecked(s.fit_enabled)
        self.irf_panel.dur_spin.setValue(s.irf_record_s)
        self._refresh_irf_list()
        if s.active_irf:
            self._set_active_irf(s.active_irf)

        for tab in self.tabs:
            h = s.hist(tab.idx)
            tab.config.set_settings(h)
            tab.config.log_chk.setChecked(s.log_y)
            tab.view.set_log(s.log_y)
            self._apply_config(tab.idx, quiet=True)
            if h.enabled:
                QTimer.singleShot(300, lambda _i=tab.idx: self.tabs[_i].config.enable_btn.setChecked(True))
        if s.window_geometry:
            try:
                from PySide6.QtCore import QByteArray

                self.restoreGeometry(QByteArray.fromHex(s.window_geometry.encode()))
            except Exception:
                pass

    # ------------------------------------------------------------ counters
    @Slot(float, object)
    def _on_counts(self, t_int: float, rates: dict) -> None:
        self.counters.update_rates(rates)

    def _laser_period_ns(self) -> float | None:
        sync = self.counters.sync_rate()
        return 1e9 / sync if sync > 1e5 else None

    # ------------------------------------------------------------ channels
    def _targets(self, ch: int) -> list[int]:
        if ch != SYNC_CH and self.channels.coupled:
            return [c for c in range(len(self.names)) if c != SYNC_CH]
        return [ch]

    def _on_threshold(self, ch: int, volts: float) -> None:
        for c in self._targets(ch):
            try:
                v = self.device.set_threshold(c, volts)
                self.channels.rows[c].set_threshold(v)
                self.settings.thresholds_v[str(c)] = v
            except Exception as exc:
                self._warn(f"threshold on {self.names[c]}: {exc}")
        self._save_settings()

    def _on_delay(self, ch: int, ns: float) -> None:
        for c in self._targets(ch):
            try:
                d = self.device.set_delay(c, ns * 1e-9)
                self.channels.rows[c].set_delay_ns(d * 1e9)
                self.settings.delays_ns[str(c)] = d * 1e9
            except Exception as exc:
                self._warn(f"delay on {self.names[c]}: {exc}")
        self._save_settings()

    def _on_auto_threshold(self, ch: int) -> None:
        self.statusBar().showMessage(f"auto threshold on {self.names[ch]} ...")
        try:
            res = self.device.auto_threshold([ch])
            for c, v in res.items():
                self.channels.rows[c].set_threshold(v)
                self.settings.thresholds_v[str(c)] = v
            self.statusBar().showMessage(f"auto threshold on {self.names[ch]}: {res.get(ch, float('nan')):.3f} V", 8000)
            self._save_settings()
        except Exception as exc:
            self._warn(f"auto threshold failed: {exc}")

    # ---------------------------------------------------------- histograms
    def _apply_config(self, idx: int, quiet: bool = False) -> None:
        tab = self.tabs[idx]
        h = tab.config.get_settings()
        if h.right_ns <= h.left_ns:
            self._warn("Right bound must be larger than the left bound.")
            return
        if h.meas_ch == h.ref_ch:
            self._warn("Reference and measure channel must differ.")
            return
        try:
            self.device.hist_configure(
                idx,
                h.ref_ch,
                h.meas_ch,
                h.left_ns * 1e-9,
                h.right_ns * 1e-9,
                h.bin_ps * 1e-12 if h.bin_ps > 0 else None,
                h.multi_hit,
                h.accumulate,
            )
            left, right, bw = self.device.hist_bounds(idx)
        except Exception as exc:
            log.warning("histogram %d configuration failed: %s", idx + 1, exc)
            if not quiet:
                self._warn(f"histogram {idx + 1} configuration failed: {exc}")
            return
        tab.left_ns = left * 1e9
        tab.bin_ns = bw * 1e9
        tab.n = int(round((right - left) / bw)) if bw > 0 else 0
        tab.t_ns = tab.left_ns + (np.arange(tab.n) + 0.5) * tab.bin_ns
        tab.bins = None
        tab.last_fit = None
        tab.config.show_bounds(tab.left_ns, right * 1e9, tab.bin_ns * 1000, tab.n)
        self._update_tab_irf(tab)
        self.settings.set_hist(idx, h)
        self._save_settings()

    def _toggle_hist(self, idx: int, on: bool) -> None:
        tab = self.tabs[idx]
        h = tab.config.get_settings()
        if on:
            # safety: no reference signal, no acquisition
            ref_rate = self.counters.rates.get(h.ref_ch, 0.0)
            if ref_rate <= 0:
                QMessageBox.warning(
                    self,
                    "No reference signal",
                    f"{self.names[h.ref_ch]} shows no counts. Check the SYNC cable, "
                    "the laser mode selector and the threshold before enabling.",
                )
                tab.config.set_enabled_state(False)
                return
            self._apply_config(idx, quiet=True)
            try:
                actual = self.device.hist_enable(idx, h.refresh_s)
            except Exception as exc:
                self._warn(f"could not enable histogram {idx + 1}: {exc}")
                tab.config.set_enabled_state(False)
                return
            tab.started = time.monotonic()
            tab.config.set_running(True)
            tab.config.show_status(f"running, refresh {actual:.2f} s")
        else:
            try:
                self.device.hist_disable(idx)
            except Exception as exc:
                self._warn(f"could not disable histogram {idx + 1}: {exc}")
            tab.config.set_running(False)
            tab.config.show_status("stopped")
        h.enabled = on
        self.settings.set_hist(idx, h)
        self._save_settings()

    def _reset_hist(self, idx: int) -> None:
        tab = self.tabs[idx]
        try:
            self.device.hist_reset(idx)
        except Exception as exc:
            self._warn(f"reset failed: {exc}")
            return
        tab.started = time.monotonic()
        tab.last_fit = None
        tab.view.clear_fit()
        if self.tabw.currentIndex() == idx:
            self.fit_panel.show_result(None)

    def _set_log(self, idx: int, on: bool) -> None:
        self.tabs[idx].view.set_log(on)
        self.settings.log_y = on
        self._save_settings()

    @Slot(int, object, object)
    def _on_histogram(self, idx: int, bins, status: HistStatus) -> None:
        if idx >= len(self.tabs) or self._closing:
            return
        tab = self.tabs[idx]
        tab.status = status
        if not status.running or bins is None:
            tab.config.set_enabled_state(False)
            tab.config.show_status("stopped")
            return
        bins = np.asarray(bins, float)
        if tab.t_ns is not None and bins.size > tab.n > 0:
            # the Felix always sends its full 8192-bin accumulator; only the first
            # (right - left) / binWidth entries belong to the configured window
            bins = bins[: tab.n]
        elif tab.t_ns is None or tab.t_ns.size != bins.size:
            tab.n = bins.size
            tab.t_ns = tab.left_ns + (np.arange(tab.n) + 0.5) * tab.bin_ns
            self._update_tab_irf(tab)
        tab.bins = bins
        tab.view.set_data(tab.t_ns, bins)
        elapsed = time.monotonic() - tab.started if tab.started else 0.0
        flags = []
        if status.lost:
            flags.append("windows lost")
        if status.overflow:
            flags.append("bin overflow")
        tab.config.show_status(
            f"running | {int(bins.sum()):,} events | window {status.win_id} | {elapsed:.0f} s"
            + (" | " + ", ".join(flags) if flags else ""),
            warning=bool(flags),
        )
        self._maybe_fit(tab)

    # ----------------------------------------------------------------- fit
    def _maybe_fit(self, tab: HistTab) -> None:
        if not self.fit_panel.enabled or tab.fit_busy or tab.bins is None:
            return
        now = time.monotonic()
        if now - tab.last_fit_time < FIT_MIN_INTERVAL_S:
            return
        if tab.bins.sum() < 50:
            return
        tab.fit_busy = True
        tab.last_fit_time = now
        self.fit_requested.emit(
            tab.idx,
            tab.bins.copy(),
            float(tab.bin_ns),
            int(self.fit_panel.n_exp),
            None if tab.irf_axis is None else tab.irf_axis.copy(),
            self._laser_period_ns(),
        )

    @Slot(int, object)
    def _on_fit_done(self, idx: int, r: FitResult) -> None:
        if idx >= len(self.tabs):
            return
        tab = self.tabs[idx]
        tab.fit_busy = False
        tab.last_fit = r
        if tab.bins is not None and r.model.size == tab.bins.size:
            if r.mode != "none":
                tab.view.set_fit(r.model, tab.left_ns + r.resid_t_ns, r.resid)
            else:
                tab.view.clear_fit()
        if self.tabw.currentIndex() == idx:
            self.fit_panel.show_result(r)

    def _on_fit_params(self, *_) -> None:
        self.settings.n_exp = self.fit_panel.n_exp
        self.settings.fit_enabled = self.fit_panel.enabled
        for tab in self.tabs:
            tab.last_fit_time = 0.0
            if not self.fit_panel.enabled:
                tab.view.clear_fit()
        if not self.fit_panel.enabled:
            self.fit_panel.show_result(None)
        self._save_settings()

    def _on_tab_changed(self, idx: int) -> None:
        self.config_stack.setCurrentIndex(idx)
        if 0 <= idx < len(self.tabs):
            self.fit_panel.show_result(self.tabs[idx].last_fit)
            self._show_active_irf_info(self.tabs[idx])

    # ----------------------------------------------------------------- IRF
    def _refresh_irf_list(self) -> None:
        self.irf_panel.set_names(self.irf_store.names(), self.active_irf.name if self.active_irf else "")

    def _set_active_irf(self, name: str) -> None:
        if not name:
            self.active_irf = None
        else:
            try:
                self.active_irf = self.irf_store.load(name)
            except Exception as exc:
                self._warn(f"could not load IRF '{name}': {exc}")
                self.active_irf = None
        self.settings.active_irf = self.active_irf.name if self.active_irf else ""
        for tab in self.tabs:
            self._update_tab_irf(tab)
            tab.last_fit_time = 0.0
        self._refresh_irf_list()
        self._show_active_irf_info(self.tabs[self.tabw.currentIndex()])
        self._save_settings()

    def _update_tab_irf(self, tab: HistTab) -> None:
        if self.active_irf is None or tab.t_ns is None or tab.n == 0:
            tab.irf_axis = None
            tab.irf_resampled = False
        else:
            tab.irf_axis, exact = self.active_irf.on_axis(tab.left_ns, tab.bin_ns, tab.n)
            tab.irf_resampled = not exact
            if tab.irf_axis.sum() <= 0:
                tab.irf_axis = None
        tab.view.set_irf(tab.irf_axis)

    def _show_active_irf_info(self, tab: HistTab) -> None:
        note = ""
        if self.active_irf is not None and tab.irf_axis is None and tab.n > 0:
            note = "<span style='color:#b00000'>IRF does not overlap this histogram range; tail fit used</span>"
        self.irf_panel.set_active_info(self.active_irf, tab.irf_resampled, note)

    def _delete_irf(self, name: str) -> None:
        if QMessageBox.question(self, "Delete IRF", f"Delete stored IRF '{name}'?") != QMessageBox.Yes:
            return
        self.irf_store.delete(name)
        if self.active_irf and self.active_irf.name == name:
            self._set_active_irf("")
        self._refresh_irf_list()

    def _start_irf_record(self, duration: float, name: str) -> None:
        idx = self.tabw.currentIndex()
        tab = self.tabs[idx]
        if not tab.config.enable_btn.isChecked():
            QMessageBox.information(
                self,
                "Histogram not running",
                "Enable the histogram on this tab first, with the scatterer in place and the emission filter removed.",
            )
            return
        if not tab.config.acc_chk.isChecked():
            QMessageBox.information(
                self,
                "Accumulate is off",
                "Switch Accumulate on for this histogram, so the IRF integrates over the whole recording time.",
            )
            return
        self._rec_tab = tab
        self._rec_name = name
        self._rec_duration = duration
        self._reset_hist(idx)
        self._rec_t0 = time.monotonic()
        self.settings.irf_record_s = duration
        self.irf_panel.set_recording(True, 0.0)
        if self._rec_timer is None:
            self._rec_timer = QTimer(self)
            self._rec_timer.timeout.connect(self._tick_irf_record)
        self._rec_timer.start(200)

    def _tick_irf_record(self) -> None:
        if self._rec_tab is None:
            self._rec_timer.stop()
            return
        frac = (time.monotonic() - self._rec_t0) / self._rec_duration
        self.irf_panel.set_recording(True, min(frac, 1.0))
        if frac >= 1.0:
            self._rec_timer.stop()
            self._finish_irf_record()

    def _cancel_irf_record(self) -> None:
        if self._rec_timer is not None:
            self._rec_timer.stop()
        self._rec_tab = None
        self.irf_panel.set_recording(False, 0.0)
        self.statusBar().showMessage("IRF recording cancelled", 5000)

    def _finish_irf_record(self) -> None:
        tab = self._rec_tab
        self._rec_tab = None
        self.irf_panel.set_recording(False, 1.0)
        if tab is None or tab.bins is None or tab.bins.sum() == 0:
            self._warn("No counts were recorded; IRF not saved.")
            return
        h = tab.config.get_settings()
        meta = {
            "device": self.device.name,
            "duration_s": self._rec_duration,
            "ref_ch": h.ref_ch,
            "meas_ch": h.meas_ch,
            "sync_rate_hz": self.counters.sync_rate(),
            "meas_rate_hz": self.counters.rates.get(h.meas_ch, 0.0),
            "thresholds_v": dict(self.settings.thresholds_v),
            "delays_ns": dict(self.settings.delays_ns),
            "multi_hit": h.multi_hit,
        }
        irf = IRF(
            name=self._rec_name, bins=tab.bins.astype(np.int64), left_ns=tab.left_ns, bin_ns=tab.bin_ns, meta=meta
        )
        path = self.irf_store.save(irf)
        self.statusBar().showMessage(f"IRF saved: {path}", 10000)
        self._set_active_irf(irf.name)
        self.irf_panel.name_edit.clear()

    # --------------------------------------------------------------- misc
    def _warn(self, text: str) -> None:
        log.warning(text)
        self.statusBar().showMessage(text, 10000)

    def _on_connection_lost(self, text: str) -> None:
        for tab in self.tabs:
            tab.config.set_enabled_state(False)
            tab.config.show_status("connection lost", warning=True)
        QMessageBox.critical(
            self,
            "Connection lost",
            text + "\nRestart the application after checking the USB cable and the power supply.",
        )

    def _save_settings(self) -> None:
        try:
            self.settings.save()
        except Exception as exc:
            log.warning("settings save failed: %s", exc)

    def closeEvent(self, event) -> None:
        self._closing = True
        if self._rec_timer is not None:
            self._rec_timer.stop()
        for tab in self.tabs:
            h = tab.config.get_settings()
            self.settings.set_hist(tab.idx, h)
        for ch, row in self.channels.rows.items():
            self.settings.thresholds_v[str(ch)] = row.threshold()
            self.settings.delays_ns[str(ch)] = row.delay_ns()
        try:
            self.settings.window_geometry = bytes(self.saveGeometry().toHex()).decode()
        except Exception:
            pass
        self._save_settings()
        try:
            self.device.close()
        except Exception as exc:
            log.warning("device close: %s", exc)
        self._fit_thread.quit()
        self._fit_thread.wait(2000)
        super().closeEvent(event)
