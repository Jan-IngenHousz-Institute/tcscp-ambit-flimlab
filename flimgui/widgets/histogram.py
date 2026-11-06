"""Histogram plot (data + IRF + fit + residuals) and histogram configuration panel."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ..settings import HistSettings

LOG_FLOOR = 0.5


class HistogramView(QWidget):
    def __init__(self):
        super().__init__()
        self.log = True
        self._t = None
        self._counts = None
        self._irf = None
        self._model = None
        self._resid = None
        self._mask = None

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.gl = pg.GraphicsLayoutWidget()
        lay.addWidget(self.gl)
        self.p_main = self.gl.addPlot(row=0, col=0)
        self.p_res = self.gl.addPlot(row=1, col=0)
        self.gl.ci.layout.setRowStretchFactor(0, 4)
        self.gl.ci.layout.setRowStretchFactor(1, 1)
        self.p_res.setXLink(self.p_main)
        self.p_main.setLabel("left", "counts")
        self.p_main.showGrid(x=True, y=True, alpha=0.25)
        self.p_res.setLabel("left", "w. resid.")
        self.p_res.setLabel("bottom", "time", units="ns")
        self.p_res.showGrid(x=True, y=True, alpha=0.25)
        self.p_main.addLegend(offset=(-10, 10))

        self.c_data = self.p_main.plot(pen=pg.mkPen("#4a90d9", width=1), name="data")
        self.c_irf = self.p_main.plot(pen=pg.mkPen("#2e9e4f", width=1, style=Qt.DashLine), name="IRF (scaled)")
        self.c_fit = self.p_main.plot(pen=pg.mkPen("#d9534f", width=2), name="fit", connect="finite")
        self.c_res = self.p_res.plot(pen=None, symbol="o", symbolSize=3, symbolBrush="#555", symbolPen=None)
        self.p_res.addLine(y=0, pen=pg.mkPen("#999"))
        self.set_log(True)

    def set_log(self, on: bool) -> None:
        self.log = bool(on)
        self.p_main.setLogMode(x=False, y=self.log)
        self._redraw()

    def _y(self, v):
        if v is None:
            return None
        v = np.asarray(v, float)
        return np.where(v > 0, v, LOG_FLOOR) if self.log else v

    def set_data(self, t_ns: np.ndarray, counts: np.ndarray) -> None:
        self._t = np.asarray(t_ns, float)
        self._counts = np.asarray(counts, float)
        self.c_data.setData(self._t, self._y(self._counts))
        if self._irf is not None:
            self._draw_irf()

    def set_irf(self, irf_counts: np.ndarray | None) -> None:
        self._irf = None if irf_counts is None else np.asarray(irf_counts, float)
        self._draw_irf()

    def _draw_irf(self) -> None:
        if self._irf is None or self._t is None or self._irf.size != self._t.size:
            self.c_irf.setData([], [])
            return
        peak_d = np.nanmax(self._counts) if self._counts is not None and self._counts.size else 1.0
        peak_i = self._irf.max() if self._irf.size else 1.0
        scale = peak_d / peak_i if peak_i > 0 else 1.0
        self.c_irf.setData(self._t, self._y(self._irf * scale))

    def set_fit(self, model: np.ndarray | None, resid_t_ns: np.ndarray | None, resid: np.ndarray | None) -> None:
        """model is on the data grid; residuals come with their own (possibly coarser) axis."""
        self._model, self._resid, self._mask = model, resid, resid_t_ns
        self._draw_fit()

    def _draw_fit(self) -> None:
        if self._model is None or self._t is None or self._model.size != self._t.size:
            self.c_fit.setData([], [])
            self.c_res.setData([], [])
            return
        m = np.asarray(self._model, float)
        y = self._y(m)
        y = np.where(np.isnan(m), np.nan, y)
        self.c_fit.setData(self._t, y)
        if self._resid is not None and self._mask is not None and len(self._resid) == len(self._mask):
            self.c_res.setData(np.asarray(self._mask, float), np.asarray(self._resid, float))
        else:
            self.c_res.setData([], [])

    def clear_fit(self) -> None:
        self.set_fit(None, None, None)

    def _redraw(self) -> None:
        if self._t is not None and self._counts is not None:
            self.c_data.setData(self._t, self._y(self._counts))
        self._draw_irf()
        self._draw_fit()


class HistogramConfigPanel(QGroupBox):
    enable_toggled = Signal(bool)
    reset_requested = Signal()
    apply_requested = Signal()
    log_toggled = Signal(bool)

    def __init__(self, channel_names: list[str], title: str = "Histogram"):
        super().__init__(title)
        self.channel_names = channel_names
        form = QFormLayout()
        self.ref_combo = QComboBox()
        self.meas_combo = QComboBox()
        for i, n in enumerate(channel_names):
            self.ref_combo.addItem(n, i)
            if i != 0:
                self.meas_combo.addItem(n, i)
        self.ref_combo.setToolTip("Reference channel. The SYNC input can only be a reference.")
        form.addRow("Reference", self.ref_combo)
        form.addRow("Measure", self.meas_combo)

        self.left_spin = QDoubleSpinBox()
        self.left_spin.setRange(-322e6, 322e6)
        self.left_spin.setDecimals(3)
        self.left_spin.setSuffix(" ns")
        self.right_spin = QDoubleSpinBox()
        self.right_spin.setRange(-322e6, 322e6)
        self.right_spin.setDecimals(3)
        self.right_spin.setSuffix(" ns")
        form.addRow("Left bound", self.left_spin)
        form.addRow("Right bound", self.right_spin)
        self.bin_spin = QDoubleSpinBox()
        self.bin_spin.setRange(0.0, 1e9)
        self.bin_spin.setDecimals(1)
        self.bin_spin.setSuffix(" ps")
        self.bin_spin.setSpecialValueText("smallest")
        self.bin_spin.setToolTip("0 = the smallest bin width the hardware allows for this range.")
        form.addRow("Bin width", self.bin_spin)
        self.refresh_spin = QDoubleSpinBox()
        self.refresh_spin.setRange(0.01, 60.0)
        self.refresh_spin.setDecimals(2)
        self.refresh_spin.setSuffix(" s")
        form.addRow("Refresh", self.refresh_spin)
        self.multi_chk = QCheckBox("Multi-hit (all photons per sync, not only the first)")
        self.acc_chk = QCheckBox("Accumulate (keep integrating until reset)")
        form.addRow(self.multi_chk)
        form.addRow(self.acc_chk)

        btns = QHBoxLayout()
        self.apply_btn = QPushButton("Apply")
        self.apply_btn.setToolTip("Send bounds and channels to the device (histogram must be stopped).")
        self.enable_btn = QPushButton("Enable")
        self.enable_btn.setCheckable(True)
        self.reset_btn = QPushButton("Reset")
        btns.addWidget(self.apply_btn)
        btns.addWidget(self.enable_btn)
        btns.addWidget(self.reset_btn)
        self.log_chk = QCheckBox("Log Y")
        self.log_chk.setChecked(True)

        self.bounds_lbl = QLabel("device bounds: --")
        self.bounds_lbl.setStyleSheet("color: #555;")
        self.status_lbl = QLabel("stopped")
        self.status_lbl.setWordWrap(True)

        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addLayout(btns)
        lay.addWidget(self.log_chk)
        lay.addWidget(self.bounds_lbl)
        lay.addWidget(self.status_lbl)

        self.apply_btn.clicked.connect(self.apply_requested)
        self.enable_btn.toggled.connect(self._on_enable)
        self.reset_btn.clicked.connect(self.reset_requested)
        self.log_chk.toggled.connect(self.log_toggled)
        self._config_widgets = [
            self.ref_combo,
            self.meas_combo,
            self.left_spin,
            self.right_spin,
            self.bin_spin,
            self.multi_chk,
            self.acc_chk,
            self.refresh_spin,
            self.apply_btn,
        ]

    def _on_enable(self, on: bool) -> None:
        self.enable_btn.setText("Disable" if on else "Enable")
        self.enable_toggled.emit(on)

    def set_enabled_state(self, on: bool) -> None:
        """Reflect the running state without emitting the toggle signal."""
        self.enable_btn.blockSignals(True)
        self.enable_btn.setChecked(on)
        self.enable_btn.setText("Disable" if on else "Enable")
        self.enable_btn.blockSignals(False)
        self.set_running(on)

    def set_running(self, running: bool) -> None:
        for w in self._config_widgets:
            w.setEnabled(not running)

    def get_settings(self) -> HistSettings:
        return HistSettings(
            ref_ch=int(self.ref_combo.currentData()),
            meas_ch=int(self.meas_combo.currentData()),
            left_ns=self.left_spin.value(),
            right_ns=self.right_spin.value(),
            bin_ps=self.bin_spin.value(),
            multi_hit=self.multi_chk.isChecked(),
            accumulate=self.acc_chk.isChecked(),
            refresh_s=self.refresh_spin.value(),
            enabled=self.enable_btn.isChecked(),
        )

    def set_settings(self, h: HistSettings) -> None:
        i = self.ref_combo.findData(h.ref_ch)
        self.ref_combo.setCurrentIndex(max(i, 0))
        i = self.meas_combo.findData(h.meas_ch)
        self.meas_combo.setCurrentIndex(max(i, 0))
        self.left_spin.setValue(h.left_ns)
        self.right_spin.setValue(h.right_ns)
        self.bin_spin.setValue(h.bin_ps)
        self.multi_chk.setChecked(h.multi_hit)
        self.acc_chk.setChecked(h.accumulate)
        self.refresh_spin.setValue(h.refresh_s)

    def show_bounds(self, left_ns: float, right_ns: float, bin_ps: float, n: int) -> None:
        self.bounds_lbl.setText(f"device bounds: {left_ns:.3f} to {right_ns:.3f} ns, bin {bin_ps:.1f} ps, {n} bins")

    def show_status(self, text: str, warning: bool = False) -> None:
        self.status_lbl.setStyleSheet("color: #b00000; font-weight: bold;" if warning else "")
        self.status_lbl.setText(text)
