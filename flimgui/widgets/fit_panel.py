"""Live fit controls and results."""

from __future__ import annotations

import math

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QCheckBox, QFormLayout, QGroupBox, QLabel, QSpinBox, QTableWidget,
                               QTableWidgetItem, QVBoxLayout)

from ..fitting import FitResult


def _fmt(v: float, nd: int = 3) -> str:
    return "--" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{nd}f}"


class FitPanel(QGroupBox):
    n_exp_changed = Signal(int)
    enabled_changed = Signal(bool)

    def __init__(self):
        super().__init__("Lifetime fit")
        form = QFormLayout()
        self.n_spin = QSpinBox()
        self.n_spin.setRange(1, 4)
        form.addRow("Exponentials", self.n_spin)
        self.enable_chk = QCheckBox("Fit continuously")
        self.enable_chk.setChecked(True)
        form.addRow(self.enable_chk)
        self.mode_lbl = QLabel("mode: --")
        form.addRow(self.mode_lbl)

        self.table = QTableWidget(1, 3)
        self.table.setHorizontalHeaderLabels(["Amplitude", "Fraction", "Lifetime (ns)"])
        self.table.verticalHeader().setVisible(True)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setFixedHeight(140)

        self.tau_amp_lbl = QLabel("--")
        self.tau_amp_lbl.setAlignment(Qt.AlignCenter)
        self.tau_amp_lbl.setStyleSheet("font-size: 26px; font-weight: bold;")
        self.tau_amp_lbl.setToolTip("Amplitude-weighted average lifetime: sum(a*tau)/sum(a)")
        self.tau_int_lbl = QLabel("intensity-weighted: --")
        self.tau_int_lbl.setAlignment(Qt.AlignCenter)
        self.tau_int_lbl.setStyleSheet("color: #555;")
        self.stats_lbl = QLabel("chi2r: --   background: --   shift: --")
        self.stats_lbl.setWordWrap(True)
        self.msg_lbl = QLabel("")
        self.msg_lbl.setWordWrap(True)
        self.msg_lbl.setStyleSheet("color: #555;")

        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addWidget(self.table)
        lay.addWidget(QLabel("Average lifetime (amplitude-weighted)"))
        lay.addWidget(self.tau_amp_lbl)
        lay.addWidget(self.tau_int_lbl)
        lay.addWidget(self.stats_lbl)
        lay.addWidget(self.msg_lbl)

        self.n_spin.valueChanged.connect(self._on_n)
        self.enable_chk.toggled.connect(self.enabled_changed)
        self._on_n(self.n_spin.value())

    def _on_n(self, n: int) -> None:
        self.table.setRowCount(n)
        for r in range(n):
            self.table.setVerticalHeaderItem(r, QTableWidgetItem(f"{r + 1}"))
        self.n_exp_changed.emit(n)

    @property
    def n_exp(self) -> int:
        return self.n_spin.value()

    @property
    def enabled(self) -> bool:
        return self.enable_chk.isChecked()

    def show_result(self, r: FitResult | None) -> None:
        if r is None or r.mode == "none":
            for row in range(self.table.rowCount()):
                for col in range(3):
                    self.table.setItem(row, col, QTableWidgetItem("--"))
            self.tau_amp_lbl.setText("--")
            self.tau_int_lbl.setText("intensity-weighted: --")
            self.stats_lbl.setText("chi2r: --   background: --   shift: --")
            self.mode_lbl.setText("mode: --")
            self.msg_lbl.setText(r.message if r is not None else "")
            return
        self.mode_lbl.setText("mode: reconvolution with IRF" if r.mode == "reconv"
                              else "mode: tail fit (no IRF active)")
        fr = r.fractions
        self.table.setRowCount(r.n_exp)
        for i in range(r.n_exp):
            self.table.setItem(i, 0, QTableWidgetItem(_fmt(float(r.amps[i]), 1)))
            self.table.setItem(i, 1, QTableWidgetItem(f"{fr[i] * 100:.1f} %"))
            self.table.setItem(i, 2, QTableWidgetItem(_fmt(float(r.taus_ns[i]), 3)))
        self.tau_amp_lbl.setText(f"{_fmt(r.tau_amp_ns, 3)} ns")
        self.tau_int_lbl.setText(f"intensity-weighted: {_fmt(r.tau_int_ns, 3)} ns")
        shift = f"{r.shift_ns * 1000:.0f} ps" if r.mode == "reconv" else "n/a"
        self.stats_lbl.setText(f"chi2r: {_fmt(r.chi2r, 2)}   background: {_fmt(r.bg, 2)} /bin   "
                               f"shift: {shift}")
        self.msg_lbl.setText("" if r.success else f"fit did not converge: {r.message}")
