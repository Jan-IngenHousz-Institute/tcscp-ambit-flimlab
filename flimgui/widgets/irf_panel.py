"""Record, store and select instrument response functions."""

from __future__ import annotations

import time

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox, QHBoxLayout,
                               QLabel, QLineEdit, QProgressBar, QPushButton, QVBoxLayout)

from ..irf import IRF


class IRFPanel(QGroupBox):
    record_requested = Signal(float, str)   # duration s, name
    cancel_requested = Signal()
    active_changed = Signal(str)            # "" = none
    delete_requested = Signal(str)

    def __init__(self):
        super().__init__("Instrument response (IRF)")
        form = QFormLayout()
        self.dur_spin = QDoubleSpinBox()
        self.dur_spin.setRange(1.0, 3600.0)
        self.dur_spin.setValue(10.0)
        self.dur_spin.setSuffix(" s")
        form.addRow("Record for", self.dur_spin)
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("name (default: date and time)")
        form.addRow("Name", self.name_edit)

        rec = QHBoxLayout()
        self.rec_btn = QPushButton("Record IRF")
        self.rec_btn.setToolTip("Put a scatterer in place of the sample, remove the emission "
                                "filter, keep the rate below 500 kcps, then record.")
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        rec.addWidget(self.rec_btn)
        rec.addWidget(self.cancel_btn)
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setTextVisible(True)

        sel = QHBoxLayout()
        self.combo = QComboBox()
        self.combo.addItem("(none)", "")
        self.use_btn = QPushButton("Use")
        self.del_btn = QPushButton("Delete")
        sel.addWidget(self.combo, 1)
        sel.addWidget(self.use_btn)
        sel.addWidget(self.del_btn)

        self.active_lbl = QLabel("active IRF: none (tail fit will be used)")
        self.active_lbl.setWordWrap(True)

        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addLayout(rec)
        lay.addWidget(self.progress)
        lay.addWidget(QLabel("Stored IRFs"))
        lay.addLayout(sel)
        lay.addWidget(self.active_lbl)

        self.rec_btn.clicked.connect(self._on_record)
        self.cancel_btn.clicked.connect(self.cancel_requested)
        self.use_btn.clicked.connect(lambda: self.active_changed.emit(self.combo.currentData() or ""))
        self.del_btn.clicked.connect(self._on_delete)

    def _on_record(self) -> None:
        name = self.name_edit.text().strip() or time.strftime("irf_%Y%m%d_%H%M%S")
        self.record_requested.emit(self.dur_spin.value(), name)

    def _on_delete(self) -> None:
        name = self.combo.currentData() or ""
        if name:
            self.delete_requested.emit(name)

    def set_names(self, names: list[str], active: str = "") -> None:
        self.combo.blockSignals(True)
        self.combo.clear()
        self.combo.addItem("(none)", "")
        for n in names:
            self.combo.addItem(n, n)
        i = self.combo.findData(active)
        self.combo.setCurrentIndex(max(i, 0))
        self.combo.blockSignals(False)

    def set_recording(self, on: bool, fraction: float = 0.0) -> None:
        self.rec_btn.setEnabled(not on)
        self.cancel_btn.setEnabled(on)
        self.progress.setValue(int(round(fraction * 100)))
        self.progress.setFormat("recording %p%" if on else "%p%")

    def set_active_info(self, irf: IRF | None, resampled: bool = False, note: str = "") -> None:
        if irf is None:
            self.active_lbl.setText("active IRF: none (tail fit will be used)")
            return
        created = time.strftime("%Y-%m-%d %H:%M", time.localtime(irf.created)) if irf.created else "?"
        txt = (f"active IRF: <b>{irf.name}</b><br>{irf.total():,} counts, {irf.n} bins of "
               f"{irf.bin_ns * 1000:.1f} ps, from {irf.left_ns:.2f} ns, recorded {created}")
        if resampled:
            txt += "<br><span style='color:#b06000'>resampled onto the current histogram axis</span>"
        if note:
            txt += f"<br>{note}"
        self.active_lbl.setText(txt)
