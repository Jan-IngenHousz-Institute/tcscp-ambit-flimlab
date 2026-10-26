"""Per-channel threshold and delay controls."""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (QCheckBox, QDoubleSpinBox, QGridLayout, QGroupBox, QLabel,
                               QPushButton, QSlider, QVBoxLayout, QWidget)


class ChannelRow(QWidget):
    threshold_changed = Signal(int, float)   # channel, volts
    delay_changed = Signal(int, float)       # channel, nanoseconds
    auto_requested = Signal(int)

    def __init__(self, ch: int, name: str, thr_range: tuple[float, float],
                 delay_range_ns: tuple[float, float]):
        super().__init__()
        self.ch = ch
        self._thr_range = thr_range
        self._delay_range = delay_range_ns
        self._debounce_thr = QTimer(self, singleShot=True, interval=150)
        self._debounce_dly = QTimer(self, singleShot=True, interval=150)
        self._debounce_thr.timeout.connect(lambda: self.threshold_changed.emit(self.ch, self.threshold()))
        self._debounce_dly.timeout.connect(lambda: self.delay_changed.emit(self.ch, self.delay_ns()))

        grid = QGridLayout(self)
        grid.setContentsMargins(0, 0, 0, 0)
        title = QLabel(f"<b>{name}</b>")
        grid.addWidget(title, 0, 0, 1, 3)

        grid.addWidget(QLabel("Threshold"), 1, 0)
        self.thr_spin = QDoubleSpinBox()
        self.thr_spin.setRange(*thr_range)
        self.thr_spin.setDecimals(3)
        self.thr_spin.setSingleStep(0.01)
        self.thr_spin.setSuffix(" V")
        self.thr_slider = QSlider(Qt.Horizontal)
        self.thr_slider.setRange(int(thr_range[0] * 1000), int(thr_range[1] * 1000))
        grid.addWidget(self.thr_spin, 1, 1)
        grid.addWidget(self.thr_slider, 1, 2)

        grid.addWidget(QLabel("Delay"), 2, 0)
        self.dly_spin = QDoubleSpinBox()
        self.dly_spin.setRange(*delay_range_ns)
        self.dly_spin.setDecimals(3)
        self.dly_spin.setSingleStep(0.1)
        self.dly_spin.setSuffix(" ns")
        self.dly_slider = QSlider(Qt.Horizontal)
        self.dly_slider.setRange(int(delay_range_ns[0] * 100), int(delay_range_ns[1] * 100))  # 10 ps steps
        grid.addWidget(self.dly_spin, 2, 1)
        grid.addWidget(self.dly_slider, 2, 2)

        self.auto_btn = QPushButton("Auto threshold")
        self.auto_btn.setToolTip("Let the Felix scan the input and put the threshold in the middle "
                                 "of the pulse amplitude.")
        grid.addWidget(self.auto_btn, 3, 1, 1, 2)
        grid.setColumnStretch(2, 1)

        self.thr_spin.valueChanged.connect(self._thr_from_spin)
        self.thr_slider.valueChanged.connect(self._thr_from_slider)
        self.dly_spin.valueChanged.connect(self._dly_from_spin)
        self.dly_slider.valueChanged.connect(self._dly_from_slider)
        self.auto_btn.clicked.connect(lambda: self.auto_requested.emit(self.ch))
        self._block = False

    # values ---------------------------------------------------------------
    def threshold(self) -> float:
        return self.thr_spin.value()

    def delay_ns(self) -> float:
        return self.dly_spin.value()

    def set_threshold(self, volts: float) -> None:
        self._block = True
        self.thr_spin.setValue(volts)
        self.thr_slider.setValue(int(round(volts * 1000)))
        self._block = False

    def set_delay_ns(self, ns: float) -> None:
        self._block = True
        self.dly_spin.setValue(ns)
        self.dly_slider.setValue(int(round(ns * 100)))
        self._block = False

    # internal sync ----------------------------------------------------------
    def _thr_from_spin(self, v: float) -> None:
        if self._block:
            return
        self._block = True
        self.thr_slider.setValue(int(round(v * 1000)))
        self._block = False
        self._debounce_thr.start()

    def _thr_from_slider(self, mv: int) -> None:
        if self._block:
            return
        self._block = True
        self.thr_spin.setValue(mv / 1000.0)
        self._block = False
        self._debounce_thr.start()

    def _dly_from_spin(self, v: float) -> None:
        if self._block:
            return
        self._block = True
        self.dly_slider.setValue(int(round(v * 100)))
        self._block = False
        self._debounce_dly.start()

    def _dly_from_slider(self, v: int) -> None:
        if self._block:
            return
        self._block = True
        self.dly_spin.setValue(v / 100.0)
        self._block = False
        self._debounce_dly.start()


class ChannelPanel(QGroupBox):
    threshold_changed = Signal(int, float)
    delay_changed = Signal(int, float)
    auto_requested = Signal(int)

    def __init__(self, names: list[str], thr_ranges: dict[int, tuple[float, float]],
                 delay_ranges_ns: dict[int, tuple[float, float]]):
        super().__init__("Channel configuration")
        self.rows: dict[int, ChannelRow] = {}
        lay = QVBoxLayout(self)
        for ch, name in enumerate(names):
            row = ChannelRow(ch, name, thr_ranges.get(ch, (0.0, 2.5)),
                             delay_ranges_ns.get(ch, (-614.0, 614.0)))
            row.threshold_changed.connect(self.threshold_changed)
            row.delay_changed.connect(self.delay_changed)
            row.auto_requested.connect(self.auto_requested)
            self.rows[ch] = row
            lay.addWidget(row)
        self.couple_chk = QCheckBox("Couple detector channels (same threshold and delay)")
        lay.addWidget(self.couple_chk)

    @property
    def coupled(self) -> bool:
        return self.couple_chk.isChecked()
