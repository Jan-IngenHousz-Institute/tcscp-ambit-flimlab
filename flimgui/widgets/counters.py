"""Count-rate tiles with a 15 s sparkline and SPAD protection colours."""

from __future__ import annotations

import time
from collections import deque

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QFrame, QHBoxLayout, QLabel, QVBoxLayout, QWidget

SPAD_WARN_HZ = 1.0e6      # amber: pile-up territory and getting close to the limit
SPAD_ALARM_HZ = 3.0e6     # red + beep: rated maximum count rate of the SPAD
PILEUP_WARN_RATIO = 0.05  # detector / sync above 5 % distorts the decay


def fmt_rate(hz: float) -> str:
    if hz is None or hz != hz:
        return "--"
    if hz >= 1e6:
        return f"{hz / 1e6:.3f} Mcps"
    if hz >= 1e3:
        return f"{hz / 1e3:.1f} kcps"
    return f"{hz:.0f} cps"


class CounterTile(QFrame):
    def __init__(self, name: str, warn_hz: float | None, alarm_hz: float | None, span_s: float = 15.0):
        super().__init__()
        self.name = name
        self.warn_hz = warn_hz
        self.alarm_hz = alarm_hz
        self.span_s = span_s
        self.history: deque[tuple[float, float]] = deque()
        self.level = "ok"
        self._last_beep = 0.0

        self.setFrameShape(QFrame.StyledPanel)
        self.setObjectName("tile")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 6)
        top = QHBoxLayout()
        self.name_lbl = QLabel(name)
        self.name_lbl.setStyleSheet("font-weight: bold; font-size: 13px;")
        self.value_lbl = QLabel("--")
        self.value_lbl.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.value_lbl.setStyleSheet("font-size: 22px; font-weight: bold;")
        top.addWidget(self.name_lbl)
        top.addStretch(1)
        top.addWidget(self.value_lbl)
        lay.addLayout(top)

        self.plot = pg.PlotWidget()
        self.plot.setFixedHeight(56)
        self.plot.hideAxis("left")
        self.plot.hideAxis("bottom")
        self.plot.setMouseEnabled(x=False, y=False)
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.plot.setXRange(-span_s, 0, padding=0)
        self.curve = self.plot.plot(pen=pg.mkPen("#4a90d9", width=2))
        lay.addWidget(self.plot)
        self._apply_level()

    def rate(self) -> float:
        return self.history[-1][1] if self.history else float("nan")

    def set_rate(self, hz: float) -> None:
        now = time.monotonic()
        self.history.append((now, hz))
        while self.history and now - self.history[0][0] > self.span_s:
            self.history.popleft()
        xs = np.array([t - now for t, _ in self.history])
        ys = np.array([r for _, r in self.history])
        self.curve.setData(xs, ys)
        self.value_lbl.setText(fmt_rate(hz))

        level = "ok"
        if self.alarm_hz is not None and hz >= self.alarm_hz:
            level = "alarm"
        elif self.warn_hz is not None and hz >= self.warn_hz:
            level = "warn"
        if level != self.level:
            self.level = level
            self._apply_level()
        if level == "alarm" and now - self._last_beep > 5.0:
            QApplication.beep()
            self._last_beep = now

    def _apply_level(self) -> None:
        colours = {"ok": "", "warn": "background-color: #ffe9a8;", "alarm": "background-color: #ff9a9a;"}
        self.setStyleSheet(f"QFrame#tile {{ border: 1px solid #999; border-radius: 6px; {colours[self.level]} }}")
        self.name_lbl.setToolTip(
            {"ok": "",
             "warn": "Above 1 Mcps: pile-up distorts the decay. Reduce the light reaching the SPAD.",
             "alarm": "Above 3 Mcps: the SPAD is at its rated maximum. Reduce the light now."}[self.level])


class CounterPanel(QWidget):
    """One tile per channel plus a pile-up indicator for the measurement channel."""

    def __init__(self, channel_names: list[str], sync_ch: int = 0):
        super().__init__()
        self.sync_ch = sync_ch
        self.tiles: dict[int, CounterTile] = {}
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        for ch, name in enumerate(channel_names):
            if ch == sync_ch:
                tile = CounterTile(name, None, None)
            else:
                tile = CounterTile(name, SPAD_WARN_HZ, SPAD_ALARM_HZ)
            self.tiles[ch] = tile
            lay.addWidget(tile)
        self.pileup_lbl = QLabel("pile-up: --")
        self.pileup_lbl.setToolTip("Detector rate divided by sync rate. Keep it below 5 % for a "
                                   "clean decay (single-hit histogram).")
        lay.addWidget(self.pileup_lbl)
        self.rates: dict[int, float] = {}

    def update_rates(self, rates: dict[int, float]) -> None:
        self.rates = dict(rates)
        for ch, hz in rates.items():
            if ch in self.tiles:
                self.tiles[ch].set_rate(hz)
        sync = rates.get(self.sync_ch, 0.0)
        others = [hz for ch, hz in rates.items() if ch != self.sync_ch]
        if sync > 0 and others:
            ratio = max(others) / sync
            colour = "#b00000" if ratio > PILEUP_WARN_RATIO else ""
            self.pileup_lbl.setStyleSheet(f"color: {colour}; font-weight: bold;" if colour else "")
            self.pileup_lbl.setText(f"pile-up: {ratio * 100:.2f} % of sync rate"
                                    + ("  (too high, lower the light)" if ratio > PILEUP_WARN_RATIO else ""))
        else:
            self.pileup_lbl.setStyleSheet("")
            self.pileup_lbl.setText("pile-up: -- (no sync)")

    def sync_rate(self) -> float:
        return self.rates.get(self.sync_ch, 0.0)
