"""Persistent user settings, stored as JSON next to the repository root."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
SETTINGS_PATH = ROOT / "flimgui_settings.json"
IRF_DIR = ROOT / "irf-library"


@dataclass
class HistSettings:
    ref_ch: int = 0
    meas_ch: int = 1
    left_ns: float = 0.0
    right_ns: float = 12.5  # one period of the laser at 80 MHz
    bin_ps: float = 0.0  # 0 = smallest the hardware allows
    multi_hit: bool = False
    accumulate: bool = True
    refresh_s: float = 0.5
    enabled: bool = False

    @classmethod
    def from_dict(cls, d: dict) -> HistSettings:
        known = {k: d[k] for k in cls.__dataclass_fields__ if k in d}
        return cls(**known)


@dataclass
class Settings:
    thresholds_v: dict[str, float] = field(default_factory=dict)  # key: channel index as str
    delays_ns: dict[str, float] = field(default_factory=dict)
    histograms: list[dict] = field(default_factory=list)
    active_irf: str = ""
    n_exp: int = 1
    fit_enabled: bool = True
    log_y: bool = True
    irf_record_s: float = 10.0
    window_geometry: str = ""
    sim_mode: bool = False

    def hist(self, idx: int) -> HistSettings:
        if idx < len(self.histograms):
            return HistSettings.from_dict(self.histograms[idx])
        h = HistSettings(meas_ch=min(idx + 1, 2))
        return h

    def set_hist(self, idx: int, h: HistSettings) -> None:
        while len(self.histograms) <= idx:
            self.histograms.append(asdict(HistSettings()))
        self.histograms[idx] = asdict(h)

    @classmethod
    def load(cls, path: Path = SETTINGS_PATH) -> Settings:
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except FileNotFoundError:
            return cls()
        except Exception as exc:  # corrupt file: start fresh, keep a copy
            log.warning("settings unreadable (%s); starting with defaults", exc)
            try:
                Path(path).rename(Path(path).with_suffix(".bad.json"))
            except OSError:
                pass
            return cls()
        known = {k: data[k] for k in cls.__dataclass_fields__ if k in data}
        return cls(**known)

    def save(self, path: Path = SETTINGS_PATH) -> None:
        tmp = Path(path).with_suffix(".tmp")
        tmp.write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")
        tmp.replace(path)
