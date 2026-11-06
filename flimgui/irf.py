"""Storage of instrument response functions (IRF) as JSON files."""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .settings import IRF_DIR


@dataclass
class IRF:
    name: str
    bins: np.ndarray  # counts per bin
    left_ns: float
    bin_ns: float
    created: float = field(default_factory=time.time)
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.bins.size)

    @property
    def right_ns(self) -> float:
        return self.left_ns + self.n * self.bin_ns

    @property
    def centers_ns(self) -> np.ndarray:
        return self.left_ns + (np.arange(self.n) + 0.5) * self.bin_ns

    def total(self) -> int:
        return int(self.bins.sum())

    def on_axis(self, left_ns: float, bin_ns: float, n: int) -> tuple[np.ndarray, bool]:
        """Return the IRF resampled onto another histogram axis.

        Returns (irf_counts, exact) where exact is True when no resampling was needed.
        """
        same = abs(left_ns - self.left_ns) < 1e-9 and abs(bin_ns - self.bin_ns) < 1e-6 and n == self.n
        if same:
            return self.bins.astype(float), True
        centers = left_ns + (np.arange(n) + 0.5) * bin_ns
        out = np.interp(centers, self.centers_ns, self.bins.astype(float), left=0.0, right=0.0)
        # keep the integral, so that the scale does not jump with the bin width
        if out.sum() > 0:
            out *= self.bins.sum() * (bin_ns / self.bin_ns) / out.sum() * (self.bin_ns / bin_ns)
        return out, False

    def to_json(self) -> dict:
        return {
            "name": self.name,
            "created": self.created,
            "left_ns": self.left_ns,
            "bin_ns": self.bin_ns,
            "bins": [int(b) for b in self.bins],
            "meta": self.meta,
        }

    @classmethod
    def from_json(cls, d: dict) -> IRF:
        return cls(
            name=d["name"],
            bins=np.asarray(d["bins"], dtype=np.int64),
            left_ns=float(d["left_ns"]),
            bin_ns=float(d["bin_ns"]),
            created=float(d.get("created", 0.0)),
            meta=dict(d.get("meta", {})),
        )


def _safe_name(name: str) -> str:
    name = re.sub(r"[^\w\-. ]+", "_", name.strip())
    return name or "irf"


class IRFStore:
    def __init__(self, directory: Path = IRF_DIR):
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)

    def path(self, name: str) -> Path:
        return self.dir / f"{_safe_name(name)}.json"

    def names(self) -> list[str]:
        return sorted(p.stem for p in self.dir.glob("*.json"))

    def save(self, irf: IRF) -> Path:
        irf.name = _safe_name(irf.name)
        p = self.path(irf.name)
        p.write_text(json.dumps(irf.to_json()), encoding="utf-8")
        return p

    def load(self, name: str) -> IRF:
        return IRF.from_json(json.loads(self.path(name).read_text(encoding="utf-8")))

    def delete(self, name: str) -> None:
        p = self.path(name)
        if p.exists():
            p.unlink()
