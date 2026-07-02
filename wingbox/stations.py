"""Wing stations: combined material/section record + JSON loader.

A wing JSON lists stations along the span, ``EC`` being the section's
elastic centre::

    {
      "Wingsection": {
        "station1": {
          "EC": [x, y, z],
          "E": 70e9, "nu": 0.33,
          "A": 5e-3, "Iy": 4.16e-7, "Iz": 1.04e-6, "J": 1.46e-6,
          "ky": 0.833, "kz": 0.833
        },
        "station2": { ... }
      }
    }

:func:`load_stations` returns a plain list of station dicts, each::

    {"name": str, "EC": (3,) ndarray, "section": Section}

ordered by the trailing integer in the station name.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Section:
    """Material + cross-section properties, in local principal axes.

    Local x runs along the beam; (y, z) are the section principal axes, so
    ``Iz`` bends in the local x-y plane, ``Iy`` in the local x-z plane, and
    ``J`` is the St. Venant torsion constant. ``ky, kz`` are Timoshenko shear
    correction factors (default 5/6).
    """

    E: float
    A: float
    Iy: float
    Iz: float
    J: float
    nu: float = 0.3
    ky: float = 5.0 / 6.0
    kz: float = 5.0 / 6.0

    @property
    def G(self) -> float:
        """Shear modulus ``E / (2 (1 + nu))``."""
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def A_sy(self) -> float:
        """Effective shear area for shear in local y."""
        return self.ky * self.A

    @property
    def A_sz(self) -> float:
        """Effective shear area for shear in local z."""
        return self.kz * self.A

    @classmethod
    def from_dict(cls, data: dict) -> "Section":
        """Build from a mapping, ignoring extra keys (e.g. ``position``)."""
        names = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in names})


def _order(name: str) -> int:
    match = re.search(r"(\d+)$", name)
    return int(match.group(1)) if match else 0


def load_stations(path: str | Path) -> list[dict]:
    """Load wing stations as a list of dicts, ordered along the span."""
    wing = json.loads(Path(path).read_text(encoding="utf-8"))["Wingsection"]
    stations = [
        {
            "name": name,
            "EC": np.asarray(block["EC"], dtype=float),
            "section": Section.from_dict(block),
        }
        for name, block in wing.items()
    ]
    stations.sort(key=lambda s: _order(s["name"]))
    return stations
