"""Distributed span loads: lift ``l(y)`` and torsion ``q(y)``.

The loading is sampled at ``n`` points ``y`` between 0 and the span ``L`` and
read from a JSON file::

    {
      "Loads": {
        "y": [0.0, ..., L],
        "l": [...],   # lift per unit span (acts along global Z / local y)
        "q": [...]    # torsion per unit span (about global Y / local x)
      }
    }

:func:`load_loads` returns a :class:`Loads` whose callables ``l(y)`` and
``q(y)`` linearly interpolate the samples (clamped to the end values outside
``[y[0], y[-1]]``). ``y`` may be a scalar or an array.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Loads:
    """Piecewise-linear lift and torsion distributions along the span.

    Parameters
    ----------
    ys : (n,) sample spanwise coordinates, ascending.
    ls : (n,) lift-per-unit-span samples at ``ys``.
    qs : (n,) torsion-per-unit-span samples at ``ys``.
    """

    ys: np.ndarray
    ls: np.ndarray
    qs: np.ndarray

    def l(self, y):
        """Lift per unit span at ``y`` (linear interpolation)."""
        return np.interp(y, self.ys, self.ls)

    def q(self, y):
        """Torsion per unit span at ``y`` (linear interpolation)."""
        return np.interp(y, self.ys, self.qs)


def load_loads(path: str | Path) -> Loads:
    """Load lift/torsion distributions from a JSON file."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))["Loads"]
    ys = np.asarray(data["y"], dtype=float)
    ls = np.asarray(data["l"], dtype=float)
    qs = np.asarray(data["q"], dtype=float)
    if not (ys.shape == ls.shape == qs.shape) or ys.ndim != 1:
        raise ValueError("'y', 'l', 'q' must be 1-D arrays of equal length")
    order = np.argsort(ys)  # np.interp requires ascending sample points
    return Loads(ys[order], ls[order], qs[order])
