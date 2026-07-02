"""3D plot of the wing mesh, its bending deformation and its twist (matplotlib).

The bending is shown by the deformed centre-line; the torsion is shown by
short chord "ribs" drawn at every station, tilted about the span axis by the
section twist (at the same magnification as the deflection) and coloured by
the local twist. Bending and twist share the single ``scale`` factor.
"""

from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from wingbox.assemble import Solution, InternalForces


def plot_internal_forces(intf: InternalForces, save: str | None = None):
    """Plot the shear ``V(y)``, bending moment ``M(y)`` and torque ``T(y)``.

    Parameters
    ----------
    intf : the :class:`~wingbox.assemble.InternalForces` diagrams to draw.
    save : path to save the figure; if ``None`` the figure is shown.
    """
    panels = [
        (intf.V / 1e3, "shear  $V(y)$  [kN]", "#1f77b4"),
        (intf.M / 1e3, "bending moment  $M(y)$  [kN·m]", "#d62728"),
        (intf.T / 1e3, "torque  $T(y)$  [kN·m]", "#2ca02c"),
    ]
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(7.0, 7.2))
    for ax, (val, label, color) in zip(axes, panels):
        ax.plot(intf.y, val, color=color, lw=2)
        ax.fill_between(intf.y, val, color=color, alpha=0.15)
        ax.axhline(0, color="0.7", lw=0.8)
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        ax.margins(x=0)
    axes[0].set_title("Spanwise internal-force diagrams (root → tip)")
    axes[-1].set_xlabel("span  $y$  [m]")
    fig.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")
    else:
        plt.show()
    return axes


def _equal_aspect(ax, pts: np.ndarray) -> None:
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    center = (lo + hi) / 2.0
    radius = max((hi - lo).max() / 2.0, 1e-9)
    for set_lim, c in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), center):
        set_lim(c - radius, c + radius)


def plot_deformation(
    coords: np.ndarray,
    sol: Solution,
    scale: float = 1.0,
    chord: float | None = None,
    save: str | None = None,
):
    """Plot the undeformed and deformed wing with a twist visualisation.

    Parameters
    ----------
    coords : (n_nodes, 3) undeformed node positions (station elastic centres).
    sol : the :class:`~wingbox.assemble.Solution` to visualise.
    scale : magnification applied to *both* the translational deflection and
        the section twist (they share the same scale).
    chord : visual rib (chord) length. Defaults to ~12% of the span.
    save : path to save the figure; if ``None`` the figure is shown.
    """
    coords = np.asarray(coords, dtype=float)
    trans, rot = sol.displacements, sol.rotations
    deformed = coords + scale * trans

    span = float(np.ptp(coords[:, 1]))
    if chord is None:
        chord = 0.12 * span

    # Twist = section rotation about the span axis (global Y). The colour scale
    # spans the actual min..max twist across the stations.
    twist = rot[:, 1]
    twist_deg = np.degrees(twist)
    lo, hi = float(twist_deg.min()), float(twist_deg.max())
    if hi - lo < 1e-9:  # avoid a degenerate colour range
        lo, hi = lo - 1e-6, hi + 1e-6
    norm = Normalize(lo, hi)
    cmap = mpl.colormaps["coolwarm"]

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Centre-lines: undeformed (dashed grey) and deformed (solid).
    ax.plot(*coords.T, "--", color="0.7", lw=1.2, label="undeformed")
    ax.plot(*deformed.T, "-", color="0.25", lw=2.0, label="deformed")
    ax.scatter(*deformed.T, s=12, color="0.25", depthshade=False)
    ax.scatter(*coords[0], s=60, marker="s", color="k", label="clamped root")

    # Chord ribs: a chordwise segment (global X) tilted about the span axis
    # by the twist (same scale as the deflection), so the edges rise/fall.
    rib_pts = [coords, deformed]
    half = chord / 2.0
    for i in range(len(coords)):
        a = scale * twist[i]
        off = np.array([half * np.cos(a), 0.0, -half * np.sin(a)])  # rotate about Y
        seg = np.array([deformed[i] - off, deformed[i] + off])
        color = cmap(norm(twist_deg[i]))
        ax.plot(*seg.T, color=color, lw=2.5)
        ax.scatter(*seg[1], s=16, color=color, depthshade=False)  # leading edge
        rib_pts.append(seg)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.10)
    cbar.set_label("twist about span [deg]")

    ax.set_xlabel("X  (chord)")
    ax.set_ylabel("Y  (span)")
    ax.set_zlabel("Z  (lift)")
    ax.set_title(
        f"Wing deformation ×{scale:g}   "
        f"(tip: {deformed[-1, 2]:.3f} m, {twist_deg[-1]:+.3f}°)"
    )
    ax.legend(loc="upper left")
    ax.view_init(elev=22, azim=-62)
    _equal_aspect(ax, np.vstack(rib_pts))
    fig.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")
    else:
        plt.show()
    return ax
