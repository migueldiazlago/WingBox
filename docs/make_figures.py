"""Generate the README diagrams (schematics) and the code-output plot.

Run from the project root:  uv run python docs/make_figures.py
Outputs SVG schematics + a PNG of the solver result into docs/.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Arc, Polygon

HERE = os.path.dirname(__file__)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def arrow(ax, tail, head, color="k", lw=1.8, ls="-", mut=14):
    ax.add_patch(
        FancyArrowPatch(
            tail, head, arrowstyle="-|>", mutation_scale=mut,
            lw=lw, ls=ls, color=color, shrinkA=0, shrinkB=0, zorder=5,
        )
    )


def curved_arrow(ax, center, r, a0, a1, color="k", lw=1.8):
    """A circular (moment) arrow from angle a0 to a1 [deg] about ``center``."""
    ax.add_patch(Arc(center, 2 * r, 2 * r, theta1=min(a0, a1), theta2=max(a0, a1),
                     color=color, lw=lw, zorder=5))
    aend = np.radians(a1)
    tip = (center[0] + r * np.cos(aend), center[1] + r * np.sin(aend))
    # tangent direction (sense a0 -> a1)
    s = 1.0 if a1 > a0 else -1.0
    tang = (-s * np.sin(aend), s * np.cos(aend))
    back = (tip[0] - 0.35 * r * tang[0], tip[1] - 0.35 * r * tang[1])
    arrow(ax, back, tip, color=color, lw=lw, mut=12)


def naca4(m, p, t, n=140):
    """NACA 4-digit airfoil upper/lower coordinates (chord = 1)."""
    x = (1 - np.cos(np.linspace(0, np.pi, n))) / 2
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                  + 0.2843 * x**3 - 0.1015 * x**4)
    yc = np.where(x < p, m / p**2 * (2 * p * x - x**2),
                  m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2))
    return x, yc + yt, yc - yt, yc  # x, upper, lower, camber


# --------------------------------------------------------------------------
# 1. wing beam model: side view (lift, BC, bending) + top view (sweep, torque)
# --------------------------------------------------------------------------
def fig_wing_model():
    L = 8.5
    ys = np.linspace(0, L, 9)
    lift = np.sqrt(np.clip(1 - (ys / L) ** 2, 0, None))  # elliptical, normalised

    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.8))

    # ---- (a) side view: span (y) vs lift (z) -----------------------------
    a.axhline(0, color="k", lw=2, zorder=3)                       # beam
    a.plot(ys, 0 * ys, "o", color="k", ms=5, zorder=4)           # nodes
    a.add_patch(Rectangle((-0.9, -1.7), 0.9, 3.4, hatch="////",
                          facecolor="0.85", edgecolor="k", zorder=2))
    a.text(-0.45, -2.15, "clamped root\nall 6 DOF = 0", ha="center",
           va="top", fontsize=9)

    for yi, li in zip(ys, lift):                                  # lift arrows
        if li > 1e-3:
            arrow(a, (yi, 0), (yi, 2.0 * li), color="#1f77b4", lw=1.6)
    a.plot(ys, 2.0 * lift, "--", color="#1f77b4", lw=1.2)
    a.text(L * 0.5, 2.35, "lift  $l(y)$  (elliptical)", color="#1f77b4",
           ha="center", fontsize=10)

    defl = 1.4 * (ys / L) ** 2                                    # deformed shape
    a.plot(ys, -defl, "-", color="0.4", lw=2)
    a.plot(ys, 0 * ys, ":", color="0.6", lw=1)
    a.text(L, -defl[-1] - 0.15, "deformed", color="0.4", ha="right", va="top",
           fontsize=9)

    arrow(a, (0, -2.6), (2, -2.6), color="k", mut=12)
    a.text(2.1, -2.6, "span  $y$", va="center", fontsize=9)
    arrow(a, (0, -2.6), (0, -1.4), color="k", mut=12)
    a.text(0.15, -1.35, "lift  $z$", va="bottom", fontsize=9)
    a.set_title("(a)  side view — cantilever bending under lift")
    a.set_xlim(-1.2, L + 1.2)
    a.set_ylim(-3.0, 2.8)
    a.set_aspect("equal")
    a.axis("off")

    # ---- (b) top view: planform, sweep, elastic vs aero axis, torque -----
    sweep = np.radians(18)  # exaggerated for clarity
    croot, ctip = 2.3, 1.0
    c = croot + (ctip - croot) * ys / L
    xea = ys * np.tan(sweep)          # elastic axis (beam), ~40% chord datum
    le = xea - 0.40 * c               # leading edge (chord measured from EA)
    te = xea + 0.60 * c               # trailing edge
    xac = xea - 0.15 * c              # aero centre (¼-chord ≈ 0.15c ahead of EA)

    b.fill_between(ys, le, te, color="#eaf2fb", zorder=1)
    b.plot(ys, le, "k", lw=1.5)
    b.plot(ys, te, "k", lw=1.5)
    b.plot([ys[0], ys[0]], [le[0], te[0]], "k", lw=1.5)
    b.plot([ys[-1], ys[-1]], [le[-1], te[-1]], "k", lw=1.5)
    b.add_patch(Rectangle((-0.9, le[0] - 0.2), 0.9, (te[0] - le[0]) + 0.4,
                          hatch="////", facecolor="0.85", edgecolor="k", zorder=2))

    b.plot(ys, xea, "-", color="#d62728", lw=2, zorder=4,
           label="elastic axis (beam, through EC)")
    b.plot(ys, xac, "--", color="#1f77b4", lw=1.8, zorder=4,
           label="aerodynamic axis (¼-chord, lift acts here)")

    # offset e -> torque, drawn at a mid station
    j = 5
    b.annotate("", xy=(ys[j], xea[j]), xytext=(ys[j], xac[j]),
               arrowprops=dict(arrowstyle="<->", color="0.2"))
    b.text(ys[j] + 0.12, 0.5 * (xea[j] + xac[j]), "$e$", fontsize=11)
    for k in (3, 6):
        curved_arrow(b, (ys[k], xea[k]), 0.45, 25, 215, color="#d62728", lw=1.8)
    b.text(ys[6] + 0.35, xea[6] - 1.15,
           "torque  $q(y)=l(y)\\,e$", color="#d62728", fontsize=10)

    # sweep angle marker
    b.plot([0, L * 0.5], [xea[0], xea[0]], ":", color="0.5", lw=1)
    b.add_patch(Arc((0, 0), 3.2, 3.2, theta1=0, theta2=np.degrees(sweep),
                    color="0.4"))
    b.text(1.85, 0.28, "$\\Lambda$ (sweep)", fontsize=9, color="0.3")

    arrow(b, (0, te[0] + 0.6), (2, te[0] + 0.6), color="k", mut=12)
    b.text(2.1, te[0] + 0.6, "span  $y$", va="center", fontsize=9)
    arrow(b, (0, te[0] + 0.6), (0, te[0] - 0.4), color="k", mut=12)
    b.text(0.12, te[0] - 0.4, "chord  $x$", va="top", fontsize=9)

    b.set_title("(b)  top view — sweep, and the lift/elastic-axis offset $e$")
    b.legend(loc="lower right", fontsize=8, framealpha=0.95)
    b.set_xlim(-1.2, L + 1.4)
    b.invert_yaxis()  # chord: leading edge up
    b.set_aspect("equal")
    b.axis("off")

    fig.tight_layout()
    fig.savefig(f"{HERE}/wing_model.svg", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# 2. six DOF at a node
# --------------------------------------------------------------------------
def fig_node_dofs():
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    o = np.array([0, 0])
    # isometric-ish unit axes: X (chord), Y (span), Z (lift)
    ex, ey, ez = np.array([-0.85, -0.5]), np.array([1.0, 0.0]), np.array([0, 1.0])

    for e, name, col in [(ex, "X (chord)", "#2ca02c"),
                         (ey, "Y (span)", "#1f77b4"),
                         (ez, "Z (lift)", "#d62728")]:
        arrow(ax, o, 2.4 * e, color=col, lw=2.2)
        ax.text(*(2.95 * e), name, color=col, fontsize=10, ha="center",
                va="center", fontweight="bold")
        # rotation arc near the axis tip
        ang = np.degrees(np.arctan2(e[1], e[0]))
        curved_arrow(ax, 1.5 * e, 0.42, ang + 40, ang + 250, color=col, lw=1.4)

    ax.plot(*o, "o", color="k", ms=10, zorder=6)
    ax.text(0.12, -0.28, "node", fontsize=10)
    ax.text(0, 3.15,
            "6 DOF per node:  3 translations $[t_x,t_y,t_z]$"
            "  +  3 rotations $[r_x,r_y,r_z]$",
            ha="center", fontsize=11)
    ax.text(0, -3.0,
            "straight arrows = translations   •   curved arrows = rotations",
            ha="center", fontsize=9, color="0.4")
    ax.set_xlim(-3.4, 3.4)
    ax.set_ylim(-3.3, 3.4)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(f"{HERE}/node_dofs.svg", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# 3. real wing section: airfoil + wing box + EC/AC + bending axes + torsion
# --------------------------------------------------------------------------
def fig_wing_section():
    x, up, lo, cam = naca4(0.02, 0.4, 0.12)
    fig, ax = plt.subplots(figsize=(9.5, 5.2))

    ax.plot(x, up, "k", lw=1.6)
    ax.plot(x, lo, "k", lw=1.6)
    ax.plot([0, 1], [0, 0], "--", color="0.6", lw=1)  # chord line

    # structural wing box between front (0.15c) and rear (0.60c) spars
    xf, xr = 0.15, 0.60
    def surf(xq, arr):
        return np.interp(xq, x, arr)
    box_x = [xf, xr, xr, xf]
    box_top = [surf(xf, up), surf(xr, up)]
    box_bot = [surf(xr, lo), surf(xf, lo)]
    ax.add_patch(Polygon(list(zip(box_x, box_top + box_bot)),
                         closed=True, facecolor="#eaf2fb", edgecolor="none", zorder=1))
    for xs in (xf, xr):                                   # spar webs
        ax.plot([xs, xs], [surf(xs, lo), surf(xs, up)], color="#1f77b4", lw=4,
                solid_capstyle="butt", zorder=2)
    xx = x[(x >= xf) & (x <= xr)]                         # skins (caps)
    ax.plot(xx, surf(xx, up), color="#1f77b4", lw=4, zorder=2)
    ax.plot(xx, surf(xx, lo), color="#1f77b4", lw=4, zorder=2)

    # aerodynamic centre (¼-chord) and elastic centre / shear axis
    xac, xec = 0.25, 0.42
    ax.plot(xac, 0, "o", color="#1f77b4", ms=9, zorder=6)
    ax.annotate("aerodynamic centre (¼-chord)\nlift $L$ acts here",
                (xac, 0), (xac - 0.02, 0.22), fontsize=9, color="#1f77b4",
                ha="center", arrowprops=dict(arrowstyle="->", color="#1f77b4"))
    ax.plot(xec, surf(xec, cam), "s", color="#d62728", ms=9, zorder=6)
    ax.annotate("elastic centre EC\n(shear axis)", (xec, surf(xec, cam)),
                (xec + 0.22, -0.20), fontsize=9, color="#d62728",
                ha="center", arrowprops=dict(arrowstyle="->", color="#d62728"))

    # offset e between lift line and EC
    ax.annotate("", xy=(xec, 0.0), xytext=(xac, 0.0),
                arrowprops=dict(arrowstyle="<->", color="0.2"))
    ax.text(0.5 * (xac + xec), 0.035, "$e$", fontsize=12, ha="center")

    # bending neutral axes through EC
    zc = surf(xec, cam)
    ax.plot([xf - 0.05, xr + 0.08], [zc, zc], "-.", color="#d62728", lw=1.2)
    ax.text(xr + 0.10, zc, "flap (vertical) bending\nneutral axis  →  $EI_z$",
            fontsize=8.5, va="center", color="#d62728")
    ax.plot([xec, xec], [surf(xec, lo) - 0.05, surf(xec, up) + 0.05],
            "-.", color="#2ca02c", lw=1.2)
    ax.text(xec, surf(xec, up) + 0.07,
            "edge (chordwise) bending → $EI_y$", fontsize=8.5, ha="center",
            color="#2ca02c")

    # box height h dimension
    ax.annotate("", xy=(0.37, surf(0.37, up)), xytext=(0.37, surf(0.37, lo)),
                arrowprops=dict(arrowstyle="<->", color="0.3"))
    ax.text(0.335, 0, "$h$", fontsize=11, color="0.3", ha="right")

    ax.text(0.5, -0.16,
            "closed cell (skins + spar webs) → torsion stiffness $GJ$;"
            "  webs carry transverse shear → shear areas $A_{sy},A_{sz}$",
            ha="center", fontsize=9, color="#1f77b4")
    ax.text(0.0, 0.19, "$x$ = chord  •  $z$ = vertical (lift)  •  $y$ = span (out of page)",
            fontsize=9, color="0.35")

    ax.set_title("Wing section — where the structural properties come from")
    ax.set_xlim(-0.05, 1.15)
    ax.set_ylim(-0.22, 0.26)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(f"{HERE}/wing_section.svg", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# 4. actual solver output
# --------------------------------------------------------------------------
def fig_output():
    from wingbox import (solve_wing, internal_forces, plot_deformation,
                         plot_internal_forces, load_stations, load_loads)

    wing, loads = "wings/pc24_wing_sections.json", "wings/pc24_loads.json"
    sol = solve_wing(wing, loads)
    coords = np.array([s["EC"] for s in load_stations(wing)])
    plot_deformation(coords, sol, scale=3.0, save=f"{HERE}/deformation.png")

    intf = internal_forces(load_loads(loads), span=coords[:, 1].max())
    plot_internal_forces(intf, save=f"{HERE}/internal_forces.png")


if __name__ == "__main__":
    fig_wing_model()
    fig_node_dofs()
    fig_wing_section()
    fig_output()
    print("figures written to docs/")
