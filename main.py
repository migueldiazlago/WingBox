"""Solve a cantilevered wing under distributed lift and torsion.

The wing runs along global Y (the WingBox convention), so the local frame is
local-x = global Y, local-y = global Z (lift), local-z = global X (chord).
The root (node 0) is always clamped.
"""

import numpy as np

from wingbox import solve_wing, plot_deformation, load_stations

WING = "wings/pc24_wing_sections.json"
LOADS = "wings/pc24_loads.json"


def main():
    sol = solve_wing(WING, LOADS)
    coords = np.array([s["EC"] for s in load_stations(WING)])
    print(f"Solved {len(coords)} nodes / {(len(coords) - 1) // 2} element(s)")

    print(f"Tip vertical deflection (global Z): {sol.displacements[-1, 2]:.4f} m")
    print(f"Tip twist (about global Y):         {np.degrees(sol.rotations[-1, 1]):+.4f} deg")
    print(f"Root reaction force  [N]:   {np.round(sol.root_force, 1)}")
    print(f"Root reaction moment [N.m]: {np.round(sol.root_moment, 1)}")

    plot_deformation(coords, sol, scale=3.0)


if __name__ == "__main__":
    main()
