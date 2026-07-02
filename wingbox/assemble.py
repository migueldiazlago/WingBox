"""Assemble and solve the cantilevered-wing FE problem.

Pipeline
--------
1. Read stations and mesh them into 3-node quadratic elements: with
   ``n = 2*n_elem + 1`` stations, element ``e`` owns nodes
   ``[2e, 2e+1, 2e+2]`` (end nodes shared between neighbours).
2. Assemble the global stiffness ``K`` from each element's global stiffness.
3. Integrate the distributed loads into a consistent nodal load vector ``F``:

       F_i = integral over span of  N_i(y) * load(y)  dy

   with lift ``l(y)`` driving the global-Z translation DOF and torsion
   ``q(y)`` the global-Y rotation DOF.
4. Clamp the root (node 0, all 6 DOF) and solve ``K u = F``.

Each element uses a single, constant :class:`Section` taken at its mid
station (piecewise-constant properties along a tapered span).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from wingbox.element import BeamElement, _GAUSS, shape_functions
from wingbox.loads import Loads, load_loads
from wingbox.stations import load_stations

# Per-node global DOF indices that the distributed loads drive.
_LIFT_DOF = 2  # translation along global Z (lift)
_TORSION_DOF = 4  # rotation about global Y (span axis)


@dataclass(frozen=True)
class Solution:
    """Result of a solve, in global axes.

    Attributes
    ----------
    displacements : (n_nodes, 3) nodal translations [tx, ty, tz].
    rotations : (n_nodes, 3) nodal rotations [rx, ry, rz].
    root_force : (3,) reaction force [Fx, Fy, Fz] at the clamped root.
    root_moment : (3,) reaction moment [Mx, My, Mz] at the clamped root.
    """

    displacements: np.ndarray
    rotations: np.ndarray
    root_force: np.ndarray
    root_moment: np.ndarray


def build_elements(stations, ref_vec=(0.0, 0.0, 1.0)) -> list[BeamElement]:
    """Mesh an ordered station list into 3-node quadratic beam elements."""
    n = len(stations)
    if n < 3 or n % 2 == 0:
        raise ValueError(f"need an odd number of stations >= 3, got {n}")

    coords = np.array([s["EC"] for s in stations])
    elements = []
    for e in range((n - 1) // 2):
        node_ids = [2 * e, 2 * e + 1, 2 * e + 2]
        section = stations[2 * e + 1]["section"]  # mid-station properties
        elements.append(BeamElement(node_ids, coords[node_ids], section, ref_vec))
    return elements


def assemble_stiffness(elements, n_nodes: int) -> np.ndarray:
    """Scatter element stiffnesses into the global (6 n_nodes)^2 matrix."""
    ndof = 6 * n_nodes
    K = np.zeros((ndof, ndof))
    for el in elements:
        dofs = el.global_dofs()
        K[np.ix_(dofs, dofs)] += el.stiffness()
    return K


def assemble_loads(elements, loads: Loads, n_nodes: int, ngp: int = 3) -> np.ndarray:
    """Consistent nodal load vector from the distributed lift/torsion."""
    ndof = 6 * n_nodes
    F = np.zeros(ndof)
    pts, wts = _GAUSS[ngp]
    for el in elements:
        y_nodes = el.node_coords[:, 1]  # span coordinate = global Y
        jac = el.length / 2.0
        fe = np.zeros(BeamElement.N_DOF)
        for xi, w in zip(pts, wts):
            N, _ = shape_functions(xi)
            y = float(N @ y_nodes)  # isoparametric map to physical span
            scale = w * jac
            fe[_LIFT_DOF::6] += N * loads.l(y) * scale
            fe[_TORSION_DOF::6] += N * loads.q(y) * scale
        F[el.global_dofs()] += fe
    return F


def solve(K: np.ndarray, F: np.ndarray, fixed) -> np.ndarray:
    """Solve K u = F with the ``fixed`` DOFs constrained to zero."""
    ndof = K.shape[0]
    free = np.setdiff1d(np.arange(ndof), np.asarray(fixed, dtype=int))
    u = np.zeros(ndof)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    return u


def solve_wing(wing_path, loads_path, ref_vec=(0.0, 0.0, 1.0)) -> Solution:
    """End-to-end solve: read files, assemble, clamp the root, solve."""
    stations = load_stations(wing_path)
    loads = load_loads(loads_path)
    elements = build_elements(stations, ref_vec)
    n_nodes = len(stations)

    K = assemble_stiffness(elements, n_nodes)
    F = assemble_loads(elements, loads, n_nodes)

    root_dofs = np.arange(6)  # node 0 clamped: all 6 DOF at the wing root
    u = solve(K, F, root_dofs)

    # Reaction at the clamped root: internal force minus applied load.
    reaction = K[root_dofs] @ u - F[root_dofs]
    nodal = u.reshape(-1, 6)
    return Solution(
        displacements=nodal[:, :3],
        rotations=nodal[:, 3:],
        root_force=reaction[:3],
        root_moment=reaction[3:],
    )
