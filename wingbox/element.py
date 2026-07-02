"""3-node quadratic Timoshenko beam element (6 DOF per node, 3D).

Formulation
-----------
Each node carries 6 DOF in the *local* frame,

    q_i = [u, v, w, rx, ry, rz]   (3 translations + 3 rotations),

so the element has 18 DOF. All six fields are interpolated with the same
quadratic (second-order) Lagrange shape functions over the three nodes
``[end1, mid, end2]`` located at natural coordinates ``xi = -1, 0, +1``.

The six generalized strains are

    eps_x    = du/dx                 (axial)
    kappa_x  = d(rx)/dx              (torsion)
    kappa_y  = d(ry)/dx              (bending in local x-z plane)
    kappa_z  = d(rz)/dx              (bending in local x-y plane)
    gamma_xy = dv/dx - rz            (transverse shear, local y)
    gamma_xz = dw/dx + ry            (transverse shear, local z)

with the diagonal constitutive law

    D = diag(EA, GJ, EIy, EIz, G*A_sy, G*A_sz).

Shear locking is controlled by *selective reduced integration*: the
axial/torsion/bending terms use full 3-point Gauss, while the two
transverse-shear terms use reduced 2-point Gauss.
"""

from __future__ import annotations

import numpy as np

from wingbox.stations import Section

# Gauss-Legendre rules on the reference interval [-1, 1]: (points, weights).
_GAUSS = {
    2: (
        np.array([-1.0, 1.0]) / np.sqrt(3.0),
        np.array([1.0, 1.0]),
    ),
    3: (
        np.array([-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)]),
        np.array([5.0, 8.0, 5.0]) / 9.0,
    ),
}


def shape_functions(xi: float):
    """Quadratic Lagrange shape functions and their xi-derivatives.

    Nodes ordered ``[end1, mid, end2]`` at ``xi = -1, 0, +1``.

    Returns
    -------
    N : (3,) shape function values.
    dN : (3,) derivatives dN/dxi.
    """
    N = np.array([0.5 * xi * (xi - 1.0), 1.0 - xi * xi, 0.5 * xi * (xi + 1.0)])
    dN = np.array([xi - 0.5, -2.0 * xi, xi + 0.5])
    return N, dN


class BeamElement:
    """A single 3-node quadratic Timoshenko beam element.

    Parameters
    ----------
    node_ids : the three global node indices ``[end1, mid, end2]`` (used
        by the assembler; not required to compute the stiffness).
    node_coords : (3, 3) global coordinates of the three nodes, same order.
    section : :class:`~wingbox.stations.Section` (material + geometry).
    ref_vec : reference vector approximating the local y-direction; the
        local frame is built from it by Gram-Schmidt against the beam
        axis. Default ``(0, 0, 1)`` (global Z). For a beam running along
        global Y this yields local y = global Z, local z = global X.
    """

    N_NODES = 3
    N_DOF = 18

    def __init__(
        self,
        node_ids,
        node_coords,
        section: Section,
        ref_vec=(0.0, 0.0, 1.0),
    ):
        self.node_ids = np.asarray(node_ids, dtype=int)
        self.node_coords = np.asarray(node_coords, dtype=float)
        if self.node_coords.shape != (3, 3):
            raise ValueError("node_coords must have shape (3, 3)")
        self.section = section
        self.ref_vec = np.asarray(ref_vec, dtype=float)

    # -- geometry -----------------------------------------------------------

    @property
    def length(self) -> float:
        """Straight-line length between the two end nodes."""
        return float(np.linalg.norm(self.node_coords[2] - self.node_coords[0]))

    @property
    def rotation_matrix(self) -> np.ndarray:
        """(3, 3) matrix R whose rows are the local axes in global coords.

        Maps a global vector to local: ``v_local = R @ v_global``.
        """
        ex = self.node_coords[2] - self.node_coords[0]
        ex = ex / np.linalg.norm(ex)

        ref = self.ref_vec
        if abs(np.dot(ref, ex)) > 1.0 - 1e-8:
            # ref_vec is (nearly) parallel to the axis: pick another.
            ref = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(ref, ex)) > 1.0 - 1e-8:
                ref = np.array([0.0, 1.0, 0.0])

        ey = ref - np.dot(ref, ex) * ex
        ey = ey / np.linalg.norm(ey)
        ez = np.cross(ex, ey)
        return np.vstack((ex, ey, ez))

    def transformation_matrix(self) -> np.ndarray:
        """(18, 18) block-diagonal transformation T with 6 copies of R.

        ``q_local = T @ q_global`` for the full element DOF vector.
        """
        R = self.rotation_matrix
        T = np.zeros((self.N_DOF, self.N_DOF))
        for b in range(self.N_DOF // 3):  # 6 vector blocks (3 nodes x 2)
            T[3 * b : 3 * b + 3, 3 * b : 3 * b + 3] = R
        return T

    # -- stiffness ----------------------------------------------------------

    def _strain_displacement(self, xi: float, jac: float) -> np.ndarray:
        """(6, 18) strain-displacement matrix B at natural coordinate xi.

        Parameters
        ----------
        xi : natural coordinate in [-1, 1].
        jac : Jacobian dx/dxi = length / 2.
        """
        N, dN = shape_functions(xi)
        dNdx = dN / jac
        B = np.zeros((6, self.N_DOF))
        for i in range(self.N_NODES):
            c = 6 * i
            B[0, c + 0] = dNdx[i]          # eps_x    <- u
            B[1, c + 3] = dNdx[i]          # kappa_x  <- rx
            B[2, c + 4] = dNdx[i]          # kappa_y  <- ry
            B[3, c + 5] = dNdx[i]          # kappa_z  <- rz
            B[4, c + 1] = dNdx[i]          # gamma_xy <- dv/dx
            B[4, c + 5] = -N[i]            # gamma_xy <- -rz
            B[5, c + 2] = dNdx[i]          # gamma_xz <- dw/dx
            B[5, c + 4] = N[i]             # gamma_xz <- +ry
        return B

    def _constitutive(self):
        """Return the (bending-group, shear-group) diagonal constitutive parts."""
        s = self.section
        E, G = s.E, s.G
        # Full 6-term diagonal law, split so each group is integrated
        # with its own Gauss rule (selective reduced integration).
        d_bending = np.diag([E * s.A, G * s.J, E * s.Iy, E * s.Iz, 0.0, 0.0])
        d_shear = np.diag([0.0, 0.0, 0.0, 0.0, G * s.A_sy, G * s.A_sz])
        return d_bending, d_shear

    def _integrate(self, D: np.ndarray, ngp: int) -> np.ndarray:
        """Integrate B^T D B with an ``ngp``-point Gauss rule."""
        jac = self.length / 2.0
        pts, wts = _GAUSS[ngp]
        K = np.zeros((self.N_DOF, self.N_DOF))
        for xi, w in zip(pts, wts):
            B = self._strain_displacement(xi, jac)
            K += (B.T @ D @ B) * w * jac
        return K

    def local_stiffness(self) -> np.ndarray:
        """(18, 18) element stiffness in the local frame.

        Axial/torsion/bending use full (3-point) integration; the two
        transverse-shear terms use reduced (2-point) integration.
        """
        d_bending, d_shear = self._constitutive()
        return self._integrate(d_bending, 3) + self._integrate(d_shear, 2)

    def stiffness(self) -> np.ndarray:
        """(18, 18) element stiffness in the global frame, ``T^T K_local T``."""
        T = self.transformation_matrix()
        return T.T @ self.local_stiffness() @ T

    def global_dofs(self) -> np.ndarray:
        """(18,) global DOF indices for this element, assuming 6 DOF/node."""
        return np.concatenate([6 * n + np.arange(6) for n in self.node_ids])
