"""WingBox: a minimal 3D finite-element solver for cantilevered beams.

The beam is discretised with 1D, 3-node quadratic (Lagrange) Timoshenko
elements carrying 6 DOF per node ``[u, v, w, rx, ry, rz]``.
"""

from wingbox.stations import Section, load_stations
from wingbox.loads import Loads, load_loads
from wingbox.element import BeamElement
from wingbox.assemble import Solution, InternalForces, solve_wing, internal_forces
from wingbox.plot import plot_deformation, plot_internal_forces

__all__ = [
    "Section",
    "load_stations",
    "Loads",
    "load_loads",
    "BeamElement",
    "Solution",
    "InternalForces",
    "solve_wing",
    "internal_forces",
    "plot_deformation",
    "plot_internal_forces",
]
