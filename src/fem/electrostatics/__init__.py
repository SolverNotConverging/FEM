"""Electrostatic potentials, fields, charges and energies in 1D, 2D and 3D."""

from .exceptions import ElectrostaticsError, NotDiscretizedError
from .results import ElectrostaticResult
from .solver import ElectrostaticSolver

__all__ = ["ElectrostaticResult", "ElectrostaticSolver", "ElectrostaticsError", "NotDiscretizedError"]
