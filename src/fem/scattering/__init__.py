"""Plane-wave scattering by bounded targets in free space."""
from .results import ScatteringResult
from .solver import ScatteringSolver2D, ScatteringSolver3D

__all__ = ["ScatteringResult", "ScatteringSolver2D", "ScatteringSolver3D"]
