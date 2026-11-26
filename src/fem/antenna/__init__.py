"""Waveguide-fed 3D antennas and radiation patterns."""

from .solver import AntennaSolver3D
from .results import AntennaResult, AntennaSweep, RadiationPattern

__all__ = ["AntennaSolver3D", "AntennaResult", "AntennaSweep", "RadiationPattern"]
