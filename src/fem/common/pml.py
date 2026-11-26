"""Polynomial Cartesian coordinate stretches for time-harmonic solvers."""
from dataclasses import dataclass

import numpy as np

from .errors import ConfigurationError
from .validation import finite_real, positive_integer, positive_real


@dataclass(frozen=True)
class PMLSpec:
    thickness: float
    order: int = 3
    sigma_max: float = 5.0
    direction: str = "all"

    def __post_init__(self):
        positive_real(self.thickness, "PML thickness")
        positive_integer(self.order, "PML order")
        if finite_real(self.sigma_max, "PML sigma_max") < 0:
            raise ConfigurationError("PML sigma_max must be nonnegative.")
        self.sides("xyz")

    def sides(self, axes):
        if self.direction == "all":
            return tuple(axis + sign for axis in axes for sign in "-+")
        for axis in axes:
            if self.direction == axis:
                return axis + "-", axis + "+"
            if self.direction in (axis + "-", axis + "+"):
                return (self.direction,)
        raise ConfigurationError(f"Invalid PML direction {self.direction!r} for axes {axes}.")

    def stretch(self, depth):
        fraction = np.clip(np.asarray(depth, dtype=float) / self.thickness, 0, 1)
        return 1.0 - 1j * self.sigma_max * fraction**self.order
