"""Closed-surface electromagnetic equivalence for exp(+i omega t) phasors."""

from dataclasses import dataclass

import numpy as np
from skfem import FacetBasis

from fem.common.constants import ETA_0
from fem.common.errors import ConfigurationError
from fem.common.nf2ff import interface_facets
from fem.common.validation import real_array


@dataclass
class HuygensSurface:
    k: float
    coordinates: np.ndarray
    electric_current: np.ndarray
    magnetic_current: np.ndarray

    @classmethod
    def sample(cls, mesh, element, box, coefficients, k, quadrature_order):
        facets, inside_first = interface_facets(mesh, box)
        face = FacetBasis(mesh, element, facets=facets, intorder=quadrature_order)
        other = FacetBasis(mesh, element, facets=facets, side=1, intorder=quadrature_order)
        field = face.interpolate(coefficients)
        normal = face.normals * np.where(inside_first, 1., -1.)[None, :, None]
        # Tangential E is continuous; average the two neighbouring traces of curl E.
        h = 1j * (field.curl + other.interpolate(coefficients).curl) / (2 * k * ETA_0)
        coordinates = np.moveaxis(face.global_coordinates(), 0, -1).reshape(-1, 3)
        electric = np.moveaxis(np.cross(normal, h, axis=0) * face.dx, 0, -1).reshape(-1, 3)
        magnetic = np.moveaxis(-np.cross(normal, field, axis=0) * face.dx, 0, -1).reshape(-1, 3)
        return cls(k, coordinates, electric, magnetic)

    def amplitude(self, directions):
        """Return the complex Cartesian coefficient of exp(-ikr)/r, in volts."""
        n = np.atleast_2d(real_array(directions, "directions"))
        if n.ndim != 2 or n.shape[1] != 3 or np.any(np.linalg.norm(n, axis=1) == 0):
            raise ConfigurationError("Observation directions must be nonzero Cartesian 3-vectors.")
        n = n / np.linalg.norm(n, axis=1, keepdims=True)
        amplitude = np.empty(n.shape, dtype=complex)
        for start in range(0, len(n), 32):
            rays = n[start : start + 32]
            phase = np.exp(1j * self.k * np.einsum("ij,kj->ik", rays, self.coordinates))
            j = np.einsum("ij,jk->ik", phase, self.electric_current)
            m = np.einsum("ij,jk->ik", phase, self.magnetic_current)
            transverse = j - rays * np.sum(j * rays, axis=1)[:, None]
            amplitude[start : start + 32] = (
                -1j * self.k / (4 * np.pi) * (ETA_0 * transverse - np.cross(rays, m))
            )
        return amplitude

    def radiated_power(self):
        # Gauss-Legendre in cos(theta), periodic trapezoidal rule in azimuth.
        radius = 0.5 * np.linalg.norm(np.ptp(self.coordinates, axis=0))
        order = max(16, int(np.ceil(2 * self.k * radius)) + 10)
        mu, weights = np.polynomial.legendre.leggauss(order)
        phi = np.arange(2 * order) * np.pi / order
        z, az = np.meshgrid(mu, phi, indexing="ij")
        radial = np.sqrt(1 - z * z)
        rays = np.stack((radial * np.cos(az), radial * np.sin(az), z), axis=-1)
        intensity = np.sum(abs(self.amplitude(rays.reshape(-1, 3))) ** 2, axis=1) / (2 * ETA_0)
        return float(weights @ intensity.reshape(order, -1).sum(axis=1) * np.pi / order)
