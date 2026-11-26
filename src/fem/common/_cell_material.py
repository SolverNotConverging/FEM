"""Constitutive coefficients on an OCC-tagged modal mesh."""

import numpy as np


class CellMaterial:
    def __init__(self, geometry, tags):
        self.geometry = geometry
        self.tags = np.asarray(tags)
        materials = [geometry.background, *(region.material for region in geometry.regions)]
        self.epsilon = np.asarray([m.eps_r for m in materials])
        self.mu = np.asarray([m.mu_r for m in materials])

    def __call__(self, *coordinates):
        """Evaluate a point query in the continuous geometry."""
        return self.geometry.transformed_material_at(*coordinates)

    def bulk(self, *coordinates, cells=slice(None)):
        shape = coordinates[0].shape
        indices = self.tags[cells] - 1
        values = []
        for table in (self.epsilon, self.mu):
            diagonal = table[indices].T.reshape(3, len(indices), *([1] * (len(shape) - 1)))
            values.append(np.broadcast_to(diagonal, (3, *shape)))
        return tuple(values)

    def on_cells(self, *coordinates, cells=slice(None)):
        eps, mu = self.bulk(*coordinates, cells=cells)
        if self.geometry.pmls:
            factors = self.geometry.pml_factors(*coordinates)
            return factors * eps, factors * mu
        return eps, mu
