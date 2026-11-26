"""Plane-wave phasors and bistatic radar cross sections."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fem.common.constants import C_0
from fem.common.contracts import MeshSnapshot, ResultMixin
from fem.common.errors import ConfigurationError
from fem.common.validation import real_array


@dataclass
class ScatteringResult(ResultMixin):
    """2D nodal Hz/Ez and cell-centred Et/Ht, or 3D cell-centred vector E."""
    frequency: float
    polarization: str
    direction: np.ndarray
    _mesh_snapshot: MeshSnapshot
    incident: np.ndarray
    scattered: np.ndarray
    _positions: np.ndarray
    _current: np.ndarray
    _moment: np.ndarray
    metadata: dict
    current_facets: np.ndarray | None = None
    surface_current_vertices: np.ndarray | None = None
    incident_transverse: np.ndarray | None = None
    scattered_transverse: np.ndarray | None = None

    @property
    def total_transverse(self):
        """Complex (x, y) vectors at physical triangle centres, Et for TE or Ht for TM."""
        if self.scattered_transverse is not None:
            return self.incident_transverse + self.scattered_transverse

    @property
    def surface_current(self):
        """Mean surface-current vector on each conducting facet, in A/m."""
        if self.surface_current_vertices is not None:
            return self.surface_current_vertices.mean(axis=1)

    @property
    def total(self):
        return self.incident + self.scattered

    @property
    def dimension(self):
        return len(self.direction)

    def radar_cross_section(self, directions):
        """Bistatic width (m) or area (m²) for each Cartesian observation direction."""
        directions = np.atleast_2d(real_array(directions, "directions"))
        if directions.ndim != 2 or directions.shape[1] != self.dimension:
            raise ConfigurationError(f"Expected an array of {self.dimension}D directions.")
        norms = np.linalg.norm(directions, axis=1, keepdims=True)
        if np.any(norms == 0):
            raise ConfigurationError("Observation directions must be nonzero.")
        k = 2 * np.pi * self.frequency / C_0
        sigma = np.empty(len(directions))
        for i, n in enumerate(directions / norms):
            phase = np.exp(1j * k * np.einsum("ij,j->i", self._positions, n))
            current = np.einsum("i,ij->j", phase, self._current)
            moment = np.einsum("i,ij->j", phase, self._moment)
            if self.dimension == 2:
                amplitude = current[0] - 1j * k * (n @ moment)
                sigma[i] = abs(amplitude)**2 / (4 * k)
            else:
                amplitude = current - n * (n @ current) + 1j * k * np.cross(n, moment)
                sigma[i] = np.vdot(amplitude, amplitude).real / (4 * np.pi)
        return sigma

    def vtk_datasets(self):
        """Return volume field vectors, a radial RCS surface and any conductor currents."""
        import pyvista as pv

        if self.dimension != 3:
            raise ConfigurationError("VTK datasets require a 3D result.")
        from fem.common.viewer import _tetra_grid
        from .visualization import rcs_surface

        grid = _tetra_grid(self.mesh_data)
        for name in ("incident", "scattered", "total"):
            value = getattr(self, name)
            grid.cell_data[f"{name}_real"] = value.real
            grid.cell_data[f"{name}_imag"] = value.imag
            grid.cell_data[f"{name}_magnitude"] = np.linalg.norm(value, axis=1)
        grid.set_active_vectors("total_real", preference="cell")
        blocks = pv.MultiBlock({"fields": grid, "rcs": rcs_surface(self)})
        if self.current_facets is not None:
            points = self.mesh_data.coordinates[self.current_facets].reshape(-1, 3)
            triangles = np.arange(len(points)).reshape(-1, 3)
            skin = pv.PolyData.from_regular_faces(points, triangles)
            current = self.surface_current_vertices.reshape(-1, 3)
            skin["Js_real"], skin["Js_imag"] = current.real, current.imag
            skin.set_active_vectors("Js_real")
            blocks["surface_current"] = skin
        return blocks

    def export_vtk(self, path):
        """Write a .vtm collection with its field and angular RCS blocks for ParaView."""
        path = Path(path)
        if path.suffix.lower() != ".vtm":
            raise ConfigurationError("Use a .vtm path for the VTK multiblock dataset.")
        path.parent.mkdir(parents=True, exist_ok=True)
        self.vtk_datasets().save(path)
        return path

    def plot(self, *, component=None, quantity="real"):
        """Plot a field such as 'Total Et' (TE) or 'Incident Ez' (TM), or 'RCS'."""
        from .visualization import plot_result
        return plot_result(self, component=component, quantity=quantity)
