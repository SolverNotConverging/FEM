"""Numerical wave ports on conforming cross-sections of the antenna mesh."""

from dataclasses import dataclass, replace

import numpy as np
from scipy.sparse import csr_matrix
from skfem import FacetBasis, LinearForm, MeshTri, asm
from skfem.helpers import dot

from fem.common import shapes
from fem.common.contracts import mesh_snapshot
from fem.common.errors import ConfigurationError, GeometryError
from fem.waveguide import ModeSolver2D, SampledFields
from fem.waveguide.meshing import FEMMesh2D, MeshInfo, _edge_extrema


@dataclass(frozen=True)
class WaveguidePort:
    name: str
    axis: str
    position: float
    cross_section: object
    num_modes: int
    mode: int
    amplitude: float
    phase_deg: float

    @property
    def bounds(self):
        limits = np.asarray(self.cross_section.bounds).reshape(2, 2).tolist()
        limits.insert("xyz".index(self.axis), (self.position, self.position))
        return np.asarray(limits).ravel()

    @property
    def drive(self):
        return self.amplitude * np.exp(1j * np.deg2rad(self.phase_deg))

    def section(self, basis, candidates, length, tags, names, impedance_edges, quadrature_order):
        mesh = basis.mesh
        axis = "xyz".index(self.axis)
        uv = [i for i in range(3) if i != axis]
        vertices = mesh.p[:, mesh.facets[:, candidates]]
        on_plane = np.all(abs(vertices[axis] - self.position) < length * 1e-9, axis=0)
        inside = self.cross_section.contains(*vertices[uv].mean(axis=1))
        facets = np.sort(candidates[on_plane & inside])
        if not len(facets):
            raise GeometryError(f"No conforming reference face for port {self.name!r}.")
        face = FacetBasis(mesh, basis.elem, facets=facets, intorder=quadrature_order)
        if np.ptp(face.normals[axis]) > 1e-8:
            raise GeometryError(f"Port {self.name!r} includes faces with opposing normals.")
        nodes, local = np.unique(mesh.facets[:, facets], return_inverse=True)
        tri = MeshTri(mesh.p[uv][:, nodes], local.reshape(3, -1))
        hmin, hmax = _edge_extrema(tri.p.T, tri.t.T)
        info = MeshInfo(tri.nvertices, tri.nelements, hmin, hmax, hmax, basis.elem.maxdeg)
        boundary = tri.boundary_facets()
        sides = {}
        for edge in boundary:
            pair = tuple(sorted(nodes[tri.facets[:, edge]]))
            wall = impedance_edges.get(pair)
            if wall is not None:
                sides.setdefault(wall, []).append(edge)
        impedance_facets = [edge for edges in sides.values() for edge in edges]
        data = FEMMesh2D(
            tri,
            tags[mesh.f2t[0, facets]],
            names,
            {
                "pec": np.setdiff1d(boundary, impedance_facets),
                "sibc": np.asarray(impedance_facets, dtype=int),
            },
            info,
            0,
        )
        return PortSection(self, face, uv, data, sides, quadrature_order)


@dataclass
class PortSection:
    port: WaveguidePort
    face: FacetBasis
    uv: list
    mesh_data: FEMMesh2D
    wall_facets: dict
    quadrature_order: int

    def solve_modes(self, frequency, props):
        data = self.mesh_data
        tri = data.mesh
        eps, mu = np.asarray(props, dtype=complex)[data.element_tags].T

        # The modal solver evaluates all material samples on native element quadrature.
        def material_at(x, y):
            return (np.broadcast_to(eps[:, None], x.shape), np.broadcast_to(mu[:, None], y.shape))

        solver = ModeSolver2D(
            frequency=frequency,
            x_range=(tri.p[0].min(), tri.p[0].max()),
            y_range=(tri.p[1].min(), tri.p[1].max()),
        )
        walls = [
            (edges, wall.at_frequency(frequency=frequency))
            for wall, edges in self.wall_facets.items()
        ]
        modes = solver._solve_cross_section(
            data, material_at, self.port.num_modes, np.sqrt(np.max(np.abs(eps * mu))) * 1.001, walls,
            self.quadrature_order,
        )

        def order(m):
            if m.metadata["classification"] == "propagating":
                return (0, -m.neff.real)
            return (1, abs(m.neff.imag))

        ordered = []
        for idx, mode in enumerate(sorted(modes, key=order)):
            fields = mode.fields
            # Set the phase against spatial moments in fixed transverse coordinates.
            # Using an edge coefficient as reference would depend on face numbering.
            xy = np.asarray(fields.coordinates)
            xy = (xy - xy.min(axis=1, keepdims=True)) / np.ptp(xy, axis=1, keepdims=True) * 2 - 1
            e = np.array([fields.values[c] for c in ("Ex", "Ey")])
            weights = solver.system.basis.dx.ravel()
            moments = [
                e @ (weights * w)
                for w in (np.ones_like(xy[0]), *xy, xy[0] * xy[1], xy[0] ** 2, xy[1] ** 2)
            ]
            moments = np.asarray(moments).ravel()
            pivot = moments[np.argmax(abs(moments))]
            if abs(pivot) < 1e-9 * np.sum(weights * np.linalg.norm(e, axis=0)):
                pivot = e.ravel()[np.argmax(abs(e))]
            phase = np.exp(-1j * np.angle(pivot))
            sampled = SampledFields(
                fields.coordinates,
                {k: phase * v for k, v in fields.values.items()},
                mesh_points=fields.mesh_points,
                mesh_cells=fields.mesh_cells,
                material=fields.material,
                metadata=fields.metadata,
            )
            metadata = dict(mode.metadata)
            metadata["coefficient_vector"] = phase * metadata["coefficient_vector"]
            ordered.append(replace(mode, fields=sampled, index=idx, metadata=metadata))
        snapshot = mesh_snapshot(data, axes=tuple("xyz"[i] for i in self.uv))
        result = replace(
            modes, modes=tuple(ordered), metadata={**modes.metadata, "port_mesh": snapshot}
        )
        object.__setattr__(result, "_mesh_snapshot", snapshot)
        return result

    def operators(self, modes):
        """Return magnetic loads, reciprocal amplitude projectors, and the power overlap."""
        face = self.face
        weights = face.dx
        # Each native triangle and its 3D facet use the same symmetric quadrature.
        # Match the point permutation explicitly for either local facet orientation.
        xy = np.stack(modes[0].fields.coordinates).reshape(2, *weights.shape)
        face_xy = face.global_coordinates()[self.uv]
        delta = face_xy[:, :, :, None] - xy[:, :, None, :]
        distances = np.linalg.norm(delta, axis=0)
        index = distances.argmin(axis=-1)
        if distances.min(axis=-1).max() > self.mesh_data.info.maximum_edge * 1e-8:
            raise GeometryError("The modal and volume face quadratures do not coincide.")
        electric, magnetic, loads = [], [], []
        for mode in modes:

            def component(key):
                return np.take_along_axis(
                    mode.fields.values[key].reshape(weights.shape), index, axis=1
                )

            e = np.stack([component("Ex"), component("Ey")])
            h = np.stack([component("Hy"), -component("Hx")])
            electric.append(e)
            magnetic.append(h)
            profile = np.zeros(face.global_coordinates().shape, dtype=complex)
            profile[self.uv] = h

            @LinearForm(dtype=complex)
            def load(test, w):
                return dot(test, profile)

            loads.append(asm(load, face))
        # A_mn = integral (H_m x z) . E_n. Unconjugated reciprocity also works
        # for evanescent and lossy modes. Use the complete overlap for hybrid modes.
        gram = np.einsum("mifq,nifq,fq->mn", magnetic, electric, weights)
        flux = np.einsum("mifq,nifq,fq->mn", np.conj(magnetic), electric, weights) / 2
        load = csr_matrix(np.column_stack(loads))
        try:
            projector = csr_matrix(np.linalg.inv(gram)) @ load.T
        except np.linalg.LinAlgError as exc:
            raise ConfigurationError(f"Singular modal basis at port {self.port.name!r}.") from exc
        profiles = np.zeros((len(magnetic), 3, *weights.shape), dtype=complex)
        profiles[:, self.uv] = magnetic
        return load, projector, flux, profiles


@dataclass(frozen=True)
class DiscretePort:
    name: str
    axis: str
    gap: object
    impedance: float
    amplitude: float
    phase_deg: float
    mode = "gap"

    @property
    def bounds(self):
        if isinstance(self.gap, shapes.Box):
            return self.gap.bounds
        base = self.gap
        rotation, shift = np.eye(3), np.zeros(3)
        while isinstance(base, shapes.Transformed):
            center = np.asarray(base.center)
            shift += rotation @ (center + base.offset - base.matrix @ center)
            rotation = rotation @ base.matrix
            base = base.shape
        if not isinstance(base, shapes.Cylinder):
            raise GeometryError("A discrete port needs a Box or Cylinder spanning the feed gap.")
        axis = "xyz".index(self.axis)
        if not np.allclose(abs(rotation[:, 2]), np.eye(3)[axis], atol=1e-10, rtol=0):
            raise GeometryError("The cylinder's longitudinal axis must match the discrete port axis.")
        midpoint = rotation @ [*base.center, np.mean(base.z_range)] + shift
        half_size = np.full(3, base.radius)
        half_size[axis] = (base.z_range[1] - base.z_range[0]) / 2
        return tuple(np.column_stack((midpoint - half_size, midpoint + half_size)).ravel())

    @property
    def drive(self):
        return self.amplitude * np.exp(1j * np.deg2rad(self.phase_deg))

    def voltage_trace(self, basis, metal_facets, cells):
        """Gap voltage is the area-average of line integrals along the positive axis."""
        mesh = basis.mesh
        axis = "xyz".index(self.axis)
        bounds = np.asarray(self.bounds).reshape(3, 2)
        spans = np.diff(bounds).ravel()
        if not len(cells):
            raise GeometryError(f"Gap {self.name!r} has no volume elements.")
        sub = basis.with_elements(cells)
        # Use the faceted mesh area so a uniform axial field integrates to E*length.
        area = sub.dx.sum() / spans[axis]
        terminals, counts = np.unique(mesh.t2f[:, cells], return_counts=True)
        terminals = terminals[counts == 1]
        vertices = mesh.p[:, mesh.facets[:, terminals]]
        for end in bounds[axis]:
            facets = terminals[np.all(abs(vertices[axis] - end) < spans.max() * 1e-8, axis=0)]
            if not len(facets) or not np.isin(facets, metal_facets).all():
                raise GeometryError(f"Both end faces of gap {self.name!r} must contact conductors.")

        @LinearForm
        def voltage(test, w):
            return test[axis] / area

        return csr_matrix(asm(voltage, sub)[:, None])
