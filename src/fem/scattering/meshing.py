"""Conforming Gmsh meshes with material provenance and PEC cavities."""
from dataclasses import asdict, replace
from threading import RLock
from uuid import uuid4

import numpy as np
from skfem import MeshTet, MeshTri

from fem.common import shapes, materials
from fem.common.materials import bulk_values
from fem.common._occ import add_shape
from fem.common.contracts import ElectromagneticSolverMixin, MeshSnapshot, _context_value, bounds
from fem.common.errors import ConfigurationError, GeometryError
from fem.common.nf2ff import interface_facets
from fem.common.pml import PMLSpec
from fem.common.radiation_pml import extend_box
from fem.common.validation import positive_integer, positive_real
from fem.common.discretization import DEFAULT_QUADRATURE_ORDER
from fem.adaptivity import bulk_mark

_lock = RLock()


class RadiationBoxSolver(ElectromagneticSolverMixin):
    def _configure_pml(self, thickness, layers):
        self._mesh_orders = {"element_order": self._default_element_order, "quadrature_order": DEFAULT_QUADRATURE_ORDER}
        self._adaptive_sizes = None
        self._mesh_auto_materials = True
        self._nf_box = None
        self._nf_distance = None
        self._pml_layers = positive_integer(layers, "pml_layers")
        auto = thickness is None or (isinstance(thickness, str) and thickness == "auto")
        width = self.wavelength / 3 if auto else thickness
        self._default_pml = PMLSpec(width)
        self._explicit_pmls = {}

    @property
    def _pml_axes(self):
        return "xyz"[:len(self._ranges)]

    @property
    def _pml_profiles(self):
        return self._explicit_pmls or dict.fromkeys(self._default_pml.sides(self._pml_axes), self._default_pml)

    def _add_pml_spec(self, spec):
        sides = spec.sides(self._pml_axes)
        if self._explicit_pmls.keys() & set(sides):
            raise GeometryError("A PML already covers one of the requested boundaries.")
        profiles = self._explicit_pmls | dict.fromkeys(sides, spec)
        self._validate_interior(self._physical_ranges(profiles))
        self._explicit_pmls = profiles
        self._mesh = self.mesh_data = self.result = None
        return spec

    def _physical_ranges(self, profiles):
        limits = np.array(self._ranges, dtype=float)
        for side, profile in profiles.items():
            axis = self._pml_axes.index(side[0])
            if side[1] == "-":
                limits[axis, 0] += profile.thickness
            else:
                limits[axis, 1] -= profile.thickness
        if np.any(limits[:, 0] >= limits[:, 1]):
            raise GeometryError("The PML thickness leaves no physical interior. Enlarge the box or reduce thickness.")
        return limits

    def _validate_interior(self, limits):
        for shape, _ in self._regions:
            self._check_inside(shape.bounds, limits, "Geometry intersects the PML. Leave vacuum clearance.")
        for extent in self._enclosed_bounds():
            self._check_inside(extent, limits, "The PML requires vacuum clearance from sources and materials.")
        monitor = self._configured_monitor(limits)
        if monitor is not None:
            self._check_inside(monitor.bounds, limits, "The NF2FF surface must be strictly inside the vacuum PML interface.")
            for extent in self._enclosed_bounds():
                self._check_inside(extent, monitor.bounds, "The NF2FF surface must enclose all sources and non-vacuum material with vacuum clearance.")

    def _configured_monitor(self, limits):
        if self._nf_distance is None:
            return self._nf_box
        shape = shapes.Rectangle if len(limits) == 2 else shapes.Box
        return shape(bounds=limits + [self._nf_distance, -self._nf_distance])

    @staticmethod
    def _check_inside(extents, outer, message):
        inner = np.asarray(extents).reshape(-1, 2)
        outer = np.asarray(outer).reshape(-1, 2)
        if np.any(inner[:, 0] <= outer[:, 0]) or np.any(inner[:, 1] >= outer[:, 1]):
            raise GeometryError(message)

    def _enclosed_bounds(self):
        return [shape.bounds for shape, mat in self._regions if not self._vacuum(mat)]

    @staticmethod
    def _vacuum(mat):
        return (isinstance(mat, materials.Material)
                and tuple(materials.bulk_values(mat, form="scalar")) == (1, 1))

    def _check_geometry_clearance(self, extents, *, vacuum=False):
        limits = self._physical_ranges(self._explicit_pmls) if self._explicit_pmls else self._ranges
        self._check_inside(extents, limits, "Geometry must stay inside the box with vacuum clearance from the PML.")
        if not vacuum and (self._nf_box is not None or self._nf_distance is not None):
            monitor = self._configured_monitor(self._physical_ranges(self._pml_profiles))
            self._check_inside(extents, monitor.bounds, "The NF2FF surface must enclose all sources and non-vacuum material with vacuum clearance.")

    def add_nf2ff(self, distance="auto", *, x_range=None, y_range=None, z_range=None):
        """Place a closed vacuum NF2FF surface by its inset from the PML or by axis ranges."""
        dim = len(self._ranges)
        axes = (x_range, y_range, z_range)
        inset = None
        auto = isinstance(distance, str) and distance == "auto"
        if any(axis is not None for axis in axes):
            if not auto or any(axis is None for axis in axes[:dim]) or (dim == 2 and z_range is not None):
                raise ConfigurationError("Specify all physical axis ranges or one NF2FF distance.")
            limits = np.array([bounds(axis) for axis in axes[:dim]])
        elif auto:
            self._nf_distance = self._nf_box = self._mesh = self.mesh_data = self.result = None
            return
        else:
            inset = positive_real(distance, "NF2FF distance")
            limits = self._physical_ranges(self._pml_profiles) + [inset, -inset]
        if np.any(limits[:, 0] >= limits[:, 1]):
            raise GeometryError("NF2FF distance leaves no integration surface inside the PML.")
        self._check_inside(limits, self._physical_ranges(self._pml_profiles), "The NF2FF surface must be strictly inside the vacuum PML interface.")
        for extents in self._enclosed_bounds():
            self._check_inside(extents, limits, "The NF2FF surface must enclose all sources and non-vacuum material with vacuum clearance.")
        self._nf_box = (shapes.Rectangle if dim == 2 else shapes.Box)(bounds=limits) if inset is None else None
        self._nf_distance = inset
        self._mesh = self.mesh_data = self.result = None

    def _nf2ff_surface(self):
        limits = self._physical_ranges(self._pml_profiles)
        self._validate_interior(limits)
        monitor = self._configured_monitor(limits)
        if monitor is not None:
            return monitor
        dim = len(limits)
        objects = np.asarray(self._enclosed_bounds()).reshape(-1, dim, 2)
        center = limits.mean(axis=1)
        extent = (np.column_stack((objects[:, :, 0].min(axis=0), objects[:, :, 1].max(axis=0)))
                  if len(objects) else np.column_stack((center, center)))
        # Aim for three free-space wavelength/10 elements, keeping vacuum on both sides.
        inset = np.minimum(0.3 * self.wavelength, 0.5 * abs(extent - limits))
        return (shapes.Rectangle if dim == 2 else shapes.Box)(bounds=limits + inset * [1, -1])

    def _attach_pml(self, physical, snapshot, faces, element_order, quadrature_order):
        tags = snapshot.metadata["element_tags"]
        vacuum = np.array([True] + [self._vacuum(mat) for _, mat in self._regions])
        nf_faces, _ = interface_facets(physical, self._surface)
        if not len(nf_faces) or not vacuum[tags[physical.f2t[:, nf_faces]]].all():
            raise GeometryError("Both sides of the NF2FF surface must be in vacuum.")
        if not vacuum[tags[physical.f2t[0, faces[0]]]].all():
            raise GeometryError("The physical interface of the PML must be in vacuum.")
        profiles = self._pml_profiles
        mesh = extend_box(physical, self._physical_ranges(profiles), profiles, self._pml_layers)
        self._faces = {}
        for owner, facets in faces.items():
            if owner:
                cells = physical.f2t[0, facets]
                local = np.argmax(physical.t2f[:, cells] == facets, axis=0)
                self._faces[owner] = mesh.t2f[local, cells]
        boundary = mesh.boundary_facets()
        metals = np.concatenate(list(self._faces.values())) if self._faces else []
        self._pml_wall = np.setdiff1d(boundary, metals)
        self._mesh = mesh
        self._mesh_orders = dict(element_order=element_order, quadrature_order=quadrature_order)
        widths = {side: profile.thickness for side, profile in profiles.items()}
        self.mesh_data = replace(snapshot, info={**snapshot.info, **self._mesh_orders,
                                 "pml_thickness": next(iter(widths.values())) if len(set(widths.values())) == 1 else widths,
                                 "pml_profiles": {side: asdict(profile) for side, profile in profiles.items()},
                                 "pml_layers": self._pml_layers,
                                 "nf2ff_bounds": self._surface.bounds})

    def add_mesh_refinement(self, *, shape, max_element_size):
        """Set a local size in a rectangle, box, circle or sphere in physical coordinates."""
        allowed = (shapes.Rectangle, shapes.Circle) if len(self._ranges) == 2 else (shapes.Box, shapes.Sphere)
        if type(shape) not in allowed:
            raise GeometryError("Use a rectangle or circle in 2D and a box or sphere in 3D.")
        self._refinements.append((shape, positive_real(max_element_size, "max_element_size")))
        self._mesh = self.mesh_data = self.result = None

    _default_element_order = 2

    @staticmethod
    def _wavelength_size(wavelength, order):
        # Quadratic elements resolve propagation with fewer cells per wavelength.
        return wavelength / (6 if order == 2 else 10)

    def _refine_driven(self, indicators):
        snapshot = self.mesh_data
        cells = snapshot.elements
        coords = snapshot.coordinates[cells]
        dim = coords.shape[2]
        h = np.maximum.reduce([np.linalg.norm(coords[:, i] - coords[:, j], axis=1)
                               for i in range(dim + 1) for j in range(i)])
        # Keep the current local resolution and halve the length target on the
        # cells carrying half the integrated squared flux-jump indicator.
        h[bulk_mark(indicators)] *= 0.5
        target = np.full(len(snapshot.coordinates), self._size)
        np.minimum.at(target, cells.ravel(), np.repeat(h, dim + 1))
        self._adaptive_sizes = (coords, target[cells])
        try:
            self.mesh(max_element_size=self._size, **self._mesh_orders)
        finally:
            self._adaptive_sizes = None


def _local_sizes(gmsh, refinements, origin, scale, size, view, material_entities):
    fields = gmsh.model.mesh.field
    ids = []
    for (dim, entity), index in material_entities.items():
        if index <= 1:
            continue
        field = fields.add("Constant")
        fields.setNumber(field, "VIn", size * scale / index)
        fields.setNumber(field, "VOut", size * scale)
        fields.setNumbers(field, "VolumesList" if dim == 3 else "SurfacesList", [entity])
        ids.append(field)
    if view is not None:
        tag = fields.add("PostView")
        fields.setNumber(tag, "ViewTag", view)
        ids.append(tag)
    for shape, local in refinements:
        if isinstance(shape, shapes.Circle):
            tag = fields.add("Ball")
            for axis, center in zip("XYZ", (shape.center - origin) * scale):
                fields.setNumber(tag, axis + "Center", center)
            fields.setNumber(tag, "Radius", shape.radius * scale)
        else:
            tag = fields.add("Box")
            limits = (np.reshape(shape.bounds, (-1, 2)) - origin[:, None]) * scale
            for axis, (lo, hi) in zip("XYZ", limits):
                fields.setNumber(tag, axis + "Min", lo)
                fields.setNumber(tag, axis + "Max", hi)
        fields.setNumber(tag, "VIn", min(local, size) * scale)
        fields.setNumber(tag, "VOut", size * scale)
        fields.setNumber(tag, "Thickness", size * scale)
        ids.append(tag)
    if ids:
        combined = fields.add("Min")
        fields.setNumbers(combined, "FieldsList", ids)
        fields.setAsBackgroundMesh(combined)


def build_mesh(ranges, regions, size, *, interior_box=None, port_gaps=None, refinements=(), adaptive_sizes=None, material_aware=False):
    import gmsh

    dim = len(ranges)
    port_gaps = {} if port_gaps is None else port_gaps
    box = (shapes.Rectangle if dim == 2 else shapes.Box)(bounds=ranges)
    origin = np.array(ranges)[:, 0]
    scale = 1 / max(b - a for a, b in ranges)
    with _lock:
        owned = not gmsh.isInitialized()
        if owned:
            gmsh.initialize(interruptible=False)
        previous = gmsh.model.getCurrent()
        # The quadratic scalar solution can resolve a curved scatterer's field
        # before a 24-segment polygon resolves its backscatter.
        curvature_points = 48 if dim == 2 and material_aware else 24
        options = {"General.Terminal": 0, "Mesh.MeshSizeMax": size * scale,
                   "Mesh.MeshSizeMin": 0, "Mesh.MeshSizeFromCurvature": curvature_points,
                   "Mesh.ElementOrder": 1, "Mesh.RecombineAll": 0}
        saved = {key: gmsh.option.getNumber(key) for key in options}
        gmsh.model.add("scattering-" + uuid4().hex)
        view = None
        try:
            for key, value in options.items():
                gmsh.option.setNumber(key, value)
            domain = add_shape(gmsh, box, origin, scale)
            targets = [add_shape(gmsh, shape, origin, scale) for shape, _ in regions]
            labels_in = list(range(len(targets) + 1))
            if interior_box is not None:
                # A conforming vacuum interface for antenna Huygens integration.
                targets.insert(0, add_shape(gmsh, interior_box, origin, scale))
                labels_in.insert(0, 0)
            gap_entities = {}
            gap_inputs = {}
            for name, gap in port_gaps.items():
                gap_inputs[len(targets) + 1] = name
                targets.append(add_shape(gmsh, gap, origin, scale))
                labels_in.append(0)
            owners = {}
            if targets:
                _, provenance = gmsh.model.occ.fragment([domain], targets)
                for source, (index, fragments) in enumerate(zip(labels_in, provenance, strict=True)):
                    if source in gap_inputs:
                        for entity in fragments:
                            if entity[0] == dim:
                                gap_entities.setdefault(entity, []).append(gap_inputs[source])
                        continue
                    for entity in fragments:
                        if entity[0] == dim:
                            owners[entity] = index
            else:
                owners[domain] = 0
            gmsh.model.occ.synchronize()
            if adaptive_sizes is not None:
                coords, sizes = adaptive_sizes
                xyz = np.zeros((len(coords), 3, dim + 1))
                xyz[:, :dim] = ((coords - origin) * scale).transpose(0, 2, 1)
                data = np.column_stack((xyz.reshape(len(coords), -1), sizes * scale))
                view = gmsh.view.add("adaptive FEM element sizes")
                gmsh.view.addListData(view, "ST" if dim == 2 else "SS", len(coords), data.ravel())
            material_entities = {}
            if material_aware:
                for entity, owner in owners.items():
                    if owner:
                        material = regions[owner - 1][1]
                        if isinstance(material, materials.Material):
                            material_entities[entity] = np.sqrt(abs(np.prod(bulk_values(material, form="scalar"))))
            _local_sizes(gmsh, refinements, origin, scale, size, view, material_entities)
            # Resolve the finite terminal separation locally even when the air box
            # is meshed at a much larger fraction of a wavelength.
            for gap in port_gaps.values():
                limits = (np.asarray(gap.bounds).reshape(3, 2) - origin[:, None]) * scale
                padding = 1e-7
                points = gmsh.model.getEntitiesInBoundingBox(
                    *(limits[:, 0] - padding), *(limits[:, 1] + padding), dim=0)
                gmsh.model.mesh.setSize(points, min(size * scale, np.min(np.diff(limits)) / 2))
            gmsh.model.mesh.generate(dim)
            tags, xyz, _ = gmsh.model.mesh.getNodes()
            order = np.argsort(tags)
            cells, labels, metal_faces = [], [], {}
            gap_cells = {name: [] for name in port_gaps}
            cell_offset = 0
            for (_, entity), owner in owners.items():
                gaps = gap_entities.get((dim, entity), ())
                if owner and not isinstance(regions[owner - 1][1], materials.Material):
                    if gaps:
                        raise GeometryError(f"Gap {gaps[0]!r} must be filled with bulk material or vacuum.")
                    for _, face in gmsh.model.getBoundary([(dim, entity)], oriented=False):
                        _, nodes = gmsh.model.mesh.getElementsByType(1 if dim == 2 else 2, face)
                        for face_nodes in nodes.reshape(-1, dim):
                            metal_faces[tuple(sorted(face_nodes))] = owner
                    continue
                _, nodes = gmsh.model.mesh.getElementsByType(2 if dim == 2 else 4, entity)
                cell = nodes.reshape(-1, dim + 1)
                cells.append(cell)
                labels.append(np.full(len(cell), owner, dtype=int))
                for name in gaps:
                    gap_cells[name].extend(range(cell_offset, cell_offset + len(cell)))
                cell_offset += len(cell)
            used, inverse = np.unique(np.vstack(cells), return_inverse=True)
            xyz = xyz.reshape(-1, 3)[order[np.searchsorted(tags[order], used)], :dim]
            points = xyz / scale + origin
            elements = inverse.reshape(-1, dim + 1)
        finally:
            if view is not None:
                gmsh.view.remove(view)
            gmsh.model.remove()
            if owned:
                gmsh.finalize()
            else:
                if previous:
                    gmsh.model.setCurrent(previous)
                for key, value in saved.items():
                    gmsh.option.setNumber(key, value)
    # Canonical local vertex order fixes the tangential face frame shared by
    # the two quadratic Nédélec functions on each face.
    native = (MeshTri if dim == 2 else MeshTet)(points.T, elements.T, sort_t=True)
    boundary = native.boundary_facets()
    face_nodes = np.sort(used[native.facets[:, boundary]].T, axis=1)
    owners = np.array([metal_faces.get(tuple(nodes), 0) for nodes in face_nodes])
    surfaces = {int(owner): boundary[owners == owner] for owner in np.unique(owners)}
    boundaries = {}
    names = {0: "background"}
    context = []
    for owner, (shape, material) in enumerate(regions, start=1):
        names[owner] = material.name
        context.append({"name": material.name, "shape": _context_value(shape),
                        "material": _context_value(material)})
        if owner in surfaces:
            kind = "pec" if material == materials.PEC else "sibc"
            boundaries.setdefault(kind, []).append(native.facets[:, surfaces[owner]].T)
    metadata = {"element_tags": np.concatenate(labels), "physical_names": names,
                "boundary_facets": {kind: np.vstack(parts) for kind, parts in boundaries.items()},
                "context": {"background": _context_value(materials.vacuum), "regions": context}}
    if port_gaps:
        metadata["port_gap_cells"] = {name: np.asarray(indices, dtype=int) for name, indices in gap_cells.items()}
    snapshot = MeshSnapshot(points, native.t.T, tuple("xyz"[:dim]), {"size": size}, metadata)
    return native, snapshot, surfaces
