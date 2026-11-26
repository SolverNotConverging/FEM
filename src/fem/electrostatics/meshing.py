"""Gmsh discretization imported directly into scikit-fem meshes."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from threading import Lock
from types import MappingProxyType
from uuid import uuid4

import numpy as np
from fem.common._mesh_sizing import append_size_view, curve_size_fields, size_view
from numpy.typing import NDArray

from .exceptions import MeshError
from .geometry import Circle, GeometryModel, Interval, Polygon, Rectangle


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
_GMSH_LOCK = Lock()


@dataclass(frozen=True, slots=True)
class MeshInfo:
    nodes: int
    elements: int
    minimum_edge: float
    maximum_edge: float
    requested_maximum_edge: float
    material_aware: bool
    interface_refinement: float | None
    boundary_refinement: float | None
    material_element_sizes: dict[str, float]
    element_order: int = 1


@dataclass(frozen=True, slots=True)
class FEMMesh:
    mesh: object
    nodes: FloatArray
    elements: IntArray
    element_tags: NDArray[np.int32]
    physical_names: dict[int, str]
    info: MeshInfo
    geometry_revision: int
    has_curves: bool = False
    charge_density: np.ndarray | None = None
    potential_cells: dict[int, np.ndarray] | None = None

    @property
    def boundary_facets(self):
        """Surfaces of embedded equipotential volumes for the 3D viewer."""
        if self.nodes.shape[1] != 3 or not self.potential_cells:
            return None
        first, second = self.mesh.f2t
        faces = np.zeros(len(first), dtype=bool)
        for cells in self.potential_cells.values():
            other = (second >= 0) & cells[np.maximum(second, 0)]
            faces |= cells[first] != other
        return {"pec": np.flatnonzero(faces)}


def _positive(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise MeshError(f"{name} must be finite and positive.")
    numeric = float(value)
    if not np.isfinite(numeric) or numeric <= 0.0:
        raise MeshError(f"{name} must be finite and positive.")
    return numeric


def _factor(value: float | None, name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)):
        raise MeshError(f"{name} must be in (0, 1] or None.")
    numeric = float(value)
    if not np.isfinite(numeric) or not 0.0 < numeric <= 1.0:
        raise MeshError(f"{name} must be in (0, 1] or None.")
    return numeric


def _import_gmsh() -> object:
    try:
        import gmsh
    except Exception as exc:  # pragma: no cover - depends on host installation
        raise MeshError(
            "Gmsh could not be imported; install the package dependencies with "
            "`pip install -e fem.electrostatics[test]`."
        ) from exc
    return gmsh


def _mesh_extrema(points: FloatArray, elements: IntArray) -> tuple[float, float]:
    lengths = np.stack([np.linalg.norm(points[elements[:, i]] - points[elements[:, j]], axis=1)
                        for i, j in combinations(range(elements.shape[1]), 2)])
    return float(lengths.min()), float(lengths.max())


def _material_targets(
    geometry: GeometryModel,
    maximum: float,
    material_aware: bool,
) -> tuple[dict[int, float], dict[int, str]]:
    table = geometry.material_table
    names = {1: "background", **{index: region.name for index, region in enumerate(geometry.materials, start=2)}}
    if not material_aware:
        return {tag: maximum for tag in table}, names
    scales = {tag: material.dk_scale for tag, material in table.items()}
    reference = max(min(scales.values(), default=1.0), np.finfo(float).tiny)
    return (
        {
            tag: min(maximum, maximum * reference / max(scale, np.finfo(float).tiny))
            for tag, scale in scales.items()
        },
        names,
    )


def _gmsh_session(gmsh: object) -> tuple[bool, str, dict[str, float]]:
    owned = not bool(gmsh.isInitialized())
    if owned:
        gmsh.initialize()
    option_names = (
        "General.Terminal",
        "Mesh.ElementOrder",
        "Mesh.MeshSizeMax",
        "Mesh.MeshSizeMin",
        "Mesh.MeshSizeExtendFromBoundary",
        "Mesh.MeshSizeFromPoints",
        "Mesh.MeshSizeFromCurvature",
    )
    previous = {} if owned else {name: float(gmsh.option.getNumber(name)) for name in option_names}
    return owned, "" if owned else str(gmsh.model.getCurrent()), previous


def _finish_gmsh(
    gmsh: object,
    owned: bool,
    previous_model: str,
    previous_options: dict[str, float],
    model_name: str,
    model_added: bool,
) -> None:
    if owned:
        gmsh.finalize()
        return
    if model_added:
        try:
            gmsh.model.setCurrent(model_name)
            gmsh.model.remove()
            if previous_model:
                gmsh.model.setCurrent(previous_model)
        except Exception:
            pass
    for name, value in previous_options.items():
        try:
            gmsh.option.setNumber(name, value)
        except Exception:
            pass


def discretize_1d(
    geometry: GeometryModel,
    *,
    max_element_size: float,
    material_aware: bool = True,
    interface_refinement: float | None = 0.7,
    boundary_refinement: float | None = 0.5,
) -> FEMMesh:
    """Generate a conforming Gmsh line mesh with high-Dk local sizing."""

    if geometry.dim != 1:
        raise MeshError("discretize_1d requires 1D geometry.")
    maximum = _positive(max_element_size, "max_element_size")
    interface_factor = _factor(interface_refinement, "interface_refinement")
    boundary_factor = _factor(boundary_refinement, "boundary_refinement")
    if not isinstance(material_aware, (bool, np.bool_)):
        raise MeshError("material_aware must be a boolean.")

    cuts = {geometry.x_span[0], geometry.x_span[1]}
    for shape in geometry.all_area_shapes():
        assert isinstance(shape, Interval)
        cuts.update(shape.x)
    ordered = sorted(cuts)
    centers = np.asarray([[(left + right) * 0.5] for left, right in zip(ordered[:-1], ordered[1:], strict=True)])
    segment_tags = geometry.material_indices_at(centers)
    targets, names = _material_targets(geometry, maximum, bool(material_aware))

    constrained_points = {geometry.x_span[0], geometry.x_span[1]}
    for potential in geometry.potentials:
        if isinstance(potential.shape, Interval):
            constrained_points.update(potential.shape.x)
        elif potential.shape == "left":
            constrained_points.add(geometry.x_span[0])
        elif potential.shape == "right":
            constrained_points.add(geometry.x_span[1])

    point_sizes: list[float] = []
    for index, coordinate in enumerate(ordered):
        adjacent = []
        if index:
            adjacent.append(targets[int(segment_tags[index - 1])])
        if index < len(segment_tags):
            adjacent.append(targets[int(segment_tags[index])])
        target = min(adjacent)
        if interface_factor is not None and index not in (0, len(ordered) - 1):
            left_tag = int(segment_tags[index - 1])
            right_tag = int(segment_tags[index])
            if left_tag != right_tag:
                target *= interface_factor
        if boundary_factor is not None and coordinate in constrained_points:
            target *= boundary_factor
        point_sizes.append(target)

    gmsh = _import_gmsh()
    xmin, xmax = geometry.x_span
    scale = 1.0 / (xmax - xmin)
    model_name = f"electrostatic_1d_{uuid4().hex}"
    with _GMSH_LOCK:
        owned, previous_model, previous_options = _gmsh_session(gmsh)
        model_added = False
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add(model_name)
            model_added = True
            occ = gmsh.model.occ
            point_tags = [
                occ.addPoint((coordinate - xmin) * scale, 0.0, 0.0, point_sizes[index] * scale)
                for index, coordinate in enumerate(ordered)
            ]
            line_tags = [
                occ.addLine(point_tags[index], point_tags[index + 1])
                for index in range(len(point_tags) - 1)
            ]
            occ.synchronize()
            for tag in sorted(set(int(value) for value in segment_tags)):
                entities = [line_tags[index] for index, value in enumerate(segment_tags) if int(value) == tag]
                gmsh.model.addPhysicalGroup(1, entities, tag)
                gmsh.model.setPhysicalName(1, tag, names[tag])
            gmsh.option.setNumber("Mesh.ElementOrder", 1)
            gmsh.option.setNumber("Mesh.MeshSizeMax", maximum * scale)
            gmsh.option.setNumber("Mesh.MeshSizeMin", max(min(point_sizes) * scale * 0.1, 1e-12))
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.model.mesh.generate(1)
            node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
            points = np.asarray(coordinates, dtype=float).reshape(-1, 3)[:, :1] / scale + xmin
            element_tags, connectivity = gmsh.model.mesh.getElementsByType(1)
            if not len(element_tags):
                raise MeshError("Gmsh generated no first-order line elements.")
            lookup = {int(tag): index for index, tag in enumerate(node_tags)}
            raw = np.asarray(connectivity, dtype=np.int64).reshape(-1, 2)
            edges = np.fromiter(
                (lookup[int(tag)] for tag in raw.ravel()), dtype=np.int64, count=raw.size
            ).reshape(-1, 2)
        except MeshError:
            raise
        except Exception as exc:
            raise MeshError(f"Gmsh failed to generate the 1D electrostatic mesh: {exc}") from exc
        finally:
            _finish_gmsh(gmsh, owned, previous_model, previous_options, model_name, model_added)

    order = np.argsort(points[:, 0])
    inverse = np.empty_like(order)
    inverse[order] = np.arange(order.size)
    points = points[order]
    edges = inverse[edges]
    centers = points[edges].mean(axis=1)
    tags = geometry.material_indices_at(centers)
    try:
        from skfem import MeshLine
    except Exception as exc:  # pragma: no cover - depends on host installation
        raise MeshError("scikit-fem could not be imported; install the package dependencies.") from exc
    mesh = MeshLine(points.T, edges.T)
    hmin, hmax = _mesh_extrema(points, edges)
    return FEMMesh(
        mesh=mesh,
        nodes=points,
        elements=edges,
        element_tags=tags,
        physical_names=names,
        info=MeshInfo(
            nodes=len(points),
            elements=len(edges),
            minimum_edge=hmin,
            maximum_edge=hmax,
            requested_maximum_edge=maximum,
            material_aware=bool(material_aware),
            interface_refinement=interface_factor,
            boundary_refinement=boundary_factor,
            material_element_sizes=MappingProxyType({names[tag]: value for tag, value in targets.items()}),  # type: ignore[arg-type]
        ),
        geometry_revision=geometry.revision,
    )


def _add_occ_shape(gmsh: object, shape: object, origin: np.ndarray, scale: float) -> tuple[int, int]:
    from fem.common.shapes import Shape
    from fem.common._occ import add_shape
    if isinstance(shape, Shape):
        return add_shape(gmsh, shape, origin, scale)
    occ = gmsh.model.occ
    x0, y0 = origin
    if isinstance(shape, Circle):
        return 2, occ.addDisk(
            (shape.center[0] - x0) * scale,
            (shape.center[1] - y0) * scale,
            0.0,
            shape.radius * scale,
            shape.radius * scale,
        )
    if isinstance(shape, Polygon):
        points = [occ.addPoint((x - x0) * scale, (y - y0) * scale, 0.0) for x, y in shape.points]
        lines = [occ.addLine(points[index], points[(index + 1) % len(points)]) for index in range(len(points))]
        loop = occ.addCurveLoop(lines)
        return 2, occ.addPlaneSurface([loop])
    raise MeshError(f"unsupported curved Gmsh shape {type(shape).__name__}.")


def _add_distance_field(
    gmsh: object,
    entities: list[int],
    target: float,
    maximum: float,
    transition: float,
    scale: float,
    fields: list[int],
    facet_dim: int,
) -> None:
    if not entities:
        return
    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "CurvesList" if facet_dim == 1 else "SurfacesList", entities)
    threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMin", target * scale)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMax", maximum * scale)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", transition * scale)
    fields.append(threshold)


def _fragment_regions(gmsh, geometry, spans, origin, scale):
    """Track CAD descendants for materials, electrodes and space charge."""
    from fem.common._occ import add_shape
    from fem.common.shapes import Box
    occ = gmsh.model.occ
    dim = geometry.dim
    shapes = list(dict.fromkeys(geometry.all_area_shapes()))
    if dim == 2:
        # Preserve the rectangular grid lines used by the planar geometry.
        cuts = [set(span) for span in spans]
        for shape in shapes:
            if isinstance(shape, Rectangle):
                cuts[0].update(shape.x)
                cuts[1].update(shape.y)
        xs, ys = map(sorted, cuts)
        domain = [(2, occ.addRectangle((x-origin[0])*scale, (y-origin[1])*scale, 0.,
                                        (x1-x)*scale, (y1-y)*scale))
                  for x, x1 in zip(xs, xs[1:]) for y, y1 in zip(ys, ys[1:])]
        tools = [shape for shape in shapes if not isinstance(shape, Rectangle)]
    else:
        domain = [add_shape(gmsh, Box(bounds=spans), origin, scale)]
        tools = shapes
    inputs = domain + [_add_occ_shape(gmsh, shape, origin, scale) for shape in tools]
    descendants = [[inputs[0]]]
    if len(inputs) > 1:
        _, descendants = occ.fragment(inputs[:1], inputs[1:])
        if len(descendants) != len(inputs):
            raise MeshError("Gmsh returned incomplete OCC fragment provenance.")
    occ.synchronize()
    regions = {tag for parts in descendants[:len(domain)] for d, tag in parts if d == dim}
    if not regions or regions != {tag for _, tag in gmsh.model.getEntities(dim)}:
        raise MeshError("The fragmented CAD regions do not match the electrostatic domain.")
    membership = {shape: {tag for d, tag in parts if d == dim and tag in regions}
                  for shape, parts in zip(tools, descendants[len(domain):], strict=True)}
    if dim == 2:
        centers = {tag: np.asarray(occ.getCenterOfMass(dim, tag)[:dim])/scale + origin for tag in regions}
        for shape in shapes:
            if isinstance(shape, Rectangle):
                membership[shape] = {tag for tag, pt in centers.items() if shape.contains(*pt)}
    return sorted(regions), membership


def discretize_simplices(
    geometry: GeometryModel,
    *,
    adaptive_sizes=None,
    geometry_tolerance: float | None = None,
    max_element_size: float,
    material_aware: bool = True,
    interface_refinement: float | None = 0.7,
    boundary_refinement: float | None = 0.5,
    interface_refinement_width: float | None = None,
    boundary_refinement_width: float | None = None,
) -> FEMMesh:
    """Generate a tagged triangular or tetrahedral mesh with local size fields."""

    dim = geometry.dim
    if dim not in (2, 3):
        raise MeshError("A triangular or tetrahedral domain is required.")
    maximum = _positive(max_element_size, "max_element_size")
    interface_factor = _factor(interface_refinement, "interface_refinement")
    boundary_factor = _factor(boundary_refinement, "boundary_refinement")
    if not isinstance(material_aware, (bool, np.bool_)):
        raise MeshError("material_aware must be a boolean.")
    interface_width = None if interface_refinement_width is None else _positive(interface_refinement_width, "interface_refinement_width")
    boundary_width = None if boundary_refinement_width is None else _positive(boundary_refinement_width, "boundary_refinement_width")
    if interface_width is not None and interface_factor is None:
        raise MeshError("interface_refinement_width requires interface_refinement.")
    if boundary_width is not None and boundary_factor is None:
        raise MeshError("boundary_refinement_width requires boundary_refinement.")

    targets, names = _material_targets(geometry, maximum, bool(material_aware))
    spans = (geometry.x_span, geometry.y_span, geometry.z_span)[:dim]
    origin = np.array([lo for lo, _ in spans])
    scale = 1.0 / max(hi - lo for lo, hi in spans)
    gmsh = _import_gmsh()
    model_name = f"electrostatic_{dim}d_{uuid4().hex}"
    element_type = 2 if dim == 2 else 4

    with _GMSH_LOCK:
        owned, previous_model, previous_options = _gmsh_session(gmsh)
        model_added = False
        view = None
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add(model_name)
            model_added = True

            regions, membership = _fragment_regions(gmsh, geometry, spans, origin, scale)
            has_curves = any(gmsh.model.getType(dim - 1, tag) != ("Line" if dim == 2 else "Plane")
                             for _, tag in gmsh.model.getEntities(dim - 1))
            region_tags = dict.fromkeys(regions, 1)
            region_potentials = {tag: [] for tag in regions}
            region_charge = dict.fromkeys(regions, 0.)
            for index, material in enumerate(geometry.materials, start=2):
                for tag in membership[material.shape]:
                    region_tags[tag] = index
            for potential in geometry.potentials:
                if not isinstance(potential.shape, str):
                    for tag in membership[potential.shape]:
                        region_potentials[tag].append(potential.id)
            for charge in geometry.charges:
                for tag in membership[charge.shape]:
                    region_charge[tag] = charge.density

            grouped: dict[int, list[int]] = {}
            for entity, tag in region_tags.items():
                grouped.setdefault(tag, []).append(entity)
            for tag, entities in grouped.items():
                gmsh.model.addPhysicalGroup(dim, entities, tag)
                gmsh.model.setPhysicalName(dim, tag, names[tag])

            fields: list[int] = []
            target_values = [maximum]
            if material_aware:
                for tag, entities in grouped.items():
                    target = targets[tag]
                    if target >= maximum * (1.0 - 32.0 * np.finfo(float).eps):
                        continue
                    field = gmsh.model.mesh.field.add("Constant")
                    gmsh.model.mesh.field.setNumber(field, "VIn", target * scale)
                    gmsh.model.mesh.field.setNumber(field, "VOut", maximum * scale)
                    gmsh.model.mesh.field.setNumbers(field, "SurfacesList" if dim == 2 else "VolumesList", entities)
                    fields.append(field)
                    target_values.append(target)

            interface_facets: dict[float, list[int]] = {}
            boundary_facets: dict[float, list[int]] = {}
            region_set = set(regions)
            for _, facet in gmsh.model.getEntities(dim - 1):
                upward, _ = gmsh.model.getAdjacencies(dim - 1, facet)
                adjacent = [int(entity) for entity in upward if int(entity) in region_set]
                adjacent_tags = {region_tags[entity] for entity in adjacent}
                if interface_factor is not None and len(adjacent_tags) > 1:
                    target = interface_factor * min(targets[tag] for tag in adjacent_tags)
                    interface_facets.setdefault(target, []).append(int(facet))
                is_outer = len(adjacent) == 1
                is_potential_edge = len(adjacent) == 2 and region_potentials[adjacent[0]] != region_potentials[adjacent[1]]
                if boundary_factor is not None and (is_outer or is_potential_edge):
                    target = boundary_factor * min(targets[region_tags[entity]] for entity in adjacent)
                    boundary_facets.setdefault(target, []).append(int(facet))

            for target, facets in interface_facets.items():
                _add_distance_field(
                    gmsh,
                    facets,
                    target,
                    maximum,
                    interface_width or 2.0 * target,
                    scale,
                    fields, dim - 1,
                )
                target_values.append(target)
            for target, facets in boundary_facets.items():
                _add_distance_field(
                    gmsh,
                    facets,
                    target,
                    maximum,
                    boundary_width or 3.0 * target,
                    scale,
                    fields, dim - 1,
                )
                target_values.append(target)
            if adaptive_sizes is not None:
                target_values.append(float(np.min(adaptive_sizes[1])))
                view = size_view(gmsh, adaptive_sizes, origin, scale)
                append_size_view(gmsh, view, fields)
            fields.extend(curve_size_fields(gmsh, geometry_tolerance, scale, maximum, surfaces=dim == 3))
            if len(fields) == 1:
                gmsh.model.mesh.field.setAsBackgroundMesh(fields[0])
            elif fields:
                combined = gmsh.model.mesh.field.add("Min")
                gmsh.model.mesh.field.setNumbers(combined, "FieldsList", fields)
                gmsh.model.mesh.field.setAsBackgroundMesh(combined)

            gmsh.option.setNumber("Mesh.ElementOrder", 1)
            gmsh.option.setNumber("Mesh.MeshSizeMax", maximum * scale)
            gmsh.option.setNumber("Mesh.MeshSizeMin", 1e-12 if geometry_tolerance is not None else max(min(target_values) * scale * 0.05, 1e-12))
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
            gmsh.model.mesh.generate(dim)

            node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
            points = np.asarray(coordinates, dtype=float).reshape(-1, 3)[:, :dim] / scale
            points += origin
            element_tags, connectivity = gmsh.model.mesh.getElementsByType(element_type)
            if not len(element_tags):
                raise MeshError("Gmsh generated no first-order simplex elements.")
            lookup = {int(tag): index for index, tag in enumerate(node_tags)}
            raw = np.asarray(connectivity, dtype=np.int64).reshape(-1, dim + 1)
            simplices = np.fromiter(
                (lookup[int(tag)] for tag in raw.ravel()), dtype=np.int64, count=raw.size
            ).reshape(-1, dim + 1)
            element_index = {int(tag): i for i, tag in enumerate(element_tags)}
            tags = np.empty(len(simplices), dtype=np.int32)
            charge_density = np.zeros(len(tags))
            potential_cells = {p.id: np.zeros(len(tags), dtype=bool)
                               for p in geometry.potentials if not isinstance(p.shape, str)}
            for entity in regions:
                cell_tags, _ = gmsh.model.mesh.getElementsByType(element_type, entity)
                indices = [element_index[int(tag)] for tag in cell_tags]
                tags[indices] = region_tags[entity]
                charge_density[indices] = region_charge[entity]
                for owner in region_potentials[entity]:
                    potential_cells[owner][indices] = True
        except MeshError:
            raise
        except Exception as exc:
            raise MeshError(f"Gmsh failed to generate the {dim}D electrostatic mesh: {exc}") from exc
        finally:
            if view is not None:
                gmsh.view.remove(view)
            _finish_gmsh(gmsh, owned, previous_model, previous_options, model_name, model_added)

    try:
        from skfem import MeshTet, MeshTri
    except Exception as exc:  # pragma: no cover - depends on host installation
        raise MeshError("scikit-fem could not be imported; install the package dependencies.") from exc
    mesh = (MeshTri if dim == 2 else MeshTet)(points.T, simplices.T)
    hmin, hmax = _mesh_extrema(points, simplices)
    return FEMMesh(
        mesh=mesh,
        nodes=points,
        elements=mesh.t.T,
        element_tags=tags,
        physical_names=names,
        info=MeshInfo(
            nodes=len(points),
            elements=len(simplices),
            minimum_edge=hmin,
            maximum_edge=hmax,
            requested_maximum_edge=maximum,
            material_aware=bool(material_aware),
            interface_refinement=interface_factor,
            boundary_refinement=boundary_factor,
            material_element_sizes=MappingProxyType({names[tag]: value for tag, value in targets.items()}),  # type: ignore[arg-type]
        ),
        potential_cells=potential_cells, has_curves=has_curves, charge_density=charge_density,
        geometry_revision=geometry.revision,
    )


# Existing callers use this name for planar meshes.
discretize_2d = discretize_simplices

__all__ = ["FEMMesh", "MeshInfo", "discretize_1d", "discretize_2d", "discretize_simplices"]
