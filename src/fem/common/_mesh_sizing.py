"""Gmsh size fields for CAD geometry and solution-driven refinement."""

import numpy as np

from .validation import positive_real


def curve_size_fields(gmsh, tolerance, scale, maximum, *, surfaces=False):
    """Limit physical chord sagitta using sampled CAD curvature in model units."""
    if tolerance is None:
        return []
    sagitta = positive_real(tolerance, "geometry_tolerance") * scale
    fields = gmsh.model.mesh.field
    result = []
    entities = gmsh.model.getEntities(1)
    if surfaces:
        entities += gmsh.model.getEntities(2)
    for dim, tag in entities:
        if gmsh.model.getType(dim, tag) in ("Line", "Plane"):
            continue
        lo, hi = gmsh.model.getParametrizationBounds(dim, tag)
        if dim == 1:
            parameters = np.linspace(lo[0], hi[0], 65)
        else:
            # Interior parameter samples avoid the singular poles of spheres.
            axes = [np.linspace(a, b, 19)[1:-1] for a, b in zip(lo, hi)]
            parameters = np.stack(np.meshgrid(*axes), axis=-1).ravel()
        curvature = np.max(np.abs(gmsh.model.getCurvature(dim, tag, parameters)))
        if curvature <= 0:
            continue
        radius = 1 / curvature
        # Exact chord of an osculating circle for the requested sagitta.
        depth = min(sagitta, radius)
        length = 2 * np.sqrt(depth * (2 * radius - depth))
        if length < maximum * scale:
            field = fields.add("Constant")
            fields.setNumbers(field, "CurvesList" if dim == 1 else "SurfacesList", [tag])
            fields.setNumber(field, "VIn", length)
            fields.setNumber(field, "VOut", maximum * scale)
            result.append(field)
    return result


def nodal_size_targets(points, cells, marked, maximum, pairs=()):
    vertices = points[cells]
    n = cells.shape[1]
    diameter = np.maximum.reduce([np.linalg.norm(vertices[:, i] - vertices[:, j], axis=1)
                                  for i in range(n) for j in range(i)])
    diameter[marked] *= 0.5
    target = np.full(len(points), maximum, dtype=float)
    np.minimum.at(target, cells.ravel(), np.repeat(diameter, n))
    if len(pairs):
        slave, master = np.asarray(pairs).T
        target[master] = target[slave] = np.minimum(target[master], target[slave])
    return vertices, target[cells]


def refined_size_targets(mesh, factor, pairs=()):
    """Scale the accepted mesh's local spacings, including both periodic traces."""
    vertices, sizes = nodal_size_targets(mesh.nodes, mesh.elements, slice(0, 0),
                                         mesh.info.requested_maximum_edge, pairs)
    return vertices, sizes / factor


def size_view(gmsh, sizes, origin, scale):
    coords, values = sizes
    dim = coords.shape[-1]
    xyz = np.zeros((len(coords), 3, dim + 1))
    xyz[:, :dim] = ((coords - origin) * scale).transpose(0, 2, 1)
    data = np.column_stack((xyz.reshape(len(coords), -1), values * scale))
    view = gmsh.view.add("adaptive FEM element sizes")
    try:
        gmsh.view.addListData(view, "ST" if dim == 2 else "SS", len(coords), data.ravel())
    except Exception:
        gmsh.view.remove(view)
        raise
    return view


def append_size_view(gmsh, view, fields):
    field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(field, "ViewTag", view)
    fields.append(field)
