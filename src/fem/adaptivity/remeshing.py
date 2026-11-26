"""Select the refinement backend and adapt the interiors of conforming CAD regions."""

from contextvars import ContextVar
from itertools import combinations

import numpy as np
from scipy.spatial import cKDTree

from ._budget import RefinementLimit

# Each candidate is built in its own Gmsh model. These settings apply only to
# an adaptive candidate, never to an initial mesh or a periodic rebuild.
_request = ContextVar("fem_mmg_request", default=None)


def adapt_candidate(owner, dim, weights, build, record, *, periodic=False):
    if dim == 1 or periodic:
        record["refinement_backend"] = "native" if dim == 1 else "gmsh"
        return build()
    saved = vars(owner).copy()
    record["refinement_backend"] = f"mmg{dim}d"
    token = _request.set((weights, getattr(owner, "_anisotropic_amr", False)))
    try:
        return build()
    except Exception as error:
        if isinstance(error, RefinementLimit) and error.status != "quality_limit":
            raise
        vars(owner).clear()
        vars(owner).update(saved)
        record["fallback_reason"] = f"MMG{dim}D: {error}"
    finally:
        _request.reset(token)
    record["refinement_backend"] = "gmsh"
    return build()


def _skin(cells):
    dim = cells.shape[1] - 1
    faces = np.concatenate([cells[:, ids] for ids in combinations(range(dim+1), dim)])
    unique, count = np.unique(np.sort(faces, axis=1), axis=0, return_counts=True)
    if np.any(count > 2):
        raise ValueError("Nonmanifold simplex face in candidate mesh.")
    return unique[count == 1]


def _relative_gradients(centers, density):
    """Recover a residual gradient at each old cell and scale it by local density."""
    tree = cKDTree(centers)
    gradient = np.zeros_like(centers)
    dim = centers.shape[1]
    k = min(len(centers), 2 * (dim + 1))
    if k >= dim + 1:
        for start in range(0, len(centers), 4096):
            batch = slice(start, start + 4096)
            _, near = tree.query(centers[batch], k=k)
            patch = centers[near]
            offset = patch - patch.mean(axis=1, keepdims=True)
            values = density[near]
            delta = values - values.mean(axis=1, keepdims=True)
            slope = (np.linalg.pinv(offset) @ delta[..., None])[..., 0]
            peak = np.maximum(values.max(axis=1), np.finfo(float).tiny)
            gradient[batch] = slope / peak[:, None]
    return tree, gradient


def _directional_metric(gradient, lengths):
    """Shorten the target edge by up to a factor of two along the gradient."""
    dim = gradient.shape[1]
    norm = np.linalg.norm(gradient, axis=1)
    direction = gradient / np.where(norm > 0, norm, 1.)[:, None]
    strength = np.minimum(1., lengths * norm)
    return (np.eye(dim) + ((1 + strength)**2 - 1)[:, None, None] *
            direction[:, :, None] * direction[:, None, :]) / lengths[:, None, None]**2


def _adapt_region(mmgpy, points, cells, skin, sizes, recovery):
    dim = points.shape[1]
    # MMG expects positively oriented input simplices.
    determinant = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
    flip = determinant < 0
    cells[flip, :2] = cells[flip, 1::-1]
    lock = np.arange(len(skin), dtype=np.int32)
    if dim == 2:
        mesh = mmgpy.MmgMesh2D()
        mesh.set_mesh_size(vertices=len(points), triangles=len(cells), edges=len(skin))
        mesh.set_vertices(points)
        mesh.set_triangles(cells)
        mesh.set_edges(skin)
        mesh.set_required_edges(lock)
        get_cells = mesh.get_triangles
    else:
        mesh = mmgpy.MmgMesh3D()
        mesh.set_mesh_size(vertices=len(points), tetrahedra=len(cells), triangles=len(skin))
        mesh.set_vertices(points)
        mesh.set_tetrahedra(cells)
        mesh.set_triangles(skin)
        mesh.set_required_triangles(lock)
        get_cells = mesh.get_tetrahedra
    if recovery is None:
        mesh.set_field("metric", sizes[:, None])
    else:
        center_tree, gradients = recovery
        _, nearest = center_tree.query(points)
        metric = _directional_metric(gradients[nearest], sizes)
        i, j = np.triu_indices(dim)
        mesh.set_field("tensor", np.ascontiguousarray(metric[:, i, j]))
    report = mesh.remesh(verbose=-1)
    if report["return_code"] != 0:
        raise RuntimeError(f"MMG returned {report['return_code']}: {report['warnings']}")
    vertices = mesh.get_vertices()
    simplices = get_cells()
    if not len(simplices) or not np.isfinite(vertices).all() or simplices.min() < 0 or simplices.max() >= len(vertices):
        raise ValueError("MMG returned invalid simplex connectivity or coordinates.")
    measure = np.linalg.det(vertices[simplices[:, 1:]] - vertices[simplices[:, :1]])
    if not np.all(measure > 0) or not np.isclose(measure.sum(), abs(determinant).sum(), rtol=1e-8, atol=0):
        raise ValueError("MMG changed region measure or returned degenerate cells.")
    return vertices, simplices


def adapt_cad_interiors(gmsh, dim, sizes, origin, scale):
    """Keep the Gmsh CAD skeleton, replacing each interior with its MMG candidate.

    Both sides of every material, conductor, monitor and port surface retain
    the exact same facet connectivity and Gmsh node tags. Solver-specific
    ownership and DOFs can then be rebuilt by the usual mesh importers.
    """
    request = _request.get()
    if request is None:
        return
    import mmgpy

    weights, anisotropic = request
    positions, targets = sizes
    positions = (positions - origin) * scale
    tree = cKDTree(positions.reshape(-1, dim))
    targets = targets.ravel() * scale
    recovery = None
    if anisotropic and weights is not None:
        if hasattr(weights, "get"):
            weights = weights.get()
        volume = abs(np.linalg.det(positions[:, 1:] - positions[:, :1]))
        density = np.maximum(np.asarray(weights), 0.) / volume
        recovery = _relative_gradients(positions.mean(axis=1), density)

    api = gmsh.model.mesh
    tags, xyz, _ = api.getNodes()
    order = np.argsort(tags)
    tags = tags[order]
    xyz = xyz.reshape(-1, 3)[order, :dim]
    next_tag = int(tags[-1]) + 1
    cell_type = 2 if dim == 2 else 4
    candidates = []
    for _, entity in gmsh.model.getEntities(dim):
        _, raw = api.getElementsByType(cell_type, entity)
        if not len(raw):
            continue
        used, inverse = np.unique(raw, return_inverse=True)
        cells = inverse.reshape(-1, dim+1).astype(np.int32)
        points = np.ascontiguousarray(xyz[np.searchsorted(tags, used)])
        skin = _skin(cells)
        fixed = np.unique(skin)
        owned = api.getNodes(dim, entity)[0]
        if not np.array_equal(used[fixed], np.setdiff1d(used, owned)):
            raise ValueError(f"CAD entity {entity} contains embedded constraints.")
        _, nearest = tree.query(points)
        # Respect finer user, material and curvature sizing in the CAD mesh.
        diameter = np.maximum.reduce([np.linalg.norm(points[cells[:, i]] - points[cells[:, j]], axis=1)
                                      for i, j in combinations(range(dim+1), 2)])
        lengths = targets[nearest].copy()
        np.minimum.at(lengths, cells.ravel(), np.repeat(diameter, dim+1))
        vertices, simplices = _adapt_region(mmgpy, points, cells, skin, lengths, recovery)
        distance, boundary = cKDTree(vertices).query(points[fixed])
        if np.any(distance > 1e-12) or len(np.unique(boundary)) != len(fixed):
            raise ValueError(f"MMG moved a protected vertex of CAD entity {entity}.")
        new_tags = np.zeros(len(vertices), dtype=np.uint64)
        new_tags[boundary] = used[fixed]
        actual = np.sort(new_tags[_skin(simplices)], axis=1)
        expected = used[skin]
        if not np.array_equal(np.unique(actual, axis=0), expected):
            raise ValueError(f"MMG changed protected facets of CAD entity {entity}.")
        interior = new_tags == 0
        count = np.count_nonzero(interior)
        new_tags[interior] = np.arange(next_tag, next_tag+count, dtype=np.uint64)
        next_tag += count
        coords = np.zeros((count, 3))
        coords[:, :dim] = vertices[interior]
        candidates.append((entity, new_tags[interior], coords, new_tags[simplices]))

    for entity, node_tags, coordinates, connectivity in candidates:
        api.clear([(dim, entity)])
        api.addNodes(dim, entity, node_tags, coordinates.ravel())
        api.addElementsByType(entity, cell_type, [], connectivity.ravel())
