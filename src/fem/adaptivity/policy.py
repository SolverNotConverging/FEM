"""Shared solve/estimate/refine policy and Maxwell interface residuals."""


from dataclasses import fields, replace
from math import factorial

import numpy as np
from skfem import InteriorFacetBasis
from ._budget import RefinementLimit, field_dofs, resource_limits
from ._volume import derivative_operators, material_from_volume, modal_volume
from fem.common.discretization import DEFAULT_QUADRATURE_ORDER
from fem.common.progress import report_refinement


def validate_controls(max_refinements, adaptive_tolerance, error_type=ValueError):
    try:
        valid = (not isinstance(max_refinements, (bool, np.bool_, str, bytes))
                 and np.isscalar(max_refinements)
                 and int(max_refinements) == max_refinements and max_refinements >= 0)
    except (TypeError, ValueError, OverflowError):
        valid = False
    if not valid:
        raise error_type("max_refinements must be a nonnegative integer.")
    try:
        valid = (not isinstance(adaptive_tolerance, (bool, np.bool_, str, bytes))
                 and np.isscalar(adaptive_tolerance)
                 and np.isrealobj(adaptive_tolerance)
                 and np.isfinite(adaptive_tolerance) and adaptive_tolerance > 0)
    except (TypeError, ValueError):
        valid = False
    if not valid:
        raise error_type("adaptive_tolerance must be finite and positive.")
    return int(max_refinements), float(adaptive_tolerance)


def run_adaptive(owner, solve_once, estimate, refine, *, max_refinements,
                 adaptive_tolerance, error_type=ValueError, retry_errors=(), max_elements=None, max_dofs=None):
    """Perform one initial solve and at most max_refinements mesh updates.

    An exhausted budget returns the last solution with converged=False.  A
    failed eigensolve may retry on a finer mesh, but is never counted as a
    converged pass.  Failed remeshing restores the last usable owner state.
    """
    limit, tolerance = validate_controls(max_refinements, adaptive_tolerance, error_type)
    limits = resource_limits(max_elements, max_dofs)

    def update_mesh(weights):
        saved = owner.__dict__.copy()
        owner._adaptive_limits = limits
        try:
            refine(weights)
        except Exception:
            owner.__dict__.clear()
            owner.__dict__.update(saved)
            raise
        finally:
            owner.__dict__.pop("_adaptive_limits", None)

    history = []
    for iteration in range(limit + 1):
        try:
            result = solve_once()
        except retry_errors as exc:
            history.append({"refinement": iteration, "residual": None,
                            "status": "solve_failed", "message": str(exc)})
            if iteration == limit:
                raise
            update_mesh(None)
            continue
        indicators, residual, *components = estimate(result)
        if not np.isfinite(residual):
            raise error_type("Adaptive discretization residual is non-finite.")
        converged = residual <= tolerance
        status = "tolerance" if converged else "refinement_limit" if iteration == limit else "refine"
        history.append({"refinement": iteration, "elements": len(indicators),
                        "residual": float(residual), "status": status})
        if components:
            history[-1]["residual_components"] = components[0]
        if len(components) > 1:
            history[-1].update(components[1])
        info = getattr(result, "metadata", None) or getattr(result, "solve_info", {})
        data = getattr(owner, "mesh_data", None)
        if data is not None:
            mesh = getattr(data, "mesh", None)
            if mesh is None:
                mesh = getattr(owner, "_mesh", None)
            if mesh is not None:
                history[-1].update(mesh_quality(mesh))
                order = data.info.element_order if hasattr(data.info, "element_order") else data.info["element_order"]
                mixed = hasattr(getattr(owner, "_system", None), "transverse_indices")
                history[-1]["field_dofs"] = field_dofs(mesh, order, mixed)
                history[-1]["mesh_elements"] = mesh.nelements
        modes = getattr(result, "modes", ())
        pml_fractions = [m.pml_fraction for m in modes if getattr(m, "pml_fraction", None) is not None]
        if pml_fractions:
            history[-1]["max_mode_pml_fraction"] = max(pml_fractions)
        if info.get("ndofs") is not None:
            history[-1]["ndofs"] = info["ndofs"]
        if status == "refine":
            try:
                update_mesh(indicators)
            except RefinementLimit as exc:
                history[-1].update(status=exc.status, message=str(exc))
                status = exc.status
        report_refinement(owner, iteration, limit, len(indicators), residual, tolerance, status,
                          info.get("ndofs"))
        if status != "refine":
            break
    key = "metadata" if any(field.name == "metadata" for field in fields(result)) else "solve_info"
    updated = replace(result, **{key: {
        **getattr(result, key), "adaptive_history": tuple(history),
        "adaptive_residual": float(residual), "adaptive_converged": bool(converged),
        "max_refinements": limit, "adaptive_tolerance": tolerance,
        "max_elements": limits[0], "max_dofs": limits[1],
    }})
    for name, value in tuple(vars(owner).items()):
        if value is result:
            setattr(owner, name, updated)
    return updated


def mesh_quality(mesh):
    """Simplex mean-ratio quality (one for a regular simplex, zero at collapse)."""
    dim = mesh.dim()
    points = mesh.p[:, mesh.t].transpose(2, 1, 0)
    edges = [np.sum((points[:, i] - points[:, j])**2, axis=1)
             for i in range(dim + 1) for j in range(i)]
    squared = np.array(edges)
    determinant = abs(np.linalg.det(points[:, 1:] - points[:, :1]))
    ratio = 2 * (determinant / np.sqrt(dim + 1))**(2 / dim) / squared.mean(axis=0)
    return {"minimum_quality": float(ratio.min()), "quality_p05": float(np.percentile(ratio, 5)),
            "minimum_edge": float(np.sqrt(squared.min())), "maximum_edge": float(np.sqrt(squared.max()))}


def bulk_mark(indicators, fraction=0.5):
    indices = np.argsort(-np.asarray(indicators), kind="stable")
    count = np.searchsorted(np.cumsum(indicators[indices]), fraction * np.sum(indicators)) + 1
    return indices[:count]


def cell_geometry(mesh):
    vertices = mesh.p[:, mesh.t].transpose(2, 1, 0)
    dim = mesh.dim()
    measure = np.abs(np.linalg.det(vertices[:, 1:] - vertices[:, :1])) / factorial(dim)
    h = np.maximum.reduce([np.linalg.norm(vertices[:, i] - vertices[:, j], axis=1)
                           for i in range(dim + 1) for j in range(i)])
    return h, measure


def periodic_facet_pairs(mesh, node_pairs):
    """Pair seam facets by their translated vertices in matching row order."""
    slave, master = np.asarray(node_pairs, dtype=int).T
    boundary = mesh.boundary_facets()
    vertices = mesh.facets[:, boundary]
    representatives = np.arange(mesh.nvertices)
    representatives[slave] = master
    def facet_map(nodes):
        selected = boundary[np.isin(vertices, nodes).all(axis=0)]
        return {tuple(sorted(representatives[mesh.facets[:, f]])): f for f in selected}
    left, right = facet_map(slave), facet_map(master)
    if left.keys() != right.keys():
        raise ValueError("Periodic seam facets do not match by vertex correspondence.")
    return np.array([(f, right[key]) for key, f in left.items()], dtype=int).reshape(-1, 2)


def seam_permutation(left, right):
    """Align quadrature points in tangential coordinates on translated z faces."""
    distance = np.sum((left[:-1, :, :, None] - right[:-1, :, None, :])**2, axis=0)
    return distance.argmin(axis=-1)


def maxwell_jump_residual(basis, fields, material, *, axes, excluded_facets=(),
                          periodic_node_pairs=None, quadrature_order=DEFAULT_QUADRATURE_ORDER):
    """Evaluate interface jumps for a point evaluator or a mesh-owned CellMaterial."""
    if hasattr(material, "on_cells"):
        def coefficients(cells, *coords):
            return material.on_cells(*coords, cells=cells)
    else:
        coefficients = material_from_volume(basis, material)

    indicators, relative, _ = maxwell_residual(
        basis, fields, coefficients, axes=axes, neff=None, excluded_facets=excluded_facets,
        periodic_node_pairs=periodic_node_pairs, quadrature_order=quadrature_order)
    return indicators, relative


def maxwell_residual(basis, fields, cell_material, *, axes, neff, natural=(), excluded_facets=(),
                          periodic_node_pairs=None, quadrature_order=DEFAULT_QUADRATURE_ORDER):
    """Estimate the Maxwell volume equations and interface and natural-boundary traces.

    `fields(basis, material)` returns E and impedance-scaled H, each shaped
    (3, cells/facets, quadrature points). Coordinates are the basis's own
    computational coordinates. Each trace uses its adjacent cell's bulk material and the PML
    stretch evaluated at the trace quadrature point.
    Periodic envelope traces are compared across corresponding seam facets.
    This estimator is not an algebraic residual or a certified error bound.
    """
    mesh = basis.mesh
    h, _ = cell_geometry(mesh)

    coefficients = cell_material(slice(None), *basis.global_coordinates())
    electric, magnetic = fields(basis, lambda *coords: coefficients)
    eps, _ = coefficients
    dfield = eps * electric
    scale = float(np.sum(basis.dx * (np.sum(abs(dfield)**2, axis=0)
                                     + np.sum(abs(magnetic)**2, axis=0))))
    indicators = np.zeros(mesh.nelements)
    facets = np.setdiff1d(np.flatnonzero(mesh.f2t[1] >= 0), excluded_facets)

    def trace(facet_ids, side):
        fb = InteriorFacetBasis(mesh, basis.elem, facets=facet_ids, side=side,
                                dofs=basis.dofs, intorder=quadrature_order)
        coords = np.asarray(fb.global_coordinates())
        coefficients = cell_material(fb.tind, *coords)

        e, hfield = fields(fb, lambda *unused: coefficients)
        epsilon, _ = coefficients
        return fb, epsilon * e, hfield, e

    def add_jump(fb, dd, dh, left, right):
        normal = np.zeros((3, *fb.normals.shape[1:]))
        normal[list(axes)] = fb.normals
        normal_d = np.sum(normal * dd, axis=0)
        tangent_h = dh - normal * np.sum(normal * dh, axis=0)
        contribution = 0.5 * (h[left] + h[right]) * np.sum(
            fb.dx * (abs(normal_d)**2 + np.sum(abs(tangent_h)**2, axis=0)), axis=1)
        np.add.at(indicators, left, contribution / 2)
        np.add.at(indicators, right, contribution / 2)

    for start in range(0, len(facets), 2048):
        selected = facets[start:start + 2048]
        fb, d0, h0, _ = trace(selected, 0)
        _, d1, h1, _ = trace(selected, 1)
        add_jump(fb, d0 - d1, h0 - h1, *mesh.f2t[:, selected])
    if periodic_node_pairs is not None and len(periodic_node_pairs):
        seam = periodic_facet_pairs(mesh, periodic_node_pairs)
        if len(seam):
            fa, da, ha, _ = trace(seam[:, 0], 0)
            fb, db, hb, _ = trace(seam[:, 1], 0)
            order = seam_permutation(fa.global_coordinates(), fb.global_coordinates())[None]
            add_jump(fa, da - np.take_along_axis(db, order, axis=2),
                     ha - np.take_along_axis(hb, order, axis=2), fa.tind, fb.tind)
    if neff is None:
        return indicators, float(np.sqrt(indicators.sum() / max(scale, np.finfo(float).tiny))), {}
    components = {"interface": indicators.sum()}
    pde, divergence = modal_volume(dfield, magnetic, neff, axes, derivative_operators(basis))
    for name, density in (("element_pde", pde), ("divergence", divergence)):
        local = h**2 * np.sum(basis.dx * density, axis=1)
        indicators += local
        components[name] = local.sum()
    components["natural_boundary"] = 0.
    for selected, admittance in natural:
        if not len(selected):
            continue
        fb, _, hfield, e = trace(selected, 0)
        normal = np.zeros_like(e.real)
        normal[list(axes)] = fb.normals
        tangent_e = e - normal * np.sum(normal * e, axis=0)
        mismatch = np.cross(normal, hfield, axis=0) + admittance * tangent_e
        local = h[fb.tind] * np.sum(fb.dx * np.sum(abs(mismatch)**2, axis=0), axis=1)
        np.add.at(indicators, fb.tind, local)
        components["natural_boundary"] += local.sum()
    scale = max(scale, np.finfo(float).tiny)
    components = {name: float(np.sqrt(value / scale)) for name, value in components.items()}
    return indicators, float(np.sqrt(indicators.sum() / scale)), components


def mixed_mode_fields(basis, coefficients, neff, material):
    et, ez = basis.interpolate(coefficients)
    eps, mu = material(*basis.global_coordinates())
    electric = np.stack((et[0], et[1], np.asarray(ez)))
    curl = np.stack((ez.grad[1] + 1j * neff * et[1],
                     -ez.grad[0] - 1j * neff * et[0], et.curl))
    return electric, 1j * curl / mu
