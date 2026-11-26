"""Shared solve/estimate/refine policy and Maxwell interface residuals."""


from dataclasses import fields, replace
from math import factorial

import numpy as np
from skfem import InteriorFacetBasis
from .remeshing import adapt_candidate
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


def stagnating(owner, history):
    window, reduction = getattr(owner, "_progress_policy", (3, .01))
    if window is None or len(history) <= window:
        return False
    recent = history[-window-1:]
    if any(r["status"] not in ("refine", "refined") for r in recent[:-1]):
        return False
    return all(b["residual"] >= (1 - reduction) * a["residual"] for a, b in zip(recent, recent[1:]))


def run_adaptive(owner, solve_once, estimate, refine, *, max_refinements,
                 adaptive_tolerance, error_type=ValueError, retry_errors=(), max_elements=None, max_dofs=None, absorber=None, periodic=False):
    """Perform one initial solve and at most max_refinements mesh updates.

    An exhausted budget returns the last solution with converged=False.  A
    failed eigensolve may retry on a finer mesh, but is never counted as a
    converged pass.  Failed remeshing restores the last usable owner state.
    """
    limit, tolerance = validate_controls(max_refinements, adaptive_tolerance, error_type)
    limits = resource_limits(max_elements, max_dofs)

    def restore(state):
        owner.__dict__.clear()
        owner.__dict__.update(state)

    def update_mesh(weights, *, shell=False):
        saved = owner.__dict__.copy()
        owner._adaptive_limits = limits
        data = getattr(owner, "mesh_data", None)
        native = getattr(data, "mesh", getattr(owner, "_mesh", None))
        owner._previous_cell_count = getattr(native, "nelements", None)
        try:
            if shell:
                owner._refine_absorber()
                history[-1]["refinement_backend"] = "pml_layers"
            elif native is not None:
                adapt_candidate(owner, native.dim(), weights, lambda: refine(weights), history[-1], periodic=periodic)
            else:
                refine(weights)
        except Exception:
            restore(saved)
            raise
        finally:
            owner.__dict__.pop("_adaptive_limits", None)
            owner.__dict__.pop("_previous_cell_count", None)

    history = []
    last_usable = None
    for iteration in range(limit + 1):
        try:
            result = solve_once()
        except retry_errors as exc:
            history.append({"refinement": iteration, "residual": None,
                            "status": "solve_failed", "message": str(exc)})
            try:
                if iteration == limit:
                    raise
                update_mesh(None)
                continue
            except retry_errors + (RefinementLimit,) as failure:
                if last_usable is None:
                    raise
                saved, result, residual = last_usable
                restore(saved)
                history[-1].update(status=getattr(failure, "status", "solve_failed"),
                                   message=str(failure), restored_last_solution=True)
                converged = False
                break
            except Exception:
                if last_usable is not None:
                    restore(last_usable[0])
                raise
        except Exception:
            if last_usable is not None:
                restore(last_usable[0])
            raise
        try:
            indicators, residual, *components = estimate(result)
            if not np.isfinite(residual):
                raise error_type("Adaptive discretization residual is non-finite.")
        except Exception:
            if last_usable is not None:
                restore(last_usable[0])
            raise
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
        if modes:
            history[-1]["effective_indices"] = np.array([m.neff for m in modes])
        if pml_fractions:
            history[-1]["max_mode_pml_fraction"] = max(pml_fractions)
        if info.get("ndofs") is not None:
            history[-1]["ndofs"] = info["ndofs"]
        if status == "refine" and stagnating(owner, history):
            status = "stagnation"
            history[-1]["status"] = status
            if components and components[0]:
                history[-1]["dominant_component"] = max(components[0], key=components[0].get)
                window, reduction = getattr(owner, "_progress_policy", (3, .01))
                previous = history[-window-1].get("residual_components", {})
                history[-1]["stagnating_components"] = [
                    key for key, value in components[0].items()
                    if value > 0 and key in previous and value >= previous[key] * (1-reduction)**window
                ]
        if absorber is not None:
            record = history[-1]
            record["physical_converged"] = converged
            try:
                passed = absorber.measure(result, record)
                record["pml_converged"] = passed
            except Exception:
                if last_usable is not None:
                    restore(last_usable[0])
                raise
            if converged and not passed:
                converged = False
                status = "pml_refine" if iteration < limit and absorber.remaining else "pml_limit"
                record["status"] = status
            if status == "refine":
                # Changes on a different physical mesh cannot certify the shell.
                absorber.previous = None
        if status in ("refine", "pml_refine"):
            last_usable = owner.__dict__.copy(), result, residual
            try:
                update_mesh(indicators, shell=status == "pml_refine")
                if status == "pml_refine":
                    absorber.remaining -= 1
            except RefinementLimit as exc:
                history[-1].update(status=exc.status, message=str(exc))
                status = exc.status
        report_refinement(owner, iteration, limit, len(indicators), residual, tolerance, status,
                          info.get("ndofs"))
        if status not in ("refine", "pml_refine"):
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


def cached(cache, key, build):
    if cache is None:
        return build()
    if key not in cache:
        cache[key] = build()
    return cache[key]


def maxwell_residual(basis, fields, cell_material, *, axes, neff, natural=(), excluded_facets=(),
                          periodic_node_pairs=None, quadrature_order=DEFAULT_QUADRATURE_ORDER, cache=None):
    """Estimate the Maxwell volume equations and interface and natural-boundary traces.

    `fields(basis, material)` returns E and impedance-scaled H, each shaped
    (3, cells/facets, quadrature points). Coordinates are the basis's own
    computational coordinates. Each trace uses its adjacent cell's bulk material and the PML
    stretch evaluated at the trace quadrature point.
    Periodic envelope traces are compared across corresponding seam facets.
    This estimator is not an algebraic residual or a certified error bound.
    """
    mesh = basis.mesh
    h, _ = cached(cache, "geometry", lambda: cell_geometry(mesh))

    coefficients = cached(cache, "volume_material", lambda: cell_material(slice(None), *basis.global_coordinates()))
    electric, magnetic = fields(basis, lambda *coords: coefficients)
    eps, _ = coefficients
    dfield = eps * electric
    scale = float(np.sum(basis.dx * (np.sum(abs(dfield)**2, axis=0)
                                     + np.sum(abs(magnetic)**2, axis=0))))
    indicators = np.zeros(mesh.nelements)
    facets = cached(cache, "interior", lambda: np.setdiff1d(np.flatnonzero(mesh.f2t[1] >= 0), excluded_facets))

    def trace(facet_ids, side):
        def build():
            fb = InteriorFacetBasis(mesh, basis.elem, facets=facet_ids, side=side,
                                    dofs=basis.dofs, intorder=quadrature_order, disable_doflocs=True)
            return fb, cell_material(fb.tind, *fb.global_coordinates())

        # Preserve batched working memory on large meshes.
        trace_cache = cache if mesh.nfacets <= 4096 else None
        fb, coefficients = cached(trace_cache, (side, np.asarray(facet_ids).tobytes()), build)

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
        seam = cached(cache, "seam", lambda: periodic_facet_pairs(mesh, periodic_node_pairs))
        if len(seam):
            fa, da, ha, _ = trace(seam[:, 0], 0)
            fb, db, hb, _ = trace(seam[:, 1], 0)
            order = cached(cache, "permutation", lambda: seam_permutation(
                fa.global_coordinates(), fb.global_coordinates())[None])
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
