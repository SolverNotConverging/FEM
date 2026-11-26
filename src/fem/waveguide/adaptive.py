"""Discretization residuals and adaptive entry points for waveguide modes."""

from dataclasses import replace
import numpy as np
from skfem import MeshLine

from fem.adaptivity._volume import domain_residuals, material_from_volume
from fem.adaptivity.policy import maxwell_residual
from fem.adaptivity import bulk_mark, mixed_mode_fields, run_adaptive, validate_controls
from fem.common._cell_material import CellMaterial
from fem.common.constants import ETA_0
from fem.common._mesh_sizing import nodal_size_targets
from fem.adaptivity._budget import check_candidate
from .assembly import evaluate_material
from .exceptions import ConfigurationError, SolverError


def estimate_1d(solver, result):
    nodes = solver.mesh_data.nodes
    widths = solver.k0 * np.diff(nodes)
    centers = (nodes[:-1] + nodes[1:]) / 2
    quadratic = solver.mesh_data.info.element_order == 2
    points, weights = np.polynomial.legendre.leggauss(max(2, (solver._quadrature_order + 2) // 2))
    weights = weights / 2
    dx = widths[:, None] * points / 2
    x = centers[:, None] + dx / solver.k0
    step = np.diff(nodes)[:, None] * 1e-5
    tags = solver.mesh_data.element_tags
    material = (solver.geometry.transformed_material_at if tags is None
                else CellMaterial(solver.geometry, tags).on_cells)
    _, pec, pmc = solver._boundary_masks(nodes)
    (pec if solver.geometry.outer_boundary == "pec" else pmc)[[0, -1]] = True

    def coefficients(polarization, coordinates):
        eps, mu = material(coordinates)
        if polarization == "TE":
            return eps[1], 1 / mu[0], 1 / mu[2]
        return mu[1], 1 / eps[0], 1 / eps[2]

    combined = np.zeros(len(widths))
    maximum = 0.0
    worst_components = {}
    domain_report = {}
    for mode in result:
        u = np.asarray(mode.metadata["nodal_primary"])
        primary = np.asarray(mode.metadata["midpoint_primary"])
        active = np.asarray(mode.fields.metadata["active_elements"], dtype=bool)
        curvature = 4 * (u[:-1] + u[1:] - 2 * primary) / widths**2 if quadratic else np.zeros_like(primary)
        derivative = np.diff(u) / widths
        uq = primary[:, None] + derivative[:, None] * dx + curvature[:, None] * dx**2 / 2
        slope = derivative[:, None] + curvature[:, None] * dx
        a, b, stiffness = coefficients(mode.polarization, x)
        # Include the spatial derivative of stretched material coefficients inside a PML.
        change = coefficients(mode.polarization, x + step)[2] - coefficients(mode.polarization, x - step)[2]
        residual = stiffness * curvature[:, None] + change / (2 * solver.k0 * step) * slope + (a - mode.neff**2 * b) * uq
        squared = widths**3 * (abs(residual)**2 @ weights) * active
        volume = squared.sum()
        left_flux = coefficients(mode.polarization, nodes[:-1])[2] * (derivative - widths * curvature / 2)
        right_flux = coefficients(mode.polarization, nodes[1:])[2] * (derivative + widths * curvature / 2)
        jumps = abs(right_flux[:-1] - left_flux[1:])**2 * (active[:-1] & active[1:])
        squared[:-1] += widths[:-1] * jumps / 2
        squared[1:] += widths[1:] * jumps / 2
        interface = squared.sum() - volume
        fixed = pec if mode.polarization == "TE" else pmc
        left_wall = active & ~np.r_[False, active[:-1]] & ~fixed[:-1]
        right_wall = active & ~np.r_[active[1:], False] & ~fixed[1:]
        boundary = widths * (left_wall * abs(left_flux)**2 + right_wall * abs(right_flux)**2)
        squared += boundary
        scale = float(np.sum(widths * ((abs(stiffness * slope)**2 + abs(b * uq)**2) @ weights) * active))
        normalized = squared / max(scale, np.finfo(float).tiny)
        combined = np.maximum(combined, normalized)
        relative = float(np.sqrt(normalized.sum()))
        if relative >= maximum:
            maximum = relative
            sums = {"element_pde": volume, "interface": interface, "natural_boundary": boundary.sum()}
            worst_components = {key: float(np.sqrt(value / max(scale, np.finfo(float).tiny)))
                                for key, value in sums.items()}
            if solver.geometry.pmls:
                physical, absorbing = domain_residuals(normalized, relative, ~solver._non_pml_elements(centers))
                domain_report = {"physical_residual": float(physical), "pml_residual": float(absorbing)}
    return combined, maximum, worst_components, domain_report


def solve_1d(solver, args, options, max_refinements, adaptive_tolerance, max_elements, max_dofs):
    validate_controls(max_refinements, adaptive_tolerance, ConfigurationError)
    if not solver.discretized:
        solver.mesh(**getattr(solver, "_mesh_settings", {"resolution": 24, "wavelength_elements": 4}))

    def refine(indicators):
        old = solver.mesh_data
        marked = np.arange(len(old.nodes) - 1) if indicators is None else bulk_mark(indicators)
        nodes = np.sort(np.concatenate((old.nodes, (old.nodes[marked] + old.nodes[marked + 1]) / 2)))
        info = replace(old.info, nodes=len(nodes), elements=len(nodes) - 1,
                       minimum_edge=float(np.diff(nodes).min()), maximum_edge=float(np.diff(nodes).max()))
        counts = np.ones(old.mesh.nelements, dtype=int)
        counts[marked] += 1
        tags = None if old.element_tags is None else np.repeat(old.element_tags, counts)
        native = MeshLine(nodes)
        check_candidate(solver, native, info.element_order)
        solver.mesh_data = replace(old, mesh=native, nodes=nodes, element_tags=tags, info=info)
        # A later explicit refine() must not rebuild from the original coarse
        # resolution and discard the smaller intervals created by adaptation.
        solver._discretization_settings = {
            **solver._discretization_settings, "resolution": None,
            "max_element_size": info.minimum_edge,
        }
        solver._clear_result_views()

    return run_adaptive(solver, lambda: solver._solve_once(*args, **options),
                        lambda result: estimate_1d(solver, result), refine,
                        max_refinements=max_refinements, adaptive_tolerance=adaptive_tolerance,
                        max_elements=max_elements, max_dofs=max_dofs,
                        error_type=ConfigurationError, retry_errors=(SolverError,))


def solve_2d(solver, args, options, max_refinements, adaptive_tolerance, max_elements, max_dofs):
    validate_controls(max_refinements, adaptive_tolerance, ConfigurationError)
    if not solver.discretized:
        solver.mesh(**getattr(solver, "_mesh_settings", {"wavelength_elements": 4}))

    def estimate(result):
        system = solver.system
        combined = np.zeros(system.physical_mesh.nelements)
        maximum = 0.0
        groups = solver.mesh_data.boundary_facets
        natural = [(groups.get(name, ()), 0.) for name in ("pmc", "outer_pmc")]
        natural.extend((ids, ETA_0 / z) for ids, z in system.impedance_boundaries)
        worst_components = {}
        domain_report = {}
        if solver.geometry.pmls:
            centers = solver.mesh_data.nodes[solver.mesh_data.elements].mean(axis=1)
            pml_cells = ~solver._physical_integration_mask(centers.T)

        if isinstance(system.material_at, CellMaterial):
            def material(cells, x, y):
                return system.material_at.on_cells(x / solver.k0, y / solver.k0, cells=cells)
        else:
            material = material_from_volume(system.basis, lambda x, y: evaluate_material(
                system.material_at, x / solver.k0, y / solver.k0))

        cache = {} if len(result) > 1 else None
        for index, mode in enumerate(result):
            full = solver.coefficients[:, index]
            fields = lambda basis, coefficient: mixed_mode_fields(basis, full, mode.neff, coefficient)
            indicator, residual, components = maxwell_residual(
                system.basis, fields, material, axes=(0, 1), neff=mode.neff,
                natural=natural, quadrature_order=system.quadrature_order, cache=cache)
            combined = np.maximum(combined, indicator / max(indicator.sum(), np.finfo(float).tiny) * residual**2)
            if residual >= maximum:
                maximum, worst_components = residual, components
                if solver.geometry.pmls:
                    physical, absorbing = domain_residuals(indicator, residual, pml_cells)
                    domain_report = {"physical_residual": float(physical), "pml_residual": float(absorbing)}
        return combined, maximum, worst_components, domain_report

    def refine(indicators):
        if indicators is None:
            return solver.refine(1.5)
        mesh = solver.mesh_data
        sizes = nodal_size_targets(mesh.nodes, mesh.elements, bulk_mark(indicators), mesh.info.requested_maximum_edge)
        return solver._discretize(**solver._discretization_settings, adaptive_sizes=sizes)

    return run_adaptive(solver, lambda: solver._solve_once(*args, **options), estimate,
                        refine,
                        max_refinements=max_refinements, adaptive_tolerance=adaptive_tolerance,
                        max_elements=max_elements, max_dofs=max_dofs,
                        error_type=ConfigurationError, retry_errors=(SolverError,))
