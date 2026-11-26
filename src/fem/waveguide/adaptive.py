"""Discretization residuals and adaptive entry points for waveguide modes."""

from dataclasses import replace
import numpy as np
from skfem import MeshLine

from fem.adaptivity import bulk_mark, maxwell_jump_residual, mixed_mode_fields, run_adaptive, validate_controls
from .exceptions import ConfigurationError, SolverError


def estimate_1d(solver, result):
    nodes = solver.mesh_data.nodes
    widths = solver.k0 * np.diff(nodes)
    centers = (nodes[:-1] + nodes[1:]) / 2
    quadratic = solver.mesh_data.info.element_order == 2
    points, weights = np.polynomial.legendre.leggauss((solver._quadrature_order + 2) // 2 if quadratic else 1)
    weights = weights / 2
    dx = widths[:, None] * points / 2
    x = centers[:, None] + dx / solver.k0
    step = np.diff(nodes)[:, None] * 1e-5
    seam_step = np.minimum(np.diff(nodes)[:-1], np.diff(nodes)[1:]) * 1e-8
    material = solver.geometry.transformed_material_at

    def coefficients(polarization, coordinates):
        eps, mu = material(coordinates)
        if polarization == "TE":
            return eps[1], 1 / mu[0], 1 / mu[2]
        return mu[1], 1 / eps[0], 1 / eps[2]

    combined = np.zeros(len(widths))
    maximum = 0.0
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
        sl = coefficients(mode.polarization, nodes[1:-1] - seam_step)[2]
        sr = coefficients(mode.polarization, nodes[1:-1] + seam_step)[2]
        right_slope = derivative + widths * curvature / 2
        left_slope = derivative - widths * curvature / 2
        jumps = abs(sl * right_slope[:-1] - sr * left_slope[1:])**2
        jumps *= active[:-1] & active[1:]
        squared[:-1] += widths[:-1] * jumps / 2
        squared[1:] += widths[1:] * jumps / 2
        scale = float(np.sum(widths * ((abs(stiffness * slope)**2 + abs(b * uq)**2) @ weights) * active))
        normalized = squared / max(scale, np.finfo(float).tiny)
        combined = np.maximum(combined, normalized)
        maximum = max(maximum, float(np.sqrt(normalized.sum())))
    return combined, maximum


def solve_1d(solver, args, options, max_refinements, adaptive_tolerance):
    validate_controls(max_refinements, adaptive_tolerance, ConfigurationError)
    if not solver.discretized:
        solver.mesh(**getattr(solver, "_mesh_settings", {"resolution": 24, "wavelength_elements": 4}))

    def refine(indicators):
        old = solver.mesh_data
        marked = np.arange(len(old.nodes) - 1) if indicators is None else bulk_mark(indicators)
        nodes = np.sort(np.concatenate((old.nodes, (old.nodes[marked] + old.nodes[marked + 1]) / 2)))
        info = replace(old.info, nodes=len(nodes), elements=len(nodes) - 1,
                       minimum_edge=float(np.diff(nodes).min()), maximum_edge=float(np.diff(nodes).max()))
        solver.mesh_data = replace(old, mesh=MeshLine(nodes), nodes=nodes, info=info)
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
                        error_type=ConfigurationError, retry_errors=(SolverError,))


def solve_2d(solver, args, options, max_refinements, adaptive_tolerance):
    validate_controls(max_refinements, adaptive_tolerance, ConfigurationError)
    if not solver.discretized:
        solver.mesh(**getattr(solver, "_mesh_settings", {"wavelength_elements": 4}))

    def estimate(result):
        system = solver.system
        combined = np.zeros(system.physical_mesh.nelements)
        maximum = 0.0
        # All exterior conditions (PEC, PMC, SIBC) are already imposed by
        # the weak/essential forms. Only interior continuity is estimated.
        material = lambda x, y: system.material_at(x / solver.k0, y / solver.k0)
        for index, mode in enumerate(result):
            full = solver.coefficients[:, index]
            fields = lambda basis, coefficient: mixed_mode_fields(basis, full, mode.neff, coefficient)
            indicator, residual = maxwell_jump_residual(system.basis, fields, material, axes=(0, 1),
                                                       quadrature_order=system.quadrature_order)
            combined = np.maximum(combined, indicator / max(indicator.sum(), np.finfo(float).tiny) * residual**2)
            maximum = max(maximum, residual)
        return combined, maximum

    return run_adaptive(solver, lambda: solver._solve_once(*args, **options), estimate,
                        lambda indicators: solver.refine(1.5),
                        max_refinements=max_refinements, adaptive_tolerance=adaptive_tolerance,
                        error_type=ConfigurationError, retry_errors=(SolverError,))
