"""Residual-controlled remeshing with regenerated periodic constraints."""

import numpy as np
from fem.common._mesh_sizing import nodal_size_targets
from fem.adaptivity._volume import domain_residuals
from fem.adaptivity.policy import maxwell_residual
from fem.adaptivity import bulk_mark, run_adaptive, validate_controls
from .exceptions import ConfigurationError, SolverError


def solve_periodic(solver, dimension, args, options, max_refinements, adaptive_tolerance, max_elements, max_dofs):
    validate_controls(max_refinements, adaptive_tolerance, ConfigurationError)
    ready = solver.discretized if dimension == 2 else solver.system is not None
    if not ready:
        solver.mesh(**getattr(solver, "_mesh_settings", {"wavelength_elements": 4}))

    def estimate(result):
        resident = getattr(solver, '_gpu_resident_estimate', None)
        if resident is not None and resident[0] is result:
            return resident[1:]
        combined = np.zeros(solver.mesh_data.info.elements)
        maximum = 0.0
        worst_components = {}
        domain_report = {}
        if solver.geometry.pmls:
            centers = solver.mesh_data.nodes[solver.mesh_data.elements].mean(axis=1)
            pml_cells = solver.geometry.pml_mask(*centers[:, :dimension-1].T)
        caches = {}
        for mode in result:
            system = solver._systems[mode.polarization] if dimension == 2 else solver.system
            magnetic_primary = dimension == 2 and system.polarization == "TM"

            def estimator_material(cells, *coordinates):
                epsilon, mu = system.material_at.on_cells(
                    *(coordinate / solver.k0 for coordinate in coordinates), cells=cells)
                return (mu, epsilon) if magnetic_primary else (epsilon, mu)

            def fields(basis, coefficient):
                field = basis.interpolate(mode.coefficients)
                eps, mu = coefficient(*basis.global_coordinates())
                if dimension == 3:
                    electric = np.asarray(field)
                    cross = np.stack((-electric[1], electric[0], np.zeros_like(electric[0])))
                    return electric, 1j * (field.curl - 1j * mode.neff * cross) / mu
                primary = np.asarray(field)
                derivative = field.grad
                zero = np.zeros_like(primary)
                longitudinal = derivative[1] - 1j * mode.neff * primary
                # Electromagnetic duality uses (H, -E) for a TM primary.
                return np.stack((zero, primary, zero)), np.stack(
                    (-1j * longitudinal / mu[0], zero, 1j * derivative[0] / mu[2]))

            natural_kind = "pec" if magnetic_primary else "pmc"
            groups = solver.mesh_data.boundary_facets
            natural = [(groups.get(name, ()), 0.) for name in (natural_kind, "outer_" + natural_kind)]
            indicator, residual, components = maxwell_residual(
                system.basis, fields, estimator_material, axes=(0, 2) if dimension == 2 else (0, 1, 2),
                neff=mode.neff, natural=natural, quadrature_order=system.quadrature_order,
                cache=caches.setdefault(id(system), {}) if len(result) > 1 else None,
                periodic_node_pairs=(np.column_stack((solver.mesh_data.slave_nodes, solver.mesh_data.master_nodes))
                                     if dimension == 2 else solver.mesh_data.periodic_node_pairs),
            )
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
        if dimension == 2:
            pairs = np.column_stack((mesh.slave_nodes, mesh.master_nodes))
            settings = solver._discretization_settings
        else:
            pairs = mesh.periodic_node_pairs
            settings = solver._discretize_options
        if hasattr(indicators, "get"):  # Transfer cell weights only when CAD remeshing is needed.
            indicators = indicators.get()
        sizes = nodal_size_targets(mesh.nodes, mesh.elements, bulk_mark(indicators),
                                   mesh.info.requested_maximum_edge, pairs)
        return solver._mesh_impl(**settings, adaptive_sizes=sizes)

    return run_adaptive(solver, lambda: solver._solve_once(*args, **options), estimate,
                        refine,
                        max_refinements=max_refinements, adaptive_tolerance=adaptive_tolerance,
                        max_elements=max_elements, max_dofs=max_dofs,
                        error_type=ConfigurationError, retry_errors=(SolverError,), periodic=True)
