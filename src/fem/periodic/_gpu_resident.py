"""GPU 3D periodic filters, fields and estimator with hybrid or resident iteration.

Geometry callbacks and scikit-fem trace geometry are prepared once on CPU.
All mode-dependent numerics then run on CUDA; host result objects are exported
only after the independent polynomial/Gauss/direction/PML filters and estimator.
"""
from dataclasses import asdict
from time import perf_counter
import weakref

import numpy as np
from skfem import InteriorFacetBasis

from fem.adaptivity.policy import cell_geometry, periodic_facet_pairs, seam_permutation
from fem.adaptivity._volume import derivative_operators, domain_residuals, modal_volume
from fem.common.gpu import _runtime, quadratic_operator
from fem.common.constants import ETA_0
from fem.common._gpu_resident_arnoldi import iterate_device, check_resident_available
from fem.common._gpu_arnoldi import iterate_hybrid_device, probe_hybrid_available
from .assembly_3d import evaluate_material_3d
from .exceptions import SolverError
from .results import PeriodicMode, PeriodicModeSet, PeriodicSampledFields

class _Basis:
    def __init__(self, basis, solver):
        import cupy as cp
        self.dofs = cp.asarray(basis.element_dofs)
        self.values = cp.asarray(np.stack([np.asarray(item[0]) for item in basis.basis]))
        self.curls = cp.asarray(np.stack([np.asarray(item[0].curl) for item in basis.basis]))
        self.dx = cp.asarray(basis.dx)
        coordinates = np.asarray(basis.global_coordinates()) / solver.k0
        cells = slice(None) if basis.tind is None else basis.tind
        eps, mu = solver.system.material_at.on_cells(*coordinates, cells=cells)
        self.epsilon, self.mu = cp.asarray(eps), cp.asarray(mu)

    def fields(self, full, neff, array_module=None):
        if array_module is None:
            import cupy as cp
        else:
            cp = array_module
        local = full[self.dofs]
        electric = cp.einsum('lceq,le->ceq', self.values, local)
        curl = cp.einsum('lceq,le->ceq', self.curls, local)
        cross = cp.stack((-electric[1], electric[0], cp.zeros_like(electric[0])))
        magnetic = 1j*(curl-1j*neff*cross)/self.mu
        return electric, magnetic  # impedance-scaled H


class ResidentSystem3D:
    def __init__(self, solver, system):
        cp, csr, _ = _runtime()
        self.solver, self.system = weakref.proxy(solver), system
        self.device_index = cp.cuda.Device().id
        self.a0, self.a1, self.a2 = (csr(a) for a in (system.A0, system.A1, system.A2))
        self.prolongation = csr(system.prolongation.matrix)
        self.g0, self.g1 = csr(system.gauss0), csr(system.gauss1)
        self.g0_norm, self.g1_norm = cp.linalg.norm(self.g0.data), cp.linalg.norm(self.g1.data)
        self.volume = _Basis(system.basis, solver)
        coordinates = np.asarray(system.basis.global_coordinates()/solver.k0)
        eps, mu = evaluate_material_3d(system.material_at.bulk, *coordinates)
        self.physical_epsilon, self.physical_mu = cp.asarray(eps), cp.asarray(mu)
        self.pml = cp.asarray(solver.geometry.pml_mask(*coordinates[:2]))
        self.weights = self.volume.dx/solver.k0**3
        self.physical_weights = self.weights*(~self.pml)
        self.weight_sum = cp.sum(self.weights, axis=1)
        weights_cpu = np.asarray(system.basis.dx/solver.k0**3)
        sums = weights_cpu.sum(axis=1)
        self.sample_points = (np.sum(coordinates*weights_cpu[None], axis=2)/sums).T
        self.sample_epsilon = (np.sum(eps*weights_cpu[None], axis=2)/sums).T
        self.sample_mu = (np.sum(mu*weights_cpu[None], axis=2)/sums).T
        self.sample_pml = solver.geometry.pml_mask(*self.sample_points[:, :2].T).astype(float)
        mesh = system.basis.mesh
        h, _ = cell_geometry(mesh)
        self.h = cp.asarray(h)
        self.derivatives = tuple(cp.asarray(a) for a in derivative_operators(system.basis))
        groups = solver.mesh_data.boundary_facets
        natural = np.union1d(groups.get("pmc", ()), groups.get("outer_pmc", ())).astype(int)
        self.natural = None
        if len(natural):
            self.natural = (*self._trace(natural, 0), cp.asarray(mesh.f2t[0, natural]))
        interior = np.flatnonzero(mesh.f2t[1] >= 0)
        self.interior = None
        if len(interior):
            left = self._trace(interior, 0)
            right = self._trace(interior, 1)
            self.interior = (left, right, cp.asarray(mesh.f2t[:, interior]))
        pairs = np.asarray(solver.mesh_data.periodic_node_pairs)
        seam = periodic_facet_pairs(mesh, pairs)
        self.seam = None
        if len(seam):
            left, coords_a = self._trace(seam[:, 0], 0, with_coordinates=True)
            right, coords_b = self._trace(seam[:, 1], 0, with_coordinates=True)
            order = cp.asarray(seam_permutation(coords_a, coords_b)[None])
            self.seam = (left, right, cp.asarray(mesh.f2t[0, seam].T), order)
        cp.cuda.get_current_stream().synchronize()

    def _trace(self, facets, side, with_coordinates=False):
        import cupy as cp
        basis, solver = self.system.basis, self.solver
        mesh = basis.mesh
        fb = InteriorFacetBasis(mesh, basis.elem, facets=facets, side=side,
                                dofs=basis.dofs, intorder=self.system.quadrature_order)
        trace = _Basis(fb, solver), cp.asarray(fb.normals)
        return (trace, fb.global_coordinates()) if with_coordinates else trace

    def physical_data(self, full, neff):
        import cupy as cp
        electric, scaled_h = self.volume.fields(full, neff)
        magnetic = scaled_h/ETA_0
        power = .5*cp.sum(self.physical_weights*(electric[0]*magnetic[1].conj()-
                                                electric[1]*magnetic[0].conj())) / self.solver.period
        density = cp.sum(cp.abs(self.physical_epsilon)*cp.abs(electric)**2, axis=0)+cp.sum(
            cp.abs(self.physical_mu)*cp.abs(scaled_h)**2, axis=0)
        total = cp.sum(self.weights*density).real
        energy = cp.maximum(cp.sum(self.physical_weights*density).real, 0.)
        fraction = cp.clip(cp.sum(self.weights*self.pml*density).real/cp.maximum(total, 1e-300), 0., 1.)
        ee = cp.sum(cp.abs(electric)**2, axis=(1, 2))
        hh = cp.sum(cp.abs(magnetic)**2, axis=(1, 2))
        polarization = cp.where((ee.sum() > 0) & (ee[2]/cp.maximum(ee.sum(), 1e-300) <= 1e-2), 0,
                                cp.where((hh.sum() > 0) & (hh[2]/cp.maximum(hh.sum(), 1e-300) <= 1e-2), 1, 2))
        return power, energy, fraction, polarization, electric, magnetic

    def estimate(self, full, neff, array_module=None):
        if array_module is None:
            import cupy as cp
        else:
            cp = array_module
        electric, magnetic = self.volume.fields(full, neff)
        dfield = self.volume.epsilon*electric
        scale = cp.sum(self.volume.dx*(cp.sum(cp.abs(dfield)**2, axis=0)+
                                       cp.sum(cp.abs(magnetic)**2, axis=0)))
        indicator = cp.zeros(self.system.basis.mesh.nelements)
        if self.interior is not None:
            (left, normal), (right, _), cells = self.interior
            e0, h0 = left.fields(full, neff)
            e1, h1 = right.fields(full, neff)
            jump = h0-h1
            normal_d = cp.sum(normal*(left.epsilon*e0-right.epsilon*e1), axis=0)
            tangent_h = jump-normal*cp.sum(normal*jump, axis=0)
            value = .5*(self.h[cells[0]]+self.h[cells[1]])*cp.sum(left.dx*(
                cp.abs(normal_d)**2+cp.sum(cp.abs(tangent_h)**2, axis=0)), axis=1)
            cp.add.at(indicator, cells[0], value/2)
            cp.add.at(indicator, cells[1], value/2)
        if self.seam is not None:
            (left, normal), (right, _), cells, order = self.seam
            e0, h0 = left.fields(full, neff)
            e1, h1 = right.fields(full, neff)
            dd = left.epsilon*e0 - cp.take_along_axis(right.epsilon*e1, order, axis=2)
            dh = h0 - cp.take_along_axis(h1, order, axis=2)
            nd = cp.sum(normal * dd, axis=0)
            tangent = dh - normal * cp.sum(normal * dh, axis=0)
            value = cp.mean(self.h[cells], axis=0) * cp.sum(
                left.dx * (cp.abs(nd)**2 + cp.sum(cp.abs(tangent)**2, axis=0)), axis=1)
            cp.add.at(indicator, cells[0], value/2)
            cp.add.at(indicator, cells[1], value/2)
        interface = indicator.sum()
        pde, divergence = modal_volume(dfield, magnetic, neff, (0, 1, 2), self.derivatives, cp)
        pde = self.h**2 * cp.sum(self.volume.dx * pde, axis=1)
        divergence = self.h**2 * cp.sum(self.volume.dx * divergence, axis=1)
        indicator += pde + divergence
        wall = cp.asarray(0.)
        if self.natural is not None:
            face, normal, cells = self.natural
            _, hfield = face.fields(full, neff)
            tangent = cp.cross(normal, hfield, axis=0)
            local = self.h[cells] * cp.sum(face.dx * cp.sum(abs(tangent)**2, axis=0), axis=1)
            cp.add.at(indicator, cells, local)
            wall = local.sum()
        scale = cp.maximum(scale, np.finfo(float).tiny)
        components = cp.sqrt(cp.stack((interface, pde.sum(), divergence.sum(), wall)) / scale)
        return indicator, cp.sqrt(indicator.sum() / scale), components

    def solve(self, *, target, requested, direction, eigensolver_tolerance, residual_tolerance,
              divergence_tolerance, propagation_ratio_tolerance, pml_limit, ncv,
              max_iterations, random_seed, prepare_seconds, backend="gpu"):
        import cupy as cp
        system, solver = self.system, self.solver
        count = min(max(4*requested, requested+12), 2*system.ndofs-2)
        shift = target+1e-6j*max(1., abs(target))
        with quadratic_operator(self.a0, self.a1, self.a2, shift) as operator:
            iteration = iterate_device if backend == "gpu" else iterate_hybrid_device
            values, companion, iteration_residuals, execution = iteration(
                operator, target, shift, count, eigensolver_tolerance, ncv, max_iterations, random_seed)
            vectors = companion[:system.ndofs]
            vectors /= cp.linalg.norm(vectors, axis=0)[None]
            qep = operator.residual(companion, values)
        full = self.prolongation @ vectors
        norm = cp.linalg.norm(vectors, axis=0)
        gauss_denominator = ((self.g0_norm+cp.abs(values)*self.g1_norm)*norm)**2
        gauss = cp.linalg.norm(self.g0@vectors+(self.g1@vectors)*values[None], axis=0)**2/cp.where(
            gauss_denominator == 0, 1., gauss_denominator)
        powers = cp.empty(count, dtype=complex)
        energy, fraction = cp.empty(count), cp.empty(count)
        polarization = cp.empty(count, dtype=int)
        for i in range(count):
            powers[i], energy[i], fraction[i], polarization[i], _, _ = self.physical_data(full[:, i], values[i])
        propagating = (cp.abs(powers) > np.finfo(float).tiny) & (
            cp.abs(powers.real)/cp.maximum(cp.abs(powers), 1e-300) >= propagation_ratio_tolerance)
        threshold = 1e-10*cp.maximum(1., cp.abs(values))
        codes = cp.where(propagating, cp.where(powers.real > 0, 1, -1),
                         cp.where(values.imag < -threshold, 2, cp.where(values.imag > threshold, -2, 0)))
        valid = cp.isfinite(values) & cp.isfinite(qep) & cp.isfinite(gauss) & (
            qep <= residual_tolerance) & (gauss <= divergence_tolerance) & (propagating | (energy > 1e-300))
        valid &= cp.isfinite(iteration_residuals) & (iteration_residuals <= eigensolver_tolerance)
        if pml_limit is not None:
            valid &= fraction <= pml_limit
        if direction != 'all':
            valid &= codes > 0 if direction == 'forward' else codes < 0
        accepted = cp.zeros(count, dtype=bool)
        for i in range(count):
            overlap = cp.abs(vectors.conj().T@vectors[:, i])
            close = cp.abs(values-values[i]) <= 1e-8*cp.maximum(1., cp.maximum(cp.abs(values), cp.abs(values[i])))
            duplicate = cp.any(accepted & close & (overlap > 1.-1e-8))
            accepted[i] = valid[i] & ~duplicate & (cp.count_nonzero(accepted) < requested)
        chosen = cp.argsort(cp.where(accepted, cp.arange(count), count))[:requested]
        selected_full, selected_fields, selected_metadata = [], [], []
        combined = cp.zeros(system.basis.mesh.nelements)
        maximum = cp.asarray(0.)
        worst_components = cp.zeros(4)
        worst_domains = cp.zeros(2)
        for j in range(requested):
            i = chosen[j]
            scale = cp.where(propagating[i], 1/cp.sqrt(cp.maximum(cp.abs(powers[i].real), 1e-300)),
                             1/cp.sqrt(cp.maximum(energy[i], 1e-300)))
            coefficient = full[:, i]*scale
            coefficient *= cp.exp(-1j*cp.angle(coefficient[cp.argmax(cp.abs(coefficient))]))
            power, _, pml, pol, electric, magnetic = self.physical_data(coefficient, values[i])
            e = cp.sum(electric*self.weights[None], axis=2)/self.weight_sum
            h = cp.sum(magnetic*self.weights[None], axis=2)/self.weight_sum
            indicator, residual, components = self.estimate(coefficient, values[i])
            combined = cp.maximum(combined, indicator/cp.maximum(indicator.sum(), np.finfo(float).tiny)*residual**2)
            if solver.geometry.pmls:
                domains = cp.stack(domain_residuals(indicator, residual, self.pml.any(axis=1), cp))
                worst_domains = cp.where(residual >= maximum, domains, worst_domains)
            worst_components = cp.where(residual >= maximum, components, worst_components)
            maximum = cp.maximum(maximum, residual)
            selected_full.append(coefficient)
            selected_fields.append(cp.concatenate((e, h), axis=0).T)
            selected_metadata.append(cp.stack((values[i], qep[i], gauss[i], pml, power, codes[i], pol, propagating[i])))
        # Final export boundary: no vectors, field traces or Ritz matrices have
        # been downloaded while iterating, filtering, reconstructing or estimating.
        accepted_count = int(cp.asnumpy(cp.count_nonzero(accepted)))
        if accepted_count < requested:
            raise SolverError(f'Requested {requested} {direction} resident GPU modes; '
                              f'only {accepted_count} passed all physical filters.')
        metadata_host = cp.asnumpy(cp.stack(selected_metadata))
        coefficients_host = cp.asnumpy(cp.stack(selected_full))
        fields_host = cp.asnumpy(cp.stack(selected_fields))
        residual_host = float(cp.asnumpy(maximum))
        modes = []
        for i, row in enumerate(metadata_host):
            neff, residual, gauss_residual, pml, power, code, pol, prop = row
            fields = PeriodicSampledFields(
                self.sample_points, {name: fields_host[i, :, j] for j, name in enumerate(
                    ('Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'))}, dimension=3,
                mesh_points=solver.mesh_data.nodes, mesh_cells=solver.mesh_data.elements,
                sample_element_indices=np.arange(solver.mesh_data.elements.shape[0]),
                material=self.sample_epsilon,
                metadata=dict(sampling='element-barycentre', time_convention='exp(+1j*omega*t - 1j*k0*neff*z)',
                              cell_epsilon_r=self.sample_epsilon, cell_mu_r=self.sample_mu,
                              cell_pml_fraction=self.sample_pml))
            modes.append(PeriodicMode(neff=complex(neff), k0=solver.k0, period=solver.period,
                fields=fields, coefficients=coefficients_host[i], index=i,
                polarization=('TE', 'TM', 'hybrid')[int(pol.real)], power=complex(power),
                direction={1:'forward', -1:'backward', 2:'right-decaying', -2:'left-decaying', 0:'indeterminate'}[int(code.real)],
                normalization='unit-power' if prop.real else 'energy-like', residual=float(residual.real),
                gauss_residual=float(gauss_residual.real), pml_fraction=float(pml.real),
                metadata=dict(classification='propagating' if prop.real else 'evanescent',
                              eigensolver='cuda-resident-refined-arnoldi' if backend == 'gpu' else 'cuda-hybrid-refined-arnoldi')))
        execution.update(device_index=self.device_index, prepare_seconds=prepare_seconds, final_numeric_download_bytes=(
            metadata_host.nbytes+coefficients_host.nbytes+fields_host.nbytes+16))
        result = PeriodicModeSet(modes, frequency=solver.frequency, period=solver.period,
            dimension=3, metadata=dict(backend=backend, eigensolver='cuda-resident-refined-arnoldi' if backend == 'gpu' else 'cuda-hybrid-refined-arnoldi',
                gpu_resident=backend == "gpu", gpu_postprocessing=True, gpu_execution=execution, neff_guess=target, requested_direction=direction,
                mesh=asdict(solver.mesh_data.info), mesh_element_tags=solver.mesh_data.element_tags,
                periodic_node_pairs=solver.mesh_data.periodic_node_pairs,
                periodic_edge_pairs=solver.mesh_data.periodic_edge_pairs, edge_nodes=solver.mesh_data.edge_nodes,
                element_order=solver.mesh_data.info.element_order, cell_edges=solver.mesh_data.cell_edges, cell_edge_signs=solver.mesh_data.cell_edge_signs,
                physical_names=dict(solver.mesh_data.physical_names), boundary_facets={
                    name: np.asarray(solver.mesh_data.mesh.facets[:, facets].T, dtype=np.int64)
                    for name, facets in solver.mesh_data.boundary_facets.items() if len(facets)},
                time_convention='exp(+1j*omega*t - 1j*k0*neff*z)'))
        components = dict(zip(("interface", "element_pde", "divergence", "natural_boundary"),
                              cp.asnumpy(worst_components).tolist(), strict=True))
        domain_report = dict(zip(("physical_residual", "pml_residual"), cp.asnumpy(worst_domains).tolist())) if solver.geometry.pmls else {}
        solver._gpu_resident_estimate = (result, combined, residual_host, components, domain_report)
        solver._result, solver.modes = result, result
        solver.neff, solver.beta, solver.gamma = result.neff, result.beta, result.gamma
        return result


def solve_resident(solver, system, **options):
    import cupy as cp
    (check_resident_available if options.get("backend", "gpu") == "gpu" else probe_hybrid_available)()
    session = getattr(solver, '_gpu_resident_session', None)
    prepared = 0.
    if session is None or session.system is not system or session.device_index != cp.cuda.Device().id:
        start = perf_counter()
        session = ResidentSystem3D(solver, system)
        prepared = perf_counter()-start
        solver._gpu_resident_session = session
    return session.solve(prepare_seconds=prepared, **options)
