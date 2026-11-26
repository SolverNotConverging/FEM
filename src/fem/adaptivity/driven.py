"""Physical-domain and absorbing-shell residuals for driven fields."""

import numpy as np
from skfem import Basis, FacetBasis, InteriorFacetBasis
from ._volume import derivative_operators, partial

from .policy import cached, cell_geometry


def driven_residual(basis, coeff, primary_scale, derivative_scale, k, quadrature_order,
                    incident=None, incident_derivative=None, *, walls=(), port_traces=(), gap_sources=(), absorber=None, cache=None):
    mesh = basis.mesh
    count = len(primary_scale)
    vector = mesh.dim() == 3

    def fluxes(space, cells, include_incident=False):
        field = space.interpolate(coeff)
        primary = np.asarray(field)
        derivative = field.curl if vector else field.grad
        if include_incident and incident is not None:
            x = space.global_coordinates()
            primary = primary + incident(x)
            derivative = derivative + incident_derivative(x)
        return primary * primary_scale[cells, None], derivative * (derivative_scale[cells, None] / k)

    primary, flux = fluxes(basis, basis.tind)
    density = abs(primary)**2
    if vector:
        density = density.sum(axis=0)
    energy = np.sum(basis.dx * (density + np.sum(abs(flux)**2, axis=0)))
    h, _ = cached(cache, "geometry", lambda: cell_geometry(mesh))
    indicators = np.zeros(count)
    left, right = mesh.f2t
    facets = cached(cache, "interior", lambda: np.flatnonzero((left < count) & (right >= 0) & (right < count)))
    for start in range(0, len(facets), 2048):
        selected = facets[start:start + 2048]
        traces = []
        for side in (0, 1):
            face = cached(cache if len(facets) <= 4096 else None, (start, side),
                          lambda: InteriorFacetBasis(mesh, basis.elem, facets=selected, side=side,
                              dofs=basis.dofs, intorder=quadrature_order, disable_doflocs=True))
            traces.append(fluxes(face, face.tind, include_incident=True))
        delta = traces[0][1] - traces[1][1]
        normal = face.normals
        if vector:
            # Continuity of normal D and tangential H across material interfaces.
            normal_d = np.sum(normal * (traces[0][0] - traces[1][0]), axis=0)
            tangent_h = np.cross(normal, delta, axis=0)
            jump = abs(normal_d)**2 + np.sum(abs(tangent_h)**2, axis=0)
        else:
            jump = abs(np.sum(normal * delta, axis=0))**2
        i, j = mesh.f2t[:, selected]
        weight = (h[i] + h[j]) * np.sum(face.dx * jump, axis=1) / 4
        np.add.at(indicators, i, weight)
        np.add.at(indicators, j, weight)
    sums = {"interface": indicators.sum()}
    operators = derivative_operators(basis)
    if vector:
        pde = -k * primary.copy()
        divergence = np.zeros_like(primary[0])
        for axis in range(3):
            b, c = (axis + 1) % 3, (axis + 2) % 3
            pde[c] += partial(flux[b], axis, operators)
            pde[b] -= partial(flux[c], axis, operators)
            divergence += partial(primary[axis], axis, operators)
        for cells, axis, source in gap_sources:
            pde[axis, cells] -= source
        if incident is not None:
            pde += k * (derivative_scale - primary_scale)[:, None] * incident(basis.global_coordinates())
        pde_density = np.sum(abs(pde)**2, axis=0)
        local = h[:count]**2 * np.sum(basis.dx * abs(divergence)**2, axis=1)
        indicators += local
        sums["divergence"] = local.sum()
    else:
        pde = k * primary.copy()
        for axis in range(2):
            pde += partial(flux[axis], axis, operators)
        if incident is not None:
            pde += k * (primary_scale - derivative_scale)[:, None] * incident(basis.global_coordinates())
        pde_density = abs(pde)**2
    local = h[:count]**2 * np.sum(basis.dx * pde_density, axis=1)
    sums["element_pde"] = local.sum()
    indicators += local
    sums["natural_boundary"] = 0.
    for facets, admittance in walls:
        face = FacetBasis(mesh, basis.elem, facets=facets, intorder=quadrature_order, dofs=basis.dofs)
        value, trace = fluxes(face, face.tind, include_incident=True)
        value = value / primary_scale[face.tind, None]
        normal = face.normals
        if vector:
            tangent = value - normal * np.sum(normal * value, axis=0)
            mismatch = np.cross(normal, trace, axis=0) - 1j * admittance * tangent
            density = np.sum(abs(mismatch)**2, axis=0)
        else:
            density = abs(np.sum(normal * trace, axis=0) + 1j * admittance * value)**2
        local = h[face.tind] * np.sum(face.dx * density, axis=1)
        np.add.at(indicators, face.tind, local)
        sums["natural_boundary"] += local.sum()
    for face, prescribed_flux in port_traces:
        _, trace = fluxes(face, face.tind)
        mismatch = np.cross(face.normals, trace, axis=0) - prescribed_flux
        local = h[face.tind] * np.sum(face.dx * np.sum(abs(mismatch)**2, axis=0), axis=1)
        np.add.at(indicators, face.tind, local)
        sums["port_boundary"] = sums.get("port_boundary", 0.) + local.sum()
    scale = max(energy, np.finfo(float).tiny)
    # Some tetrahedral rules have negative weights. Clamp roundoff after integration.
    indicators = np.maximum(indicators, 0.) / scale
    residual = float(np.sqrt(indicators.sum()))
    components = {
        name: float(np.sqrt(max(value, 0.) / scale)) for name, value in sums.items()}

    if absorber is None:
        return indicators, residual, components
    report = shell_residual(basis, coeff, count, k, absorber, h, scale, quadrature_order)
    report["physical_residual"] = residual
    return indicators, residual, components, report


def shell_residual(basis, coeff, count, k, absorber, h, scale, quadrature_order):
    """Strong equations in stretched vacuum and transmission at the shell entrance."""
    from fem.common.radiation_pml import material_diagonal

    mesh = basis.mesh
    vector = mesh.dim() == 3
    volume = interface = energy = 0.
    dofs = basis.dofs.element_dofs
    for start in range(count, mesh.nelements, 1024):
        cells = np.arange(start, min(start + 1024, mesh.nelements))
        if not np.any(coeff[dofs[:, cells]]):
            continue
        shell = Basis(mesh, basis.elem, elements=cells, dofs=basis.dofs, mapping=basis.mapping,
                      quadrature=basis.quadrature, disable_doflocs=True)
        tensor = material_diagonal(shell.global_coordinates(), *absorber)
        field = shell.interpolate(coeff)
        operators = derivative_operators(shell)
        if vector:
            primary, flux = tensor * field, field.curl / (k * tensor)
            equation = -k * primary.copy()
            divergence = np.zeros_like(primary[0])
            for axis in range(3):
                b, c = (axis + 1) % 3, (axis + 2) % 3
                equation[c] += partial(flux[b], axis, operators)
                equation[b] -= partial(flux[c], axis, operators)
                divergence += partial(primary[axis], axis, operators)
            density = np.sum(abs(equation)**2, axis=0) + abs(divergence)**2
            amplitude = np.sum(abs(primary)**2, axis=0)
        else:
            primary, flux = tensor[2] * field, tensor[:2] * field.grad / k
            equation = k * primary
            for axis in range(2):
                equation += partial(flux[axis], axis, operators)
            density = abs(equation)**2
            amplitude = abs(primary)**2
        volume += (h[cells]**2) @ np.maximum(0., np.sum(shell.dx * density, axis=1))
        energy += np.sum(shell.dx * (amplitude + np.sum(abs(flux)**2, axis=0)))

    # Both materials are vacuum on this interface; the stretch is unity at depth zero.
    left, right = mesh.f2t
    faces = np.flatnonzero((right >= 0) & ((left < count) != (right < count)))
    for start in range(0, len(faces), 1024):
        ids = faces[start:start + 1024]
        if not np.any(coeff[dofs[:, mesh.f2t[:, ids]]]):
            continue
        traces = []
        for side in (0, 1):
            face = InteriorFacetBasis(mesh, basis.elem, facets=ids, side=side, dofs=basis.dofs,
                                      intorder=quadrature_order,
                                      disable_doflocs=True)
            traces.append(face.interpolate(coeff))
        a, b = traces
        if vector:
            jump = (abs(np.sum(face.normals * (a-b), axis=0))**2
                    + np.sum(abs(np.cross(face.normals, (a.curl-b.curl)/k, axis=0))**2, axis=0))
        else:
            jump = abs(np.sum(face.normals * (a.grad-b.grad), axis=0) / k)**2
        interface += np.mean(h[mesh.f2t[:, ids]], axis=0) @ np.maximum(0., np.sum(face.dx * jump, axis=1))
    return {"pml_residual": float(np.sqrt(volume / scale)),
            "pml_interface_residual": float(np.sqrt(interface / scale)),
            "pml_field_norm_ratio": float(np.sqrt(max(0., energy) / scale))}


class AbsorberCheck:
    """Compare the same observable on two shells around an identical physical mesh."""

    def __init__(self, tolerance, observable=None, max_refinements=3):
        from fem.common.validation import positive_integer, positive_real
        from fem.common.errors import ConfigurationError
        self.tolerance = positive_real(tolerance, "pml_tolerance")
        if observable is not None and not callable(observable):
            raise ConfigurationError("pml_observable must be a function of the solved result.")
        self.remaining = positive_integer(max_refinements, "max_pml_refinements", minimum=0)
        self.observable = observable
        self.previous = None

    def measure(self, result, record):
        if self.observable is not None:
            name, value = "custom", self.observable(result)
        elif hasattr(result, "s_parameters"):
            name, value = "S", result.s_parameters
        else:
            if result.dimension == 2:
                angle = 2 * np.pi * np.arange(36) / 36
                directions = np.column_stack((np.cos(angle), np.sin(angle)))
            else:
                z = 1 - (2 * np.arange(36) + 1) / 36
                phi = np.arange(36) * np.pi * (3 - np.sqrt(5))
                r = np.sqrt(1 - z*z)
                directions = np.column_stack((r*np.cos(phi), r*np.sin(phi), z))
            name, value = "RCS", result.radar_cross_section(directions)
        values = np.atleast_1d(np.asarray(value, dtype=complex)).copy()
        if not values.size or not np.isfinite(values).all():
            raise ValueError("PML observable must have nonempty finite values.")
        record.update(pml_layers=result.mesh_data.info["pml_layers"], pml_observable=name,
                      observable_values=values, pml_tolerance=self.tolerance)
        status = "unchecked"
        if all(record[key] == 0 for key in ("pml_field_norm_ratio", "pml_interface_residual", "pml_residual")):
            status = "unexcited"
        elif self.previous is not None:
            if values.shape != self.previous.shape:
                raise ValueError("PML observable shape changed between shell refinements.")
            change = np.linalg.norm(values - self.previous) / max(
                np.linalg.norm(values), np.linalg.norm(self.previous), np.finfo(float).tiny)
            record["pml_observable_change"] = float(change)
            status = "tolerance" if change <= self.tolerance else "unresolved"
        record["pml_status"] = status
        self.previous = values
        return status in ("unexcited", "tolerance")
