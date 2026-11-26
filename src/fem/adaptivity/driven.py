"""Physical-domain volume, interface and boundary residuals for driven fields."""

import numpy as np
from skfem import FacetBasis, InteriorFacetBasis
from ._volume import derivative_operators, partial

from .policy import cell_geometry


def driven_residual(basis, coeff, primary_scale, derivative_scale, k, quadrature_order,
                    incident=None, incident_derivative=None, *, walls=(), port_traces=(), gap_sources=()):
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
    h, _ = cell_geometry(mesh)
    indicators = np.zeros(count)
    left, right = mesh.f2t
    facets = np.flatnonzero((left < count) & (right >= 0) & (right < count))
    for start in range(0, len(facets), 2048):
        selected = facets[start:start + 2048]
        traces = []
        for side in (0, 1):
            face = InteriorFacetBasis(mesh, basis.elem, facets=selected, side=side,
                                      dofs=basis.dofs, intorder=quadrature_order, disable_doflocs=True)
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
    return indicators, float(np.sqrt(indicators.sum())), {
        name: float(np.sqrt(max(value, 0.) / scale)) for name, value in sums.items()}
