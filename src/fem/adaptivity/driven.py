"""Physical-domain interface flux estimator for driven Maxwell and Helmholtz fields."""

import numpy as np
from skfem import InteriorFacetBasis

from .policy import cell_geometry


def driven_residual(basis, coeff, primary_scale, derivative_scale, k, quadrature_order,
                    incident=None, incident_derivative=None):
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
    indicators = np.maximum(indicators, 0.) / max(energy, np.finfo(float).tiny)
    return indicators, float(np.sqrt(indicators.sum()))
