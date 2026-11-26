"""Hierarchical first-family Nédélec elements on affine tetrahedra."""

import numpy as np
from skfem import ElementTetN1, FacetBasis
from skfem.element import ElementHcurl
from skfem.refdom import RefTet, RefTri

from .constants import ETA_0


class ElementTetN2(ElementHcurl):
    """20-dimensional P-minus-2 space: two functions per edge and per face.

    Local vertices follow global vertex order. Thus adjacent tetrahedra use
    the same two tangential traces on a face, even for opposite Jacobian signs.
    """

    refdom = RefTet
    edge_dofs = facet_dofs = 2
    maxdeg = 2
    dofnames = ["u^t", "u^t", "u^t", "u^t"]
    doflocs = np.repeat([RefTet.p[:, entity].mean(axis=1)
                         for entity in (*RefTet.edges, *RefTet.facets)], 2, axis=0)
    _gradients = np.vstack((-np.ones(3), np.eye(3)))

    def orient(self, mapping, i, tind=None):
        if i == 0 and np.any(np.diff(mapping.mesh.t, axis=0) < 0):
            raise ValueError("ElementTetN2 requires tetrahedra sorted by global vertex index.")
        count = mapping.mesh.nelements if tind is None else mapping.mesh.t[:, tind].shape[1]
        return np.ones(count)

    def lbasis(self, X, i):
        bary = (1 - X.sum(axis=0), *X)
        grad = self._gradients.reshape(4, 3, *([1] * (X.ndim - 1)))
        if i < 12:
            a, b = RefTet.edges[i // 2]
            weight = 1 if i % 2 == 0 else bary[a] - bary[b]
            dweight = np.zeros_like(grad[a]) if i % 2 == 0 else grad[a] - grad[b]
        else:
            # Of the three cyclic face bubbles lambda_c W_ab, two are
            # independent because their cyclic sum vanishes.
            a, b, c = RefTet.facets[(i - 12) // 2]
            if i % 2:
                a, b, c = b, c, a
            weight, dweight = 4 * bary[c], 4 * grad[c]
        whitney = bary[a] * grad[b] - bary[b] * grad[a]
        curl = 2 * np.cross(grad[a], grad[b], axis=0)
        return weight * whitney, weight * curl + np.cross(dweight, whitney, axis=0)


def tetrahedral_element(order):
    return ElementTetN2() if order == 2 else ElementTetN1()


def conductor_current(basis, coeff, mu, k, facets, incident_curl=None):
    """One-sided n × H at the three corners of each conductor facet (A/m)."""
    wall = FacetBasis(basis.mesh, basis.elem, facets=facets, dofs=basis.dofs,
                      quadrature=(RefTri.p, np.full(3, 1 / 6)))
    curl_e = wall.interpolate(coeff).curl
    if incident_curl is not None:
        curl_e = curl_e + incident_curl(wall.global_coordinates())
    h = 1j * curl_e / (k * ETA_0 * mu[wall.tind, None])
    return np.cross(-wall.normals, h, axis=0).transpose(1, 2, 0)


def surface_vector_samples(vertices, current):
    """Three distinct samples inside each face for its linear vector interpolant."""
    points = (vertices.mean(axis=1, keepdims=True) + vertices) / 2
    phasors = (current.mean(axis=1, keepdims=True) + current) / 2
    return points.reshape(-1, 3), phasors.reshape(-1, 3)
