"""Manufactured smooth field and shared-face traces of the quadratic tetrahedron."""

import numpy as np
from skfem import Basis, FacetBasis, InteriorFacetBasis, MeshTet

from fem.antenna.solver import _recover_fields
from fem.common.constants import ETA_0
from fem.common.tet_hcurl import ElementTetN2


def test_quadratic_curl_and_surface_current_on_oppositely_oriented_cells():
    raw = MeshTet.init_tensor(*([np.array([-.3, .2, 1.1])] * 3))
    mesh = MeshTet(raw.p, raw.t, sort_t=True)
    basis = Basis(mesh, ElementTetN2(), intorder=5)

    def exact(x):
        return np.array([-x[0] * x[1], x[0]**2, 0 * x[2]])

    coeff = basis.project(exact)
    np.testing.assert_allclose(basis.interpolate(coeff), exact(basis.global_coordinates()), atol=2e-12)
    facets = mesh.boundary_facets()
    _, magnetic, current, corners = _recover_fields(basis, coeff, np.ones(mesh.nelements), 1., facets)
    expected_h = np.zeros_like(magnetic)
    expected_h[:, 2] = 3j / ETA_0 * mesh.p[0, mesh.t].mean(axis=0)
    np.testing.assert_allclose(magnetic, expected_h, atol=1e-12)
    corner_h = np.zeros_like(corners)
    corner_h[:, :, 2] = 3j / ETA_0 * mesh.p[0, mesh.facets[:, facets]].T
    wall = FacetBasis(mesh, basis.elem, facets=facets)
    expected = np.cross(-wall.normals[:, :, 0].T[:, None], corner_h)
    np.testing.assert_allclose(corners, expected, atol=1e-12)
    np.testing.assert_allclose(current, expected.mean(axis=1), atol=1e-12)

    # All edge and face degrees of freedom must have a continuous tangential trace.
    random = np.random.default_rng(42).normal(size=basis.N)
    left = InteriorFacetBasis(mesh, basis.elem, side=0)
    right = InteriorFacetBasis(mesh, basis.elem, side=1)
    jump = left.interpolate(random) - right.interpolate(random)
    np.testing.assert_allclose(np.cross(left.normals, jump, axis=0), 0., atol=2e-12)
