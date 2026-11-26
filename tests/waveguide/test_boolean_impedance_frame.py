"""An opaque Boolean frame must retain and classify its inner aperture wall."""
import numpy as np

from fem import ModeSolver2D, materials, shapes
import pytest


@pytest.mark.parametrize("material", [materials.PEC, materials.copper])
def test_boolean_ring_assigns_both_circular_walls_to_conductor(material):
    solver = ModeSolver2D(frequency=299792458., x_range=1., y_range=1.)
    center = (.5, .5)
    ring = shapes.Difference(shape=shapes.Circle(center=center, radius=.25),
                             tool=shapes.Circle(center=center, radius=.14))
    solver.add_geometry(shape=ring, material=material, name="ring")
    mesh = solver.mesh(max_element_size=.1)
    endpoints = mesh.nodes[mesh.mesh.facets[:, mesh.boundary_facets["ring"]]]
    radii = np.linalg.norm(endpoints - center, axis=-1)
    inner = np.isclose(radii, .14).all(axis=0)
    outer = np.isclose(radii, .25).all(axis=0)
    assert inner.any() and outer.any() and (inner | outer).all()


def test_copper_difference_frame_keeps_only_aperture_and_impedance_walls():
    width, height, wall = .02286, .01016, .001
    bounds = ((-wall, width+wall), (-wall, height+wall))
    solver = ModeSolver2D(frequency=10e9, x_range=bounds[0], y_range=bounds[1])
    solver.add_geometry(name='copper_frame', material=materials.copper,
        shape=shapes.Difference(shape=shapes.Rectangle(bounds=bounds),
            tool=shapes.Rectangle(bounds=((0., width), (0., height)))))
    mesh = solver.mesh(resolution=(12, 6), element_order=2, material_aware=False,
                       interface_refinement=None, boundary_refinement=None)
    centers = mesh.nodes[mesh.elements].mean(axis=1)
    assert np.all((centers[:, 0] > 0) & (centers[:, 0] < width))
    assert np.all((centers[:, 1] > 0) & (centers[:, 1] < height))
    np.testing.assert_array_equal(np.sort(mesh.boundary_facets['copper_frame']),
                                  np.sort(mesh.mesh.boundary_facets()))
    assert mesh.boundary_facets['outer_pec'].size == 0
