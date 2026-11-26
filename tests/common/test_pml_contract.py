import numpy as np
import pytest

import fem


@pytest.mark.parametrize(("solver_type", "axes", "waveguide"), [
    (fem.ModeSolver1D, "x", True),
    (fem.ModeSolver2D, "xy", True),
    (fem.PeriodicModeSolver2D, "xz", False),
    (fem.PeriodicModeSolver3D, "xyz", False),
    (fem.ScatteringSolver2D, "xy", False),
    (fem.ScatteringSolver3D, "xyz", False),
    (fem.AntennaSolver3D, "xyz", False),
])
@pytest.mark.parametrize("layer_first", [True, False])
def test_material_entering_one_sided_pml(solver_type, axes, waveguide, layer_first):
    solver = solver_type(frequency=fem.C_0, **{axis + "_range": 2. for axis in axes})
    material = fem.Material(epsilon=2.25)
    limits = [(.1, .6)] + [(.5, 1.5)] * (len(axes) - 1)
    shape = (fem.shapes.Interval, fem.shapes.Rectangle, fem.shapes.Box)[len(axes) - 1]

    def geometry():
        if axes == "x":
            return solver.add_layer(x_range=limits[0], material=material)
        return solver.add_geometry(shape=shape(bounds=limits), material=material)

    def absorber():
        pml = solver.add_pml(thickness="auto", direction="x-", order=2, sigma_max=3.)
        assert pml.thickness == pytest.approx(1/3)
        assert pml.stretch(pml.thickness / 2) == pytest.approx(1 - .75j)

    first, second = (absorber, geometry) if layer_first else (geometry, absorber)
    first()
    if waveguide:
        second()
        sample = (1/6,) if axes == "x" else (1/6, 1.)
        eps, mu = solver.geometry.transformed_material_at(*sample)
        stretch = 1 - .75j
        np.testing.assert_allclose(mu, [1/stretch, stretch, stretch])
        np.testing.assert_allclose(eps, 2.25 * mu)
    else:
        with pytest.raises(fem.GeometryError, match="vacuum"):
            second()


def test_independent_directional_profiles_reach_the_specified_outer_faces():
    solver = fem.ScatteringSolver3D(frequency=fem.C_0, x_range=2., y_range=2., z_range=2.)
    solver.add_pml(thickness=.2, direction="x-")
    solver.add_pml(thickness=.3, order=2, sigma_max=2., direction="z+")
    with pytest.raises(fem.GeometryError, match="already covers"):
        solver.add_pml(direction="x")
    snapshot = solver.mesh(max_element_size=.5)
    np.testing.assert_allclose(snapshot.coordinates.min(axis=0), [.2, 0., 0.], atol=1e-12)
    np.testing.assert_allclose(snapshot.coordinates.max(axis=0), [2., 2., 1.7])
    np.testing.assert_allclose(solver._mesh.p.min(axis=1), 0., atol=1e-12)
    np.testing.assert_allclose(solver._mesh.p.max(axis=1), 2.)
    assert set(snapshot.info["pml_profiles"]) == {"x-", "z+"}
