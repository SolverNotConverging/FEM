import os

import numpy as np
import pytest

import fem
from fem.electrostatics.solver import EPSILON_0


def test_cross_component_permittivity_in_tetrahedral_stiffness():
    eps = np.array([[3., .3, -.25], [.3, 2., .4], [-.25, .4, 5.]])
    solver = fem.ElectrostaticSolver(dim=3, background_material=fem.Material(epsilon=eps))
    solver.mesh(max_element_size=.4, boundary_refinement=None)
    matrix, _, basis = solver._assemble()
    gradient = np.array([1., -2., 3.])
    potential = gradient @ basis.doflocs
    expected = EPSILON_0 * (gradient @ eps @ gradient)
    assert potential @ (matrix @ potential) == pytest.approx(expected, rel=1e-12, abs=0.)


@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("backend", [
    "cpu", pytest.param("gpu", marks=pytest.mark.skipif(
        os.environ.get("FEM_TEST_CUDA") != "1", reason="set FEM_TEST_CUDA=1 to run cuDSS checks")),
])
def test_layered_capacitor_charge_energy_and_fields(order, backend):
    gap, area, volts = .004, .006 * .008, 3.
    solver = fem.ElectrostaticSolver(dim=3, x_range=.006, y_range=.008,
                                     z_range=gap, outer_potential=None)
    solver.add_box(x_range=solver.x_range, y_range=solver.y_range, z_range=(gap/2, gap),
                   material=fem.Material(epsilon=(2., 3., 4.)), name="dielectric")
    solver.set_potential(geometry="z-", potential=0., name="ground")
    solver.set_potential(geometry="z+", potential=volts, name="drive")
    mesh = solver.mesh(max_element_size=.002, element_order=order, interface_refinement=None,
                       boundary_refinement=None, material_aware=False)
    result = solver.solve(backend=backend)
    capacitance = EPSILON_0 * area / (gap/2 + gap/8)
    assert result.energy == pytest.approx(.5 * capacitance * volts**2, rel=1e-10, abs=0.)
    assert result.conductor_charge("drive") == pytest.approx(capacitance * volts, rel=1e-10, abs=0.)
    assert result.conductor_charge("ground") == pytest.approx(-capacitance * volts, rel=1e-10, abs=0.)
    field = result.element_electric_field
    np.testing.assert_allclose(field[:, :2], 0., atol=1e-9)
    np.testing.assert_allclose(field[:, 2], -volts / (gap/2 + gap/8)
                               / np.where(mesh.element_tags == 1, 1., 4.), rtol=1e-10)
    assert result.adaptive_residual < 1e-10
    np.testing.assert_array_equal(solver.compute_electric_field(), result.electric_field.T)


@pytest.mark.parametrize("order", [1, 2])
def test_volume_charge_and_tetrahedral_adaptation(order):
    solver = fem.ElectrostaticSolver(dim=3, outer_potential=None)
    solver.set_potential(geometry="front", potential=0., name="a")
    solver.set_potential(geometry="back", potential=0., name="b")
    solver.add_charge_density(geometry=fem.shapes.Box(bounds=((0., 1.),)*3), density=EPSILON_0)
    solver.mesh(element_order=order, max_element_size=.35, boundary_refinement=None)
    result = solver.solve(backend="cpu", max_refinements=1)
    total = sum(result.conductor_charges.values())
    assert total == pytest.approx(-EPSILON_0, rel=1e-10, abs=0.)
    if order == 1:
        first, second = result.adaptive_history
        assert second["elements"] > first["elements"]
        assert second["residual"] < first["residual"]
    else:
        z = result.coordinates[:, 2]
        np.testing.assert_allclose(result.potential, z * (1. - z) / 2, atol=1e-12)
        assert result.adaptive_residual < 1e-10


def test_curved_conductor_archive_and_viewer(tmp_path):
    solver = fem.ElectrostaticSolver(dim=3)
    electrode = solver.add_sphere(center=(.5, .5, .5), radius=.22, material=fem.materials.PEC)
    solver.set_potential(geometry=electrode, potential=1.)
    solver.mesh(max_element_size=.28, geometry_tolerance=.015, boundary_refinement=None)
    result = solver.solve(max_refinements=1, backend="cpu")
    assert len(result.adaptive_history) == 2
    owned = next(iter(result.mesh.potential_cells.values()))
    indices = np.unique(result.elements[owned])
    np.testing.assert_allclose(result.potential[indices], 1.)
    assert result.conductor_charge(electrode.name) > 0
    loaded = fem.load(result.save(tmp_path / "sphere.h5"))
    np.testing.assert_array_equal(loaded.element_electric_field, result.element_electric_field)
    viewer = loaded.show(block=False)
    try:
        for name in ("potential", "E", "Ez", "|E|", "D", "Dz"):
            viewer._component.setCurrentText(name)
            assert not viewer._status.currentMessage()
            assert viewer._plotter is not None
        viewer._mesh.setChecked(True)
        viewer._material.setChecked(True)
        assert not viewer._status.currentMessage()
        assert "boundary_pec" in viewer._plotter.actors
        assert viewer._component.findText("H") == -1
    finally:
        viewer._window.close()
    plot = loaded.plot(component="E")
    try:
        assert "|E| (V/m)" in plot.scalar_bars
    finally:
        plot.close()
