"""Regression tests for selection and application of the calibrated driven defaults."""

import numpy as np
import pytest

import fem
from benchmarks.analytical.amr_observables import two_port_guide
from benchmarks.analytical.patch_marking_weights import inset_patch
from fem.scattering.mesh_policy import driven_mesh_policy


@pytest.mark.parametrize(
    "length,eps,name,density,fraction,factor,backend",
    [
        (0.55, 1.0, "air guide", 2.5, 0.5, 2.25, "gmsh"),
        (0.275, 1.0, "short air guide", 2.75, 0.5, 1.7, "gmsh"),
        (0.55, 1.44, "dielectric guide", 3.25, 0.12, 1.1, "auto"),
        (0.275, 1.44, "short dielectric guide", 3.0, 0.5, 2.1, "gmsh"),
    ],
)
def test_uniform_guide_defaults(length, eps, name, density, fraction, factor, backend):
    solver = two_port_guide(length, eps)
    data = solver.mesh()
    assert data.info["mesh_profile"] == name
    assert data.info["cad_point_size"] == data.info["size"] == pytest.approx(1 / density)
    assert (solver._marking_fraction, solver._refinement_factor, solver._refinement_backend) == (
        fraction,
        factor,
        backend,
    )


@pytest.mark.parametrize("eps,wavelength", [(2.25, 1.0), (1.44 - 0.1j, 1.0), (1.0, 0.6)])
def test_uniform_guide_uses_general_settings_outside_calibration(eps, wavelength):
    solver = two_port_guide(epsilon=eps)
    profile = driven_mesh_policy(solver._regions, tuple(solver._ports.values()), 3, wavelength, 2)
    assert profile.name == "uniform guide"
    assert profile.density is None and not profile.fitted_weights


def test_default_port_convergence_and_manual_overrides():
    solver = two_port_guide()
    solver.configure_adaptivity(minimum_quality=0.0)
    result = solver.solve(backend="cpu")
    history = result.solve_info["adaptive_history"]
    assert result.solve_info["max_refinements"] == 8
    assert len(history) == 3 and history[-1]["status"] == "tolerance"
    assert history[-1]["s_parameter_delta"] < 0.02
    assert history[0]["refinement_backend"] == "gmsh"

    solver.mesh(wavelength_elements=5)
    solver.configure_adaptivity(
        refinement_backend="auto", marking_fraction=0.3, refinement_factor=1.4, marking_weights={}
    )
    result = solver.solve(backend="cpu", max_refinements=0, adaptive_tolerance=10.0)
    (record,) = result.solve_info["adaptive_history"]
    assert record["marking_weights"] == {} and "s_parameter_delta" not in record
    assert (record["marking_fraction"], record["refinement_factor"]) == (0.3, 1.4)
    assert result.mesh_data.info["size"] == 0.2
    assert result.mesh_data.info["cad_point_size"] is None


def test_printed_antenna_and_dielectric_scattering():
    solver = inset_patch()
    initial = solver.mesh()
    assert initial.info["mesh_profile"] == "printed conductor"
    assert initial.info["automatic_local_sizes"][0] < initial.info["size"]
    assert solver._marking_weights == {"interface": 4.0, "element_pde": 4.0}
    assert solver.mesh(max_element_size=0.015).info["automatic_local_sizes"] == []
    solver.configure_adaptivity(marking_weights={})
    assert solver._marking_weights == {}

    scatter = fem.ScatteringSolver3D(
        frequency=fem.C_0, x_range=(-1, 1), y_range=(-1, 1), z_range=(-1, 1)
    )
    scatter.add_geometry(
        shape=fem.shapes.Box(bounds=((-0.3, 0.3), (-0.25, 0.25), (-0.02, 0.02))),
        material=fem.Material(epsilon=4),
    )
    data = scatter.mesh()
    assert data.info["mesh_profile"] == "dielectric slab"
    assert data.info["automatic_local_sizes"] == pytest.approx([1 / 12])
    assert scatter._marking_weights["interface"] == 4


@pytest.mark.parametrize(
    "ground,name,weighted", [(True, "monopole", True), (False, "dipole", False)]
)
def test_gap_electrodes_select_dipole_or_monopole(ground, name, weighted):
    antenna = fem.AntennaSolver3D(
        frequency=fem.C_0, x_range=(-1, 1), y_range=(-1, 1), z_range=(-1, 1)
    )

    def box(x, y, z):
        return fem.shapes.Box(bounds=(x, y, z))

    bounds = (-0.01, 0.01)
    bottom = (
        box((-0.2, 0.2), (-0.2, 0.2), (-0.02, -0.01))
        if ground
        else box(bounds, bounds, (-0.25, -0.01))
    )
    for electrode in (bottom, box(bounds, bounds, (0.01, 0.25))):
        antenna.add_geometry(shape=electrode, material=fem.materials.PEC)
    antenna.add_discrete_port(name="drive", axis="z", gap=box(bounds, bounds, (-0.01, 0.01)))
    policy = driven_mesh_policy(antenna._regions, tuple(antenna._ports.values()), 3, 1.0, 2)
    assert policy.name == name and policy.fitted_weights == weighted
    assert all(np.isclose(h, 1 / 12) for _, h in policy.local_sizes)


def test_guide_with_internal_obstacle_uses_discontinuity_policy():
    solver = two_port_guide()
    solver.add_geometry(
        shape=fem.shapes.Box(bounds=((-0.36, -0.25), (-0.15, 0.15), (-0.02, 0.02))),
        material=fem.materials.PEC,
    )
    mesh = solver.mesh()
    assert mesh.info["mesh_profile"] == "guide discontinuity" and mesh.info[
        "size"
    ] == pytest.approx(1 / 6)
    assert solver._marking_weights["interface"] == 4
