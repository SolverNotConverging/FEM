"""Adaptive API contracts and coarse-start integration for every FEM backend."""
from fem.common import SurfaceImpedance, materials, shapes

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.constants import c

from fem.adaptivity import run_adaptive, validate_controls
from fem.waveguide import ModeSolver1D, ModeSolver2D
from fem.periodic import PeriodicModeSolver2D, PeriodicModeSolver3D
from fem.common.errors import ConfigurationError as WaveConfigurationError


def modal_case(kind):
    if kind == "scalar":
        return ModeSolver1D(frequency=c, x_range=1.0), {"num_modes": 3, "neff_guess": .85}
    if kind == "vector":
        return ModeSolver2D(frequency=c, x_range=1.0, y_range=0.5), {"num_modes": 1, "neff_guess": .85}
    if kind == "periodic":
        return PeriodicModeSolver2D(frequency=c, x_range=1.0, z_range=0.25, polarization='TE'), {"num_modes": 1, "neff_guess": .85}


def diagnostics(result):
    return result.metadata if hasattr(result, "metadata") else result.solve_info


def test_periodic_tm_measures_electric_flux_jump():
    solver = PeriodicModeSolver2D(frequency=c, x_range=1., z_range=.25, polarization="TM",
                                  boundary=materials.PMC)
    solver.mesh(max_element_size=.2)
    result = solver.solve(backend="cpu", num_modes=1, neff_guess=.85,
                          max_refinements=1, adaptive_tolerance=1e-8)
    first, last = result.metadata["adaptive_history"]
    assert first["residual"] > 1e-4
    assert last["elements"] > first["elements"]
    assert last["residual"] < first["residual"]


@pytest.mark.gmsh
@pytest.mark.parametrize("kind", ["scalar", "vector", "periodic", "periodic3d"])
def test_modal_sweep_reuses_high_frequency_mesh(kind, capsys):
    if kind == "periodic3d":
        solver = PeriodicModeSolver3D(frequency=c, x_range=1., y_range=.5, z_range=.25)
        options = {"num_modes": 1, "neff_guess": .85}
    else:
        solver, options = modal_case(kind)
    low, high = solver.sweep(frequencies=[.95 * c, 1.05 * c], max_refinements=1,
                             adaptive_tolerance=1e-9, backend="cpu", **options)
    assert solver.frequency == c
    assert high.solve_info["max_refinements"] == 1
    assert low.solve_info["max_refinements"] == 0
    assert low.solve_info["mesh_adaptation_frequency"] == high.frequency
    np.testing.assert_array_equal(low.mesh_data.elements, high.mesh_data.elements)
    np.testing.assert_array_equal(low.mesh_data.coordinates, high.mesh_data.coordinates)
    assert low[0].beta != high[0].beta
    output = capsys.readouterr().err
    assert "Frequency sweep" in output and "2/2" in output
    assert "elements=" in output and "residual=" in output and "tolerance=" in output


@pytest.mark.gmsh
@pytest.mark.parametrize("kind", ["scalar", "vector", "periodic"])
def test_default_budget_refines_coarse_modes_and_reduces_residual(kind):
    solver, options = modal_case(kind)
    coarse = solver.solve(max_refinements=0, adaptive_tolerance=1e-5, **options)
    before = diagnostics(coarse)
    assert len(before["adaptive_history"]) == 1
    assert not before["adaptive_converged"]
    result = solver.solve(adaptive_tolerance=1e-5, **options)
    info = diagnostics(result)
    history = info["adaptive_history"]
    assert info["max_refinements"] == 2
    assert len(history) == 3
    assert history[-1]["status"] == "refinement_limit"
    assert history[-1]["elements"] > history[0]["elements"]
    assert info["adaptive_residual"] < before["adaptive_residual"]
    assert abs(result[0].neff - np.sqrt(.75)) < 3e-3
    if kind == "periodic":
        mesh = solver.mesh_data
        np.testing.assert_allclose(mesh.nodes[mesh.slave_nodes, 0], mesh.nodes[mesh.master_nodes, 0])
    stopped = solver.solve(adaptive_tolerance=100., **options)
    assert len(diagnostics(stopped)["adaptive_history"]) == 1
    assert diagnostics(stopped)["adaptive_converged"]


@pytest.mark.gmsh
def test_3d_periodic_coarse_start_regenerates_valid_constraints():
    solver = PeriodicModeSolver3D(frequency=10000000000.0, x_range=0.02, y_range=0.01, z_range=0.005, background_material=materials.Material(epsilon=2.25, mu=1.0))
    result = solver.solve(adaptive_tolerance=1e-05, num_modes=1, neff_guess=1.3)
    history = result.metadata["adaptive_history"]
    assert len(history) == 3
    assert history[-1]["status"] == "refinement_limit"
    assert result[0].residual < 1e-8
    expected = np.sqrt(2.25 - (np.pi / (.02 * solver.k0))**2)
    assert abs(result[0].neff - expected) < .01
    pairs = solver.mesh_data.periodic_node_pairs
    nodes = solver.mesh_data.nodes
    np.testing.assert_allclose(nodes[pairs[:, 0], :2], nodes[pairs[:, 1], :2], atol=1e-12)


@pytest.mark.parametrize("value", [-1, 1.5, np.nan, np.inf, True, np.bool_(True), "2", [], 1j])
def test_invalid_refinement_budget(value):
    with pytest.raises(ValueError, match="max_refinements"):
        validate_controls(value, .05)


@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf, True, np.bool_(True), "0.05", [], 1j, np.complex128(.05)])
def test_invalid_residual_threshold(value):
    with pytest.raises(ValueError, match="adaptive_tolerance"):
        validate_controls(2, value)


def test_manual_refinement_after_adaptation_retains_fine_intervals():
    solver, options = modal_case("scalar")
    solver.solve(max_refinements=3, adaptive_tolerance=1e-5, **options)
    previous_minimum = np.diff(solver.mesh_data.nodes).min()
    refined = solver.refine()
    assert np.diff(refined.nodes).max() <= previous_minimum / 2 * (1 + 1e-12)


@pytest.mark.parametrize("kind", ["scalar", "vector", "periodic"])
def test_public_controls_are_validated_before_solving(kind):
    solver, options = modal_case(kind)
    with pytest.raises((ValueError, WaveConfigurationError), match="max_refinements"):
        solver.solve(max_refinements=-1, **options)


def test_failed_solve_consumes_budget_and_remesh_failure_restores_state():
    @dataclass
    class Result:
        metadata: dict

    owner = SimpleNamespace(mesh="coarse", solution=None)
    def solve():
        if owner.mesh == "coarse":
            raise RuntimeError("coarse eigensolve failed")
        owner.solution = Result({})
        return owner.solution

    def refine(indicators):
        owner.mesh = "fine"

    result = run_adaptive(owner, solve, lambda result: (np.ones(2), .1), refine,
                          max_refinements=1, adaptive_tolerance=.1, retry_errors=(RuntimeError,))
    assert result is owner.solution
    assert [step["status"] for step in result.metadata["adaptive_history"]] == ["solve_failed", "tolerance"]
    def fail(indicators):
        owner.mesh = "broken"
        raise RuntimeError("remeshing failed")

    with pytest.raises(RuntimeError, match="remeshing failed"):
        run_adaptive(owner, lambda: result, lambda result: (np.ones(2), 1.), fail,
                     max_refinements=1, adaptive_tolerance=.1)
    assert owner.mesh == "fine"


@pytest.mark.parametrize("order", [1, 2])
def test_dielectric_edge_between_fixed_exterior_vertices_stays_active(order):
    from fem import ElectrostaticSolver, Material
    solver = ElectrostaticSolver(dim=2, x_range=1., y_range=1.)
    solver.add_geometry(shape=shapes.Rectangle(bounds=((.5, 1.), (0., 1.))), material=Material(epsilon=4))
    solver.set_potential(geometry="right", potential=1.)
    mesh = solver.mesh(max_element_size=2., element_order=order, material_aware=False,
                       interface_refinement=None, boundary_refinement=None)
    result = solver.solve(max_refinements=0, backend="cpu")
    interface = np.flatnonzero(np.all(np.isclose(mesh.nodes[mesh.mesh.facets, 0], .5), axis=0))
    assert len(interface) == 1 and solver.fixed_mask[mesh.mesh.facets[:, interface]].all()
    assert not solver._prescribed_facets[interface].any()
    if order == 2:
        _, _, basis = solver._assemble()
        fixed, *_ = solver._dirichlet_data(basis)
        assert not fixed[basis.facet_dofs[:, interface]].any()
    assert result.adaptive_residual > .05
    assert len(solver.solve(max_refinements=1, backend="cpu").elements) > len(mesh.elements)


def test_thin_periodic_pml_is_sampled_and_changes_the_pencil():
    from fem.periodic.assembly_2d import assemble_periodic_system_2d
    solver = PeriodicModeSolver2D(frequency=c, x_range=1., z_range=1., polarization="TE")
    solver.add_pml(thickness=.001, direction="x")
    mesh = solver.mesh()
    x = solver.system.basis.global_coordinates()[0] / solver.k0
    for a, b in zip(np.linspace(0, .001, 4)[:-1], np.linspace(0, .001, 4)[1:]):
        assert ((a < x) & (x < b)).all(axis=1).any()
        assert ((1-b < x) & (x < 1-a)).all(axis=1).any()
    vacuum = assemble_periodic_system_2d(mesh, polarization="TE", frequency=c, k0=solver.k0,
                                         material_at=solver.geometry.material_at)
    assert np.linalg.norm((solver.system.A0 - vacuum.A0).data) > 0.


def test_curved_interface_trace_uses_the_adjacent_cell_tag():
    from skfem import InteriorFacetBasis
    from fem import Material
    solver = ModeSolver2D(frequency=c, x_range=1., y_range=1.)
    solver.add_geometry(shape=shapes.Circle(center=(.5, .5), radius=.2), material=Material(epsilon=4))
    mesh = solver.mesh(max_element_size=.1)
    first, second = mesh.mesh.f2t
    ids = np.flatnonzero((second >= 0) & (mesh.element_tags[first] != mesh.element_tags[second]))
    assert len(ids)
    for side in (0, 1):
        space = InteriorFacetBasis(solver.system.basis.mesh, solver.system.basis.elem, facets=ids, side=side)
        eps, _ = solver.system.material_at.on_cells(*(space.global_coordinates()/solver.k0), cells=space.tind)
        expected = np.where(mesh.element_tags[space.tind] == 2, 4., 1.)
        assert np.all(eps == expected[None, :, None])


def test_refined_interval_tags_survive_archive(tmp_path):
    from fem import Material, load
    solver, options = modal_case("scalar")
    solver.add_geometry(shape=shapes.Interval(bounds=(.55, .8)), material=Material(epsilon=4))
    solver.mesh(max_element_size=.1)
    result = solver.solve(**options, backend="cpu", max_refinements=1, adaptive_tolerance=1e-8)
    path = tmp_path / "guide.h5"
    result.save(path)
    for mesh in (result.mesh_data, load(path).mesh_data):
        center = mesh.coordinates[mesh.elements].mean(axis=1)[:, 0]
        np.testing.assert_array_equal(mesh.metadata["element_tags"], 1 + ((center > .55) & (center < .8)))


@pytest.mark.parametrize("limit,status", [("max_elements", "element_limit"), ("max_dofs", "dof_limit")])
def test_modal_budget_retains_last_mesh_and_reports_quality(limit, status):
    solver, options = modal_case("vector")
    mesh = solver.mesh(max_element_size=.25)
    result = solver.solve(**options, backend="cpu", adaptive_tolerance=1e-7,
                          **{limit: len(mesh.elements) if limit == "max_elements" else solver.system.basis.N})
    final, = result.solve_info["adaptive_history"]
    assert final["status"] == status and solver.mesh_data is mesh
    assert 0 < final["minimum_quality"] <= final["quality_p05"] <= 1
    assert np.linalg.norm(list(final["residual_components"].values())) == pytest.approx(final["residual"])


def test_anisotropic_background_uses_cross_component_wavelength():
    from fem import Material
    solver = PeriodicModeSolver2D(frequency=c, x_range=.012, z_range=.012, polarization="TE",
                                  background_material=Material(epsilon=(1, 100, 1), mu=(1, 1, 100)))
    assert solver.mesh().info.requested_maximum_edge == pytest.approx(.0025)


def test_curvature_tolerance_reduces_area_and_capacitance_error():
    from fem import ElectrostaticSolver, EPSILON_0
    exact = 2 * np.pi * EPSILON_0 / np.log(.4 / .16)
    errors = []
    for tolerance in (.01, .0002):
        solver = ElectrostaticSolver(dim=2, x_range=1., y_range=1.)
        disk = shapes.Circle(center=(.5, .5), radius=.16)
        ground = shapes.Difference(shape=shapes.Rectangle(bounds=((0., 1.), (0., 1.))),
                                   tool=shapes.Circle(center=(.5, .5), radius=.4))
        solver.set_potential(geometry=disk, potential=1., name="core")
        solver.set_potential(geometry=ground, potential=0.)
        mesh = solver.mesh(max_element_size=.15, interface_refinement=None, boundary_refinement=None,
                           geometry_tolerance=tolerance)
        result = solver.solve(max_refinements=0, backend="cpu")
        owned = mesh.potential_cells[solver.geometry.potentials[0].id]
        vertices = mesh.nodes[mesh.elements[owned]]
        area = .5 * abs(np.linalg.det(vertices[:, 1:] - vertices[:, :1])).sum()
        errors.append([abs(area / (np.pi * .16**2) - 1), abs(2 * result.energy / exact - 1)])
    assert np.all(np.array(errors[1]) < errors[0])


def test_pml_layer_convergence_against_stretched_guide():
    solver = PeriodicModeSolver2D(frequency=c, x_range=1., z_range=.15, polarization="TE")
    thickness, sigma, order = .1, 3., 3
    solver.add_pml(thickness=thickness, direction="x", sigma_max=sigma, order=order)
    # The two polynomial stretches give this exact complex plate separation.
    separation = 1 - 2j * thickness * sigma / (order + 1)
    exact = np.sqrt(1 - (.5 / separation)**2)
    errors = []
    for n in (1, 4):
        solver.mesh(max_element_size=.14, pml_layers=n)
        result = solver.solve(num_modes=1, neff_guess=exact, backend="cpu", max_refinements=0)
        errors.append(abs(result[0].neff - exact))
        step = result.solve_info["adaptive_history"][0]
        assert 0 < step["max_mode_pml_fraction"] < 1
        assert np.hypot(step["physical_residual"], step["pml_residual"]) == pytest.approx(step["residual"])
    assert errors[1] < .1 * errors[0]


def test_smooth_field_with_no_jumps_has_nonzero_volume_residuals():
    from skfem import MeshTri, Basis, ElementTriP2
    from fem.adaptivity.policy import maxwell_residual
    basis = Basis(MeshTri().refined(), ElementTriP2(), intorder=5)
    def fields(space, coefficient):
        x, y = space.global_coordinates()
        electric = np.stack((x, y, 0*x)).astype(complex)
        return electric, np.zeros_like(electric)
    def material(cells, x, y):
        unit = np.ones((3, *x.shape))
        return unit, unit
    _, residual, parts = maxwell_residual(basis, fields, material, axes=(0, 1), neff=0)
    assert parts["interface"] < 1e-12
    assert parts["element_pde"] > 0 and parts["divergence"] > 0
    assert residual == pytest.approx(np.linalg.norm(list(parts.values())))
