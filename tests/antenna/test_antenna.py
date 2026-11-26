import os

import numpy as np
import pytest

import fem


def guide(*, open_end=False, higher_mode=False, wall=fem.materials.PEC, pml_layers=5, **mesh_options):
    u, v = (-0.36, 0.36), (-0.15, 0.15)
    s = fem.AntennaSolver3D(
        frequency=(1.6 if higher_mode else 1) * fem.C_0,
        x_range=(-1.0, 1.0),
        y_range=(-0.85, 0.85),
        z_range=(-0.95, 0.95),
        pml_layers=pml_layers,
    )
    outer = fem.shapes.Box(bounds=((-0.45, 0.45), (-0.27, 0.27), (-0.37, 0.37)))
    bore = fem.shapes.Box(bounds=(u, v, (-0.275, 0.42 if open_end else 0.275)))
    s.add_geometry(shape=fem.shapes.Difference(shape=outer, tool=bore), material=wall)
    for name, z in [("left", -0.275)] if open_end else [("left", -0.275), ("right", 0.275)]:
        s.add_waveguide_port(
            name=name,
            axis="z",
            position=z,
            cross_section=fem.shapes.Rectangle(bounds=(u, v)),
            num_modes=None if higher_mode and z < 0 else 3,
            mode=2 if higher_mode and z < 0 else 1,
            amplitude=1 if z < 0 else 0.5,
            phase_deg=0 if z < 0 else 60,
        )
    s.mesh(**{"max_element_size": 0.13, **mesh_options})
    return s


def test_two_port_guide_s_and_coherent_excitation(tmp_path, monkeypatch):
    s = guide()
    solves = []
    solve = s._solve

    def record(frequency, backend, basis, *args, **kwargs):
        solves.append((frequency, basis.mesh, kwargs.get("estimate", False)))
        return solve(frequency, backend, basis, *args, **kwargs)

    monkeypatch.setattr(s, "_solve", record)
    sweep = s.sweep(frequencies=np.array([1.0, 1.05]) * fem.C_0, backend="cpu",
                    max_refinements=1, adaptive_tolerance=1e-9)
    low, high = sweep.frequencies
    assert [(f, estimate) for f, _, estimate in solves] == [(high, True), (high, True), (low, False)]
    assert solves[1][1] is solves[2][1]
    assert sweep[0].mesh_data is sweep[1].mesh_data
    assert sweep[1].solve_info["adaptive_history"][0]["refinement_backend"] == "mmg3d"
    assert sweep[1].solve_info["adaptive_history"][-1]["status"] == "refinement_limit"
    assert "adaptive_history" not in sweep[0].solve_info
    assert sweep.s_parameters.shape == (2, 2, 2)
    for case in sweep:
        beta = np.sqrt((2 * np.pi * case.frequency / fem.C_0) ** 2 - (np.pi / 0.72) ** 2)
        matrix = case.s_parameters
        assert abs(matrix[0, 0]) < 0.08
        assert abs(matrix[1, 0] - np.exp(-1j * beta * 0.55)) < 0.008
        assert case.port_modes["left"].mesh_data.info["element_order"] == 2
        np.testing.assert_allclose(matrix, matrix.T, atol=1e-10)
        assert case.incident_power == pytest.approx(1.25)
        assert abs(case.accepted_power) < 1e-8
        assert case.radiated_power < 1e-15
        np.testing.assert_allclose(case.active_s_parameters, matrix @ case.incident / case.incident)
    loaded = fem.load(sweep.save(tmp_path / "two_ports.h5"))
    assert all(case.solve_info["mesh_adaptation_frequency"] == high for case in loaded)
    assert loaded[0].mesh_data.info["pml_layers"] == 5
    assert loaded[0].mesh_data.info["pml_thickness"] == pytest.approx(1 / 3)
    assert len(loaded[0].electric_field) == len(loaded[0].mesh_data.elements)
    np.testing.assert_array_equal(loaded[0].magnetic_field, sweep[0].magnetic_field)
    case = loaded[0]
    np.testing.assert_array_equal(case.surface_current_vertices, sweep[0].surface_current_vertices)
    np.testing.assert_allclose(case.surface_current_vertices.mean(axis=1), case.surface_current)
    np.testing.assert_array_equal(loaded.s_parameters, sweep.s_parameters)
    assert len(loaded.plot().axes) == 2
    assert loaded.plot(quantity="phase").axes[1].get_ylabel() == "Phase (degrees)"
    s.set_excitation(port="right", amplitude=0)
    result = s.solve(backend="cpu", max_refinements=0)
    assert np.isnan(result.active_s_parameters[1])
    assert result.active_s_parameters[0] == pytest.approx(result.s_parameters[0, 0])


@pytest.mark.parametrize("order, tolerance", [(1, .055), (2, .003)])
def test_through_guide_against_analytical_transmission(order, tolerance):
    solver = guide(max_element_size=.16, element_order=order)
    case = solver.solve(backend="cpu", max_refinements=0)
    beta = np.sqrt((2 * np.pi)**2 - (np.pi / .72)**2)
    assert abs(case.s_parameters[1, 0] - np.exp(-1j * beta * .55)) < tolerance
    assert case.port_modes["left"].mesh_data.info["element_order"] == order


def test_radiating_guide_nf2ff():
    s = guide(open_end=True, element_order=1)
    r = s.solve(backend="cpu")
    pattern = r.far_field([[0, 0, 1], [0, 0, -1]])
    assert 0 < r.accepted_power <= r.incident_power
    assert pattern.intensity[0] > pattern.intensity[1]
    assert r.radiated_power / r.accepted_power == pytest.approx(1, abs=0.25)
    np.testing.assert_allclose(
        pattern.realized_gain, pattern.gain * r.accepted_power / r.incident_power
    )
    np.testing.assert_allclose(
        pattern.directivity, pattern.gain * r.accepted_power / r.radiated_power
    )
    np.testing.assert_allclose(pattern.electric[:, 2], 0, atol=1e-12)
    s.set_excitation(port="left", amplitude=0.5, phase_deg=90)
    shifted = s.solve(backend="cpu", max_refinements=0)
    scaled = shifted.far_field([[0, 0, 1], [0, 0, -1]])
    np.testing.assert_allclose(scaled.electric, 0.5j * pattern.electric, atol=1e-11)
    np.testing.assert_allclose(scaled.realized_gain, pattern.realized_gain)


def test_open_guide_checks_s_matrix_on_fixed_physical_mesh():
    solver = guide(open_end=True, pml_layers=1, element_order=1, max_element_size=.18)
    result = solver.solve(backend="cpu", adaptive_tolerance=10., pml_tolerance=.5, max_refinements=1)
    first, final = result.solve_info["adaptive_history"]
    assert first["status"] == "pml_refine"
    assert [first["pml_layers"], final["pml_layers"]] == [1, 2]
    assert first["elements"] == final["elements"] and final["mesh_elements"] > first["mesh_elements"]
    assert final["pml_observable"] == "S" and final["pml_interface_residual"] > 0
    np.testing.assert_array_equal(final["observable_values"], result.s_parameters)
    change = np.linalg.norm(final["observable_values"] - first["observable_values"])
    norm = max(np.linalg.norm(row["observable_values"]) for row in (first, final))
    assert final["pml_observable_change"] == pytest.approx(change / norm)


def test_selected_mode_must_propagate():
    with pytest.raises(fem.ConfigurationError, match="below cutoff"):
        guide().sweep(frequencies=[0.4 * fem.C_0], backend="cpu")


def test_viewer_inspects_port_modes_and_animates_vectors():
    import pyvista as pv
    from PySide6.QtTest import QSignalSpy

    view = guide().solve(backend="cpu", max_refinements=0).show(block=False)
    try:
        view._inspect.click()
        port = view._port_windows[0]
        assert port._mode.value() == 1
        port._mode.setValue(2)
        assert "neff" in port._info.text()
        port._window.close()
        assert view._component.findText("Total H") == view._component.findText("Total E") + 1
        pose = lens = None
        for component in ("Total E", "Total H", "Surface current"):
            view._component.setCurrentText(component)
            if pose is not None:
                view._mesh.click()
                view._material.click()
                view._app.processEvents()
                current = view._plotter.camera
                np.testing.assert_allclose((current.position, current.focal_point, current.up), pose)
                np.testing.assert_allclose((current.view_angle, current.parallel_scale, *current.clipping_range), lens)
            assert not view._status.currentMessage()
            view._phase.setValue(0)
            actor = view._plotter.actors["antenna_vectors"]
            first = pv.wrap(actor.mapper.GetInput()).points.copy()
            initial_scalars = pv.wrap(actor.mapper.GetInput()).active_scalars.copy()
            camera = view._plotter.camera
            camera.SetParallelProjection(component == "Surface current")
            camera.Azimuth(27)
            camera.Elevation(14)
            camera.Dolly(1.5)
            camera.Zoom(2)
            pose = (camera.position, camera.focal_point, camera.up)
            lens = (camera.view_angle, camera.parallel_scale, *camera.clipping_range)
            view._phase.setValue(180)
            second = pv.wrap(actor.mapper.GetInput()).points
            assert not np.allclose(first, second)
            view._phase.setValue(90)
            assert not np.allclose(initial_scalars, pv.wrap(actor.mapper.GetInput()).active_scalars)
            change = QSignalSpy(view._phase.valueChanged)
            view._play.click()
            assert change.wait(2000)
            view._play.click()
            assert view._phase.value() != 90
            np.testing.assert_allclose((camera.position, camera.focal_point, camera.up), pose)
            np.testing.assert_allclose(
                (camera.view_angle, camera.parallel_scale, *camera.clipping_range), lens
            )
            assert view._plotter.actors["antenna_vectors"] is actor
        view._component.setCurrentText("S-parameters")
        assert not view._timer.isActive()
    finally:
        view._window.close()


def test_pml_shell_matches_physical_mesh_and_outer_box():
    from skfem import MeshTet
    from fem.common.pml import PMLSpec
    from fem.common.radiation_pml import extend_box

    box = ((-1.0, 1.0), (-0.4, 0.6), (1.0, 2.2))
    physical = MeshTet.init_tensor(*(np.linspace(lo, hi, 3) for lo, hi in box))
    width = 0.3
    profile = PMLSpec(width)
    mesh = extend_box(physical, box, dict.fromkeys(profile.sides("xyz"), profile), layers=3)
    np.testing.assert_array_equal(mesh.t[:, : physical.nelements], physical.t)
    np.testing.assert_array_equal(mesh.p[:, : physical.nvertices], physical.p)
    limits = np.asarray(box) + [-width, width]
    faces = mesh.p[:, mesh.facets[:, mesh.boundary_facets()]]
    on_outer_plane = np.isclose(faces, limits[:, 0, None, None]).all(axis=1) | np.isclose(
        faces, limits[:, 1, None, None]
    ).all(axis=1)
    assert on_outer_plane.any(axis=0).all()
    vertices = mesh.p[:, mesh.t].transpose(2, 1, 0)
    volumes = abs(np.linalg.det(vertices[:, 1:] - vertices[:, :1])) / 6
    assert np.min(volumes) > 0
    assert volumes.sum() == pytest.approx(np.prod(np.diff(limits, axis=1)))


def test_pml_material_stretches_faces_edges_and_corners():
    from fem.common.radiation_pml import material_diagonal
    from fem.common.pml import PMLSpec

    points = np.array([[0.0, 1.5, -1.5, 1.5], [0.0, 0.0, 1.5, -1.5], [0.0, 0.0, 0.0, 1.5]])
    profile = PMLSpec(1, order=2, sigma_max=4)
    s = 1 - 1j
    expected = np.array([[1, 1, 1], [1 / s, s, s], [1, 1, s * s], [s, s, s]]).T
    profiles = dict.fromkeys(profile.sides("xyz"), profile)
    np.testing.assert_allclose(material_diagonal(points, [(-1, 1)] * 3, profiles), expected)


@pytest.mark.parametrize("angle", [0, 45, 70])
def test_pml_plane_wave_reflection_at_oblique_incidence(angle):
    """Planar TE reduction of Maxwell's equations, with a PEC-backed x layer."""
    from scipy.sparse.linalg import spsolve
    from skfem import Basis, BilinearForm, ElementLineP1, MeshLine, asm
    from fem.common.radiation_pml import material_diagonal
    from fem.common.pml import PMLSpec

    k, width = 2 * np.pi, 1 / 3
    tangential_k = k * np.sin(np.deg2rad(angle))
    normal_k = k * np.cos(np.deg2rad(angle))
    basis = Basis(MeshLine(np.linspace(0, width, 81)), ElementLineP1(), intorder=4)
    x = basis.global_coordinates()[0]
    coordinates = np.stack((x, np.zeros_like(x), np.zeros_like(x)))
    profile = PMLSpec(width, direction="x+")
    eps = material_diagonal(coordinates, [(-1, 0), (-1, 1), (-1, 1)], {"x+": profile})

    @BilinearForm(dtype=complex)
    def te(u, v, w):
        return u.grad[0] * v.grad[0] / eps[2] + (tangential_k**2 / eps[0] - k**2 * eps[1]) * u * v

    matrix = asm(te, basis)
    field = np.zeros(basis.N, dtype=complex)
    field[0] = 1  # prescribe total E at the interface, E=0 on the outer wall
    field[1:-1] = spsolve(matrix[1:-1, 1:-1], -matrix[1:-1, 0].toarray().ravel())
    flux = (matrix @ field)[0]
    reflected = (1j * normal_k - flux) / (1j * normal_k + flux)
    expected = -np.exp(-2j * normal_k * width - 2 * normal_k * width * profile.sigma_max / (profile.order + 1))
    assert abs(reflected - expected) < 4e-4


def test_wave_port_inherits_impedance_walls():
    solver = guide(wall=fem.SurfaceImpedance(impedance=2 + 2j))
    mode = solver.solve_port_modes(port="left")[0]
    assert mode.neff.real > 0 and mode.neff.imag < 0
    case = solver.solve(backend="cpu", max_refinements=0)
    assert case.accepted_power > 0
    assert abs(case.s_parameters[1, 0]) < 1


def test_dielectric_coax_port_on_annulus(tmp_path):
    solver = fem.AntennaSolver3D(
        frequency=fem.C_0, x_range=(-0.83, 0.83), y_range=(-0.83, 0.83), z_range=(-0.83, 0.83)
    )
    aperture = fem.shapes.Annulus(center=(0.0, 0.0), inner_radius=0.07, outer_radius=0.16)
    bore = fem.shapes.Extrusion(shape=aperture, z_range=(-0.13, 0.13))
    block = fem.shapes.Box(bounds=((-0.20, 0.20), (-0.20, 0.20), (-0.17, 0.17)))
    solver.add_geometry(
        shape=fem.shapes.Difference(shape=block, tool=bore), material=fem.materials.PEC
    )
    solver.add_geometry(shape=bore, material=fem.materials.Material(name="Fill", epsilon=2.25))
    for name, z in (("in", -0.13), ("out", 0.13)):
        solver.add_waveguide_port(
            name=name, axis="z", position=z, cross_section=aperture, amplitude=float(name == "in")
        )
    solver.mesh(max_element_size=0.11)
    modes = solver.solve_port_modes(port="in")
    assert modes[0].neff.real == pytest.approx(1.5, rel=0.015)
    case = solver.solve(backend="cpu", max_refinements=0)
    assert abs(case.s_parameters[1, 0]) > 0.97
    assert abs(case.s_parameters[1, 0] - np.exp(-2j * np.pi * 1.5 * 0.26)) < 0.15
    loaded = fem.load(case.save(tmp_path / "coax.h5"))
    assert loaded.port_modes["out"].mesh_data is not None
    np.testing.assert_allclose(loaded.port_modes["in"].neff, modes.neff)
    assert loaded.port_modes["in"].plot(component="Ey").axes


def test_discrete_port_power_and_surface_current():
    from benchmarks.analytical.short_dipole import build_dipole

    solver = build_dipole()
    solver.mesh(max_element_size=0.085, element_order=1)
    case = solver.solve(backend="cpu")
    power = (case.port_voltage["feed"] * np.conj(case.port_current["feed"])).real / 2
    assert power == pytest.approx(case.accepted_power)
    assert 0 < power < case.incident_power
    assert case.radiated_power > 0
    triangles = case.mesh_data.coordinates[case.current_facets]
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    assert np.linalg.norm(case.surface_current) > 0
    np.testing.assert_allclose(np.sum(normals * case.surface_current, axis=1), 0, atol=1e-12)


@pytest.mark.parametrize("end_radius", [.035, .025])
def test_cylindrical_gap_terminals_and_axial_voltage(end_radius):
    solver = fem.AntennaSolver3D(frequency=fem.C_0, x_range=(-.6, .6), y_range=(-.6, .6),
                                 z_range=(-.6, .6), pml_layers=2)

    def tube(radius, limits):
        return fem.shapes.Cylinder(center=(0, 0), radius=radius, z_range=limits).rotated(
            angle=90, axis=(0, 1, 0))

    for radius, limits in ((.035, (-.075, -.015)), (end_radius, (.015, .075))):
        solver.add_geometry(shape=tube(radius, limits), material=fem.materials.PEC)
    solver.add_geometry(shape=tube(.035, (0., .015)), material=fem.Material(epsilon=2.25))
    gap = tube(.035, (-.015, .015))
    with pytest.raises(fem.GeometryError, match="axis"):
        solver.add_discrete_port(name="wrong", axis="z", gap=gap)
    solver.add_discrete_port(name="feed", axis="x", gap=gap)
    solver.mesh(max_element_size=.12, element_order=1)
    gap_cells = solver.mesh_data.metadata["port_gap_cells"]["feed"]
    assert set(solver.mesh_data.metadata["element_tags"][gap_cells]) == {0, 3}
    if end_radius < .035:
        with pytest.raises(fem.GeometryError, match="end faces"):
            solver._prepare(fem.C_0)
        return
    basis, _, sections, _ = solver._prepare(fem.C_0)
    # For a uniform Ex=1 V/m the first-order edge integrals are the x displacements.
    x = basis.mesh.p[0]
    edge_integrals = np.diff(x[basis.mesh.edges], axis=0).ravel()
    voltage = (sections["feed"].T @ edge_integrals).item()
    assert voltage == pytest.approx(.03)


def test_mixed_waveguide_and_discrete_ports():
    solver = fem.AntennaSolver3D(
        frequency=fem.C_0, x_range=(-0.73, 0.73), y_range=(-0.73, 0.73), z_range=(-0.73, 0.73)
    )
    outer, inner = (-0.1, 0.1), (-0.04, 0.04)
    cavity = fem.shapes.Box(bounds=(outer, outer, (-0.13, 0.13)))
    shield = fem.shapes.Box(bounds=((-0.13, 0.13), (-0.13, 0.13), (-0.16, 0.16)))
    solver.add_geometry(
        shape=fem.shapes.Difference(shape=shield, tool=cavity), material=fem.materials.PEC
    )
    solver.add_geometry(
        shape=fem.shapes.Box(bounds=(inner, inner, (-0.13, 0.118))), material=fem.materials.PEC
    )
    solver.add_waveguide_port(
        name="guide",
        axis="z",
        cross_section=fem.shapes.Rectangle(bounds=(outer, outer)),
        position=-0.13,
        amplitude=0.5,
        phase_deg=30,
    )
    solver.add_discrete_port(
        name="load", axis="z", gap=fem.shapes.Box(bounds=(inner, inner, (0.118, 0.13)))
    )
    solver.mesh(max_element_size=0.09)
    case = solver.solve(backend="cpu")
    assert case.channels == (("guide", 1), ("load", "gap"))
    assert abs(case.s_parameters[1, 0]) > 0.7
    np.testing.assert_allclose(case.s_parameters, case.s_parameters.T, atol=1e-9)
    assert abs(case.accepted_power) < 1e-8
    np.testing.assert_allclose(
        case.active_s_parameters, case.s_parameters @ case.incident / case.incident
    )


@pytest.mark.skipif(os.environ.get("FEM_TEST_CUDA") != "1", reason="requires CUDA/cuDSS")
def test_cudss_solves_all_port_columns():
    solver = guide()
    reference = solver.solve(backend="cpu", max_refinements=0)
    accelerated = solver.solve(backend="gpu", max_refinements=0)
    assert accelerated.solve_info["method"] == "cudss-direct"
    np.testing.assert_allclose(accelerated.s_parameters, reference.s_parameters, atol=1e-9)
    np.testing.assert_allclose(accelerated.electric_field, reference.electric_field, atol=1e-8)
    np.testing.assert_allclose(accelerated.magnetic_field, reference.magnetic_field, atol=1e-9)
    np.testing.assert_allclose(accelerated.surface_current_vertices, reference.surface_current_vertices, atol=1e-9)


def test_higher_mode_launch_and_passive_multimode_receiver():
    s = guide(higher_mode=True)
    s.set_excitation(port="right", amplitude=0)
    r = s.solve(backend="cpu")
    assert r.channels == (("left", 2), ("right", 1))
    assert abs(r.outgoing["right", 2]) > 0.95
    assert abs(r.outgoing["right", 1]) < 0.1
    assert abs(r.accepted_power) < 1e-8


def test_huygens_transform_against_hertzian_dipole():
    from fem.antenna.farfield import HuygensSurface
    from fem.common.constants import ETA_0

    # Exact fields of a unit z-directed current moment on a sphere.
    k, r = 2 * np.pi, 0.4
    mu, w = np.polynomial.legendre.leggauss(14)
    mu, phi = np.meshgrid(mu, np.arange(28) * np.pi / 14, indexing="ij")
    sin = np.sqrt(1 - mu**2)
    n = np.stack((sin * np.cos(phi), sin * np.sin(phi), mu), axis=-1)
    theta = np.stack((mu * np.cos(phi), mu * np.sin(phi), -sin), axis=-1)
    az = np.stack((-np.sin(phi), np.cos(phi), np.zeros_like(phi)), axis=-1)
    factor = 1j * k * np.exp(-1j * k * r) * sin / (4 * np.pi * r)
    h_phi = factor * (1 + 1 / (1j * k * r))
    e_theta = ETA_0 * factor * (1 + 1 / (1j * k * r) - 1 / (k * r) ** 2)
    dS = r * r * w[:, None] * np.pi / 14
    j = -theta * (h_phi * dS)[:, :, None]
    m = -az * (e_theta * dS)[:, :, None]
    h = HuygensSurface(k, r * n.reshape(-1, 3), j.reshape(-1, 3), m.reshape(-1, 3))
    amplitude = h.amplitude([[1, 0, 0], [0, 0, 1]])
    np.testing.assert_allclose(amplitude[0], [0, 0, -1j * ETA_0 * k / (4 * np.pi)], atol=1e-9)
    np.testing.assert_allclose(amplitude[1], 0, atol=1e-9)
    assert h.radiated_power() == pytest.approx(ETA_0 * k * k / (12 * np.pi), rel=1e-9)
