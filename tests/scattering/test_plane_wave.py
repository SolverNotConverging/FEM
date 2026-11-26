"""Physical normalization, field phasing, persistence and visualization regressions."""
import numpy as np
import pytest
from scipy.special import hankel2, h2vp, jv, jvp

import fem
from fem.common.constants import ETA_0


@pytest.mark.parametrize("solver_type", [fem.ScatteringSolver2D, fem.ScatteringSolver3D, fem.AntennaSolver3D])
def test_vacuum_monitor_placement_and_reconfiguration(solver_type):
    dim = 2 if solver_type is fem.ScatteringSolver2D else 3
    shape = fem.shapes.Rectangle if dim == 2 else fem.shapes.Box
    extent = (-1/3, 2 + 1/3)
    settings = dict(frequency=fem.C_0, x_range=extent, y_range=extent, pml_layers=3)
    if dim == 3:
        settings["z_range"] = extent
    solver = solver_type(**settings)
    target = shape(bounds=[(.4, 1.)] * dim)
    solver.add_geometry(shape=target, material=fem.Material(epsilon=2.))
    # The implicit absorber uses the same automatic thickness as add_pml().
    solver.add_nf2ff("auto")
    mesh = solver.mesh(max_element_size=.35)
    np.testing.assert_allclose(np.reshape(mesh.info["nf2ff_bounds"], (-1, 2)), [[.2, 1.7]] * dim)
    with pytest.raises(fem.GeometryError, match="vacuum"):
        solver.add_nf2ff(distance=.4)
    with pytest.raises(fem.GeometryError, match="PML"):
        solver.add_pml(thickness=.8)
    solver.add_nf2ff(distance=.1)
    assert solver.mesh_data is None
    with pytest.raises(fem.GeometryError, match="vacuum"):
        solver.add_geometry(shape=shape(bounds=[(.05, .11)] * dim), material=fem.materials.PEC)
    solver.add_pml(thickness=.24)
    mesh = solver.mesh(max_element_size=.35)
    np.testing.assert_allclose(np.reshape(mesh.info["nf2ff_bounds"], (-1, 2)),
                               [[extent[0] + .34, extent[1] - .34]] * dim)
    np.testing.assert_allclose(solver._mesh.p.min(axis=1), extent[0])
    assert mesh.info["pml_layers"] == 3


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("round_region", [True, False])
def test_local_refinement_concentrates_cells_and_preserves_material_tags(dim, round_region):
    center = np.array([1.4, -0.5, 2.0][:dim])
    options = dict(frequency=fem.C_0, x_range=(.9-1/3, 1.9+1/3), y_range=(-1.-1/3, 1/3))
    if dim == 2:
        solver = fem.ScatteringSolver2D(**options)
    else:
        solver = fem.ScatteringSolver3D(**options, z_range=(1.5-1/3, 2.5+1/3))
    if round_region:
        region = (fem.shapes.Circle if dim == 2 else fem.shapes.Sphere)(center=center, radius=.22)
    else:
        region = (fem.shapes.Rectangle if dim == 2 else fem.shapes.Box)(
            bounds=[(x-.22, x+.22) for x in center])
    solver.mesh(max_element_size=.16)
    solver.add_mesh_refinement(shape=region, max_element_size=.05)
    assert solver.mesh_data is None
    mesh = solver.mesh(max_element_size=.16)
    vertices = mesh.coordinates[mesh.elements]
    distances = np.linalg.norm(vertices.mean(axis=1) - center, axis=1)
    sizes = np.linalg.norm(vertices[:, 1] - vertices[:, 0], axis=1)
    assert np.median(sizes[distances < .17]) < .65 * np.median(sizes[distances > .45])
    assert np.all(mesh.metadata["element_tags"] == 0)


def cylinder(polarization, material=None):
    solver = fem.ScatteringSolver2D(frequency=fem.C_0, x_range=(-1.2, 1.2), y_range=(-1.2, 1.2),
                                    polarization=polarization)
    if material is not None:
        solver.add_geometry(shape=fem.shapes.Circle(center=(0, 0), radius=.18), material=material)
    solver.mesh(max_element_size=.035 if material is not None else .2)
    return solver.solve(backend="cpu")


@pytest.mark.parametrize("polarization", ["TM", "TE"])
@pytest.mark.parametrize("target", ["dielectric", "PEC", "impedance"])
def test_cylinder_rcs_against_partial_wave_series(polarization, target):
    k = 2*np.pi
    x = k * .18
    n = np.arange(-14, 15)
    z = .55 + .2j
    if target == "dielectric":
        material = fem.Material(epsilon=2.25)
        q = 1.5 if polarization == "TM" else 1/1.5
        alpha = q * jvp(n, 1.5*x) / jv(n, 1.5*x)
        coef = -(jvp(n, x) - alpha*jv(n, x)) / (h2vp(n, x) - alpha*hankel2(n, x))
    elif target == "PEC":
        material = fem.materials.PEC
        coef = -jv(n,x)/hankel2(n,x) if polarization == "TM" else -jvp(n,x)/h2vp(n,x)
    else:
        material = fem.SurfaceImpedance(impedance=z * ETA_0)
        alpha = 1j * (z if polarization == "TE" else 1/z)
        coef = -(jvp(n, x) - alpha*jv(n, x)) / (h2vp(n, x) - alpha*hankel2(n, x))
    r = cylinder(polarization, material)
    theta = np.linspace(0, np.pi, 19)
    exact = 4/k * abs(np.exp(1j*theta[:,None]*n) @ coef)**2
    numerical = r.radar_cross_section(np.column_stack((np.cos(theta), np.sin(theta))))
    assert np.linalg.norm(numerical-exact) / np.linalg.norm(exact) < .12


@pytest.mark.parametrize("dim", [2, 3])
def test_vacuum_plane_wave(dim):
    direction = np.ones(dim)/np.sqrt(dim)
    outer = (-1/3, 1 + 1/3)
    options = dict(frequency=fem.C_0, x_range=outer, y_range=outer, direction=direction)
    if dim == 2:
        solver = fem.ScatteringSolver2D(**options)
    else:
        solver = fem.ScatteringSolver3D(**options, z_range=outer, electric_vector=(1,-1,0))
    solver.mesh(max_element_size=.4)
    assert solver.mesh_data.info["pml_layers"] == 5
    assert solver.mesh_data.info["element_order"] == 2
    assert solver.mesh_data.info["pml_thickness"] == pytest.approx(1/3)
    np.testing.assert_allclose(solver._mesh.p.min(axis=1), -1/3)
    np.testing.assert_allclose(solver._mesh.p.max(axis=1), 1 + 1/3)
    r = solver.solve(backend="cpu")
    assert r.solve_info["adaptive_history"][0]["status"] == "tolerance"
    assert len(r.solve_info["adaptive_history"]) == 1
    points = r.mesh_data.coordinates
    if dim == 3:
        points = points[r.mesh_data.elements].mean(axis=1)
    wave = np.exp(-2j*np.pi*(points @ direction))
    expected = wave if dim == 2 else wave[:,None] * solver.electric_vector
    np.testing.assert_allclose(r.total, expected)
    np.testing.assert_array_equal(r.scattered, 0)
    np.testing.assert_array_equal(r.radar_cross_section(-direction), 0)


def test_adaptive_dielectric_mesh_budget_and_orders():
    solver = fem.ScatteringSolver2D(frequency=fem.C_0, x_range=(-1, 1), y_range=(-1, 1))
    solver.add_geometry(shape=fem.shapes.Circle(center=(0, 0), radius=.18),
                        material=fem.Material(epsilon=2.25))
    solver.mesh(max_element_size=.15, element_order=1, quadrature_order=6)
    result = solver.solve(backend="cpu", max_refinements=2, adaptive_tolerance=1e-8)
    history = result.solve_info["adaptive_history"]
    assert [step["status"] for step in history] == ["refine", "refine", "refinement_limit"]
    assert history[-1]["elements"] > history[0]["elements"]
    assert history[-1]["residual"] < history[0]["residual"]
    assert result.mesh_data.info["quadrature_order"] == 6
    assert result.mesh_data.info["element_order"] == 1


def test_scalar_plot_and_roundtrip(tmp_path):
    import matplotlib.pyplot as plt
    r = cylinder("TE", fem.materials.copper)
    path = r.save(tmp_path / "cylinder.h5")
    loaded = fem.load(path)
    np.testing.assert_array_equal(loaded.total, r.total)
    np.testing.assert_array_equal(loaded.total_transverse, r.total_transverse)
    np.testing.assert_array_equal(loaded.radar_cross_section([1,0]), r.radar_cross_section([1,0]))
    figure = loaded.plot(quantity="imag")
    assert len(figure.axes) == 7
    plt.close(figure)


@pytest.mark.parametrize("polarization", ["TE", "TM"])
def test_plane_wave_transverse_field(polarization):
    direction = np.array([.6, .8])
    solver = fem.ScatteringSolver2D(frequency=fem.C_0, polarization=polarization,
                                    x_range=(-1, 1), y_range=(-1, 1), direction=direction)
    case = solver.solve(backend="cpu", max_refinements=0)
    centers = case.mesh_data.coordinates[case.mesh_data.elements].mean(axis=1)
    vector = np.array([-direction[1], direction[0]])
    vector *= ETA_0 if polarization == "TE" else -1 / ETA_0
    expected = np.exp(-2j * np.pi * (centers @ direction))[:, None] * vector
    np.testing.assert_allclose(case.incident_transverse, expected)
    np.testing.assert_allclose(case.total_transverse, expected)
    np.testing.assert_array_equal(case.scattered_transverse, 0)


@pytest.mark.parametrize("polarization", ["TE", "TM"])
def test_transverse_total_uses_local_material(polarization):
    solver = fem.ScatteringSolver2D(frequency=fem.C_0, polarization=polarization,
                                    x_range=(-1, 1), y_range=(-1, 1))
    material = fem.Material(epsilon=2.25, mu=1.4)
    solver.add_geometry(shape=fem.shapes.Rectangle(bounds=((-0.2, 0.2), (-0.1, 0.1))),
                        material=material)
    solver.mesh(max_element_size=.13, element_order=1)
    case = solver.solve(backend="cpu", max_refinements=0)
    mesh = case.mesh_data
    inside = mesh.metadata["element_tags"] == 1
    triangles = mesh.elements[inside]
    delta = mesh.coordinates[triangles[:, 1:]] - mesh.coordinates[triangles[:, :1]]
    differences = case.scattered[triangles[:, 1:]] - case.scattered[triangles[:, :1]]
    gradient = np.linalg.solve(delta, differences[..., None])[..., 0]
    rotation = np.column_stack((gradient[:, 1], -gradient[:, 0]))
    factor, relative = (-1j * ETA_0, 2.25) if polarization == "TE" else (1j / ETA_0, 1.4)
    expected = (case.incident_transverse[inside] + factor / (2*np.pi) * rotation) / relative
    np.testing.assert_allclose(case.total_transverse[inside], expected)


def test_dielectric_sphere_rayleigh_rcs_and_vtk(tmp_path):
    import pyvista as pv
    k, radius, eps = 2*np.pi, .035, 2.
    extent = (-.28-1/3, .28+1/3)
    solver = fem.ScatteringSolver3D(frequency=fem.C_0, x_range=extent, y_range=extent,
                                    z_range=extent)
    solver.add_geometry(shape=fem.shapes.Sphere(center=(0,0,0), radius=radius),
                        material=fem.Material(epsilon=eps))
    solver.mesh(max_element_size=.12, element_order=1)
    r = solver.solve(backend="cpu")
    expected = 4*np.pi*k**4*radius**6*((eps-1)/(eps+2))**2
    back, null = r.radar_cross_section([[-1,0,0], [0,1,0]])
    assert back == pytest.approx(expected, rel=.12)
    assert null < back * .002
    path = fem.load(r.save(tmp_path/"sphere.h5")).export_vtk(tmp_path/"sphere.vtm")
    dataset = pv.read(path)
    assert dataset.keys() == ["fields", "rcs"]
    for name in ("incident", "total", "scattered"):
        vector = dataset["fields"].cell_data[name+"_real"] + 1j*dataset["fields"].cell_data[name+"_imag"]
        np.testing.assert_allclose(vector, getattr(r, name))
    sphere = dataset["rcs"]
    np.testing.assert_allclose(sphere["RCS_m2"], r.radar_cross_section(sphere["direction"]))
    np.testing.assert_allclose(sphere.points, sphere["direction"] * sphere["RCS_m2"][:, None])


@pytest.mark.parametrize(("dim", "polarization"), [(2, "TE"), (2, "TM"), (3, None)])
def test_shared_viewer_switches_fields_and_rcs(dim, polarization, tmp_path):
    outer = (-1/3, 1 + 1/3)
    options = dict(frequency=fem.C_0, x_range=outer, y_range=outer)
    if dim == 3:
        solver = fem.ScatteringSolver3D(**options, z_range=outer, pml_layers=2)
        target = fem.shapes.Box(bounds=((.35, .65),)*3)
    else:
        solver = fem.ScatteringSolver2D(**options, polarization=polarization)
        target = fem.shapes.Rectangle(bounds=((.35, .65),)*2)
    solver.add_geometry(shape=target, material=fem.materials.copper)
    solver.mesh(max_element_size=.23)
    result = solver.solve(backend="cpu", max_refinements=0)
    viewer = fem.show(result.save(tmp_path/"case.h5"), block=False)
    try:
        symbols = ("Et", "Hz") if polarization == "TE" else ("Ez", "Ht")
        assert viewer._component.currentText() == ("total" if dim == 3 else f"Total {symbols[0]}")
        viewer._material.setChecked(True)
        fields = ("incident", "scattered", "total", "RCS")
        if dim == 3:
            fields += ("Surface current",)
            skin = result.vtk_datasets()["surface_current"]
            np.testing.assert_allclose(skin["Js_real"] + 1j * skin["Js_imag"],
                                       result.surface_current_vertices.reshape(-1, 3))
        else:
            fields = [f"{name} {symbol}" for name in ("Incident", "Total", "Scattered")
                      for symbol in symbols] + ["RCS"]
            assert fields == [viewer._component.itemText(i) for i in range(viewer._component.count())]
        for component in fields:
            viewer._component.setCurrentText(component)
            assert viewer._status.currentMessage() == ""
            if component == "RCS":
                assert not viewer._quantity.isEnabled() and not viewer._mesh.isEnabled()
            else:
                viewer._quantity.setCurrentText("imag")
                assert viewer._quantity.isEnabled() and viewer._mesh.isEnabled()
            if dim == 3:
                datasets = [actor.mapper.dataset for actor in viewer._plotter.actors.values()
                            if hasattr(getattr(actor, "mapper", None), "dataset")]
                if component == "RCS":
                    assert any("RCS_m2" in data.point_data for data in datasets)
                else:
                    assert "boundary_sibc" in viewer._plotter.actors
                    label = "|Js| (A/m)" if component == "Surface current" else "|E| (V/m)"
                    assert any(label in data.point_data for data in datasets)
            else:
                from matplotlib.quiver import Quiver
                ax = viewer._canvas.figure.axes[0]
                assert (ax.name == "polar") == (component == "RCS")
                assert any(isinstance(c, Quiver) for c in ax.collections) == component.endswith("t")
    finally:
        viewer._clear_body()
        viewer._window.close()
