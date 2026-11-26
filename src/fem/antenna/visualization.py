"""Network plots and VTK radiation patterns for the shared result viewer."""

import numpy as np
from matplotlib.figure import Figure

from fem.common.errors import ConfigurationError

PATTERNS = {"Directivity": "directivity", "Gain": "gain", "Realized gain": "realized_gain"}
VECTOR_FIELDS = ("Total E", "Total H", "Surface current")
COMPONENTS = ("S-parameters", *VECTOR_FIELDS, *PATTERNS)


def _db(values):
    # Keep an exact zero visible at the plotting floor; raw results retain zero.
    return 20 * np.log10(np.maximum(abs(values), 1e-15))


def _network_scale(quantity):
    if quantity in (None, "dB"):
        return _db, "Magnitude (dB)"
    if quantity == "phase":
        return lambda z: np.angle(z, deg=True), "Phase (degrees)"
    raise ConfigurationError("Select dB or phase for network parameters.")


def plot_sweep(sweep, quantity):
    transform, label = _network_scale(quantity)
    fig = Figure(figsize=(8, 5), layout="constrained")
    channels = sweep[0].channels
    panels = 2 if len(channels) > 1 else 1
    network = fig.add_subplot(panels, 1, 1)
    f = sweep.frequencies / 1e9
    receivers = dict.fromkeys(key for case in sweep for key in case.modal_s_parameters)
    for dst, mode in receivers:
        for j, (src, launch_mode) in enumerate(channels):
            samples = np.array(
                [
                    r.modal_s_parameters[dst, mode][j]
                    if (dst, mode) in r.modal_s_parameters
                    else np.nan
                    for r in sweep
                ]
            )
            network.plot(f, transform(samples), ".-", label=f"{dst} {mode} ← {src} {launch_mode}")
    network.set(title="S-parameters", ylabel=label, xlabel="Frequency (GHz)")
    network.legend()
    if panels == 2:
        ax = fig.add_subplot(2, 1, 2)
        for i, (name, _) in enumerate(channels):
            ax.plot(f, transform(sweep.active_s_parameters[:, i]), ".-", label=name)
        ax.set(title="Active reflection", ylabel=label, xlabel="Frequency (GHz)")
        ax.legend()
    return fig


def plot_network(case, quantity):
    transform, label = _network_scale(quantity)
    fig = Figure(figsize=(8, 5), layout="constrained")
    ax = fig.subplots()
    names = [f"{p} {m}" for p, m in case.channels]
    keys, rows = zip(*case.modal_s_parameters.items())
    values = transform(np.stack(rows))
    image = ax.imshow(values, cmap="viridis", aspect="auto")
    ax.set(
        xticks=range(len(names)),
        xticklabels=names,
        xlabel="Excited port mode",
        yticks=range(len(keys)),
        yticklabels=[f"{p} {m}" for p, m in keys],
        ylabel="Outgoing mode",
        title=f"S-parameters at {case.frequency / 1e9:.5g} GHz",
    )
    for (i, j), db in np.ndenumerate(values):
        ax.text(
            j,
            i,
            f"{db:.1f}",
            ha="center",
            va="center",
            color="white" if image.norm(db) < 0.5 else "black",
        )
    fig.colorbar(image, ax=ax, label=label)
    if len(names) > 1:
        labels = []
        for (name, _), val in zip(case.channels, case.active_s_parameters):
            value = (
                f"{_db(val):.2f} dB ∠ {np.angle(val, deg=True):.1f}°"
                if np.isfinite(val)
                else "undriven"
            )
            labels.append(f"{name}  {value}")
        fig.supxlabel("Active reflection\n" + "     ".join(labels))
    return fig


def _phase_vectors(plotter, case, component, quantity, phase_deg):
    """Connect one live glyph pipeline and return its phase control."""
    import pyvista as pv
    from vtkmodules.vtkFiltersCore import vtkGlyph3D
    from fem.common.viewer import _MAX_ARROWS

    if component == "Total E":
        cells, phasors, label = case.mesh_data.elements, case.electric_field, "|E| (V/m)"
    elif component == "Total H":
        if case.magnetic_field is None:
            raise ConfigurationError("Run the antenna solve to obtain the magnetic field.")
        cells, phasors, label = case.mesh_data.elements, case.magnetic_field, "|H| (A/m)"
    else:
        cells, phasors, label = case.current_facets, case.surface_current, "|Js| (A/m)"
    vertices = case.mesh_data.coordinates[cells]
    points = vertices.mean(axis=1)
    if component == "Surface current" and case.surface_current_vertices is not None:
        # Sample inside each face so adjacent sides of a sharp edge retain
        # distinct arrows even at coincident mesh vertices.
        from fem.common.tet_hcurl import surface_vector_samples
        points, phasors = surface_vector_samples(vertices, case.surface_current_vertices)
    step = max(1, int(np.ceil(len(points) / _MAX_ARROWS)))
    points, phasors = points[::step], phasors[::step]
    peak = float(np.max(np.linalg.norm(phasors, axis=1))) or 1.0
    cloud = pv.PolyData(points)
    cloud["vectors"] = np.zeros_like(points)
    cloud[label] = np.zeros(len(points))
    cloud.set_active_vectors("vectors")
    cloud.set_active_scalars(label)
    glyphs = vtkGlyph3D()
    glyphs.SetInputData(cloud)
    glyphs.SetSourceData(pv.Arrow(tip_resolution=8, shaft_resolution=8))
    glyphs.SetVectorModeToUseVector()
    glyphs.SetScaleModeToScaleByVector()
    glyphs.SetColorModeToColorByScalar()
    glyphs.SetScaleFactor(0.08 * np.max(np.ptp(points, axis=0)) / peak)

    def set_vectors(angle):
        field = phasors * np.exp(1j * np.deg2rad(angle))
        shown = field.imag if quantity == "imag" else field.real
        cloud["vectors"][:] = shown
        cloud[label][:] = np.linalg.norm(shown, axis=1)
        cloud.Modified()

    set_vectors(phase_deg)
    actor = plotter.add_mesh(
        glyphs, scalars=label, clim=(0, peak), cmap="viridis", name="antenna_vectors"
    )

    def update_phase(angle):
        set_vectors(angle)
        actor.mapper.update()

    return update_phase


def draw_3d(
    plotter, case, component, quantity, *, phase_deg=0.0, mesh_overlay=False, material=False
):
    import pyvista as pv
    from fem.common.viewer import add_3d_overlays

    if component in VECTOR_FIELDS:
        if component == "Surface current":
            cells = case.current_facets
            faces = np.column_stack((np.full(len(cells), 3), cells))
            skin = pv.PolyData(case.mesh_data.coordinates, faces.ravel()).clean()
            plotter.add_mesh(skin, color="silver", opacity=0.22)
        update_phase = _phase_vectors(plotter, case, component, quantity, phase_deg)
        message = add_3d_overlays(
            plotter, case.mesh_data, mesh_overlay=mesh_overlay, material=material
        )
        return message, update_phase
    if component not in PATTERNS:
        raise ConfigurationError(f"Unknown antenna view {component!r}.")
    sphere = pv.Sphere(radius=1.0, theta_resolution=48, phi_resolution=25)
    directions = np.asarray(sphere.points, dtype=float)
    values = getattr(case.far_field(directions), PATTERNS[component])
    if not np.isfinite(values).all() or not np.any(values):
        plotter.add_text(f"{component}: no measurable power")
    else:
        sphere.points = directions * values[:, None]
        del sphere.point_data["Normals"]
        label = f"{component} (dBi)"
        sphere[label] = 10 * np.log10(np.maximum(values, 1e-15))
        plotter.add_mesh(sphere, scalars=label, lighting=False, cmap="viridis")
    plotter.add_axes()
    return "", None


def plot_result(case, component, quantity, phase_deg):
    if component == "S-parameters":
        return plot_network(case, quantity)
    import pyvista as pv

    plotter = pv.Plotter()
    draw_3d(plotter, case, component, quantity or "real", phase_deg=phase_deg)
    return plotter
