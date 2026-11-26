"""Matplotlib scalar fields, VTK vector glyphs and bistatic RCS patterns."""
import numpy as np

from fem.common.errors import ConfigurationError
from fem.common._vector_plot import aligned_components, draw_vectors

COMPONENTS = ("incident", "total", "scattered", "RCS")


def field_views(result):
    if result.dimension == 3:
        views = []
        for field in COMPONENTS[:-1]:
            views.append(f"{field.capitalize()} E")
            if field == "incident" or result.scattered_h is not None:
                views.append(f"{field.capitalize()} H")
        views.append("RCS")
        if result.current_facets is not None:
            views.append("Surface current")
        return views
    fields = ("Et", "Hz") if result.polarization == "TE" else ("Ez", "Ht")
    return [f"{name.capitalize()} {field}" for name in COMPONENTS[:-1] for field in fields
            if field.endswith("z") or result.scattered_transverse is not None] + ["RCS"]


def rcs_surface(result):
    """Radial distance and point colour both represent linear RCS in square metres."""
    import pyvista as pv

    surface = pv.Sphere(radius=1., theta_resolution=48, phi_resolution=25)
    directions = surface.points.astype(float)
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    sigma = result.radar_cross_section(directions)
    surface.points = directions * sigma[:, None]
    del surface.point_data["Normals"]
    surface.point_data["direction"] = directions
    surface.point_data["RCS_m2"] = sigma
    surface.set_active_scalars("RCS_m2")
    return surface


def draw_3d(plotter, result, component, quantity, *, mesh_overlay=False, material=False):
    from fem.common.viewer import add_arrows, add_3d_overlays

    if component == "RCS":
        surface = rcs_surface(result)
        if surface["RCS_m2"].any():
            plotter.add_mesh(surface, scalars="RCS_m2", cmap="viridis", lighting=False,
                             scalar_bar_args={"title": "RCS (m^2)", "fmt": "%.2g"})
        else:
            plotter.add_points(np.zeros((1, 3)), point_size=8)
            plotter.add_text("Zero RCS in all directions", position="lower_left")
        plotter.add_text("Bistatic RCS", font_size=12)
        plotter.add_axes()
        return ""
    mesh = result.mesh_data
    if component == "Surface current":
        from fem.common.tet_hcurl import surface_vector_samples
        centers, vectors = surface_vector_samples(mesh.coordinates[result.current_facets],
                                                  result.surface_current_vertices)
        values, label, title = vectors.T, "|Js| (A/m)", component
    else:
        centers = mesh.coordinates[mesh.elements].mean(axis=1)
        field, symbol = component.split()
        name = field.lower() if symbol == "E" else f"{field.lower()}_h"
        values = getattr(result, name).T
        label = "|E| (V/m)" if symbol == "E" else "|H| (A/m)"
        title = component
    if quantity == "magnitude":
        values = aligned_components(values)
    add_arrows(plotter, centers, values, quantity, label=label)
    plotter.add_text(title, font_size=12)
    return add_3d_overlays(plotter, mesh, mesh_overlay=mesh_overlay, material=material)


def plot_result(result, *, component, quantity):
    if isinstance(component, str) and component.lower() == "rcs":
        component = "RCS"
    views = field_views(result)
    longitudinal = "Hz" if result.polarization == "TE" else "Ez"
    if component in COMPONENTS[:-1]:
        component = f"{component.capitalize()} {longitudinal if result.dimension == 2 else 'E'}"
    if component is not None and component not in views:
        raise ConfigurationError(f"Available views: {', '.join(views)}.")
    if component != "RCS" and quantity not in ("real", "imag", "magnitude"):
        raise ConfigurationError("Select real, imag, or magnitude.")
    if result.dimension == 3:
        import pyvista as pv
        plotter = pv.Plotter()
        draw_3d(plotter, result, component or "Total E", quantity)
        return plotter

    from matplotlib.figure import Figure
    choices = ([f"{name.capitalize()} {longitudinal}" for name in COMPONENTS[:-1]] + ["RCS"]
               if component is None else (component,))
    layout = (2, 2) if component is None else (1, 1)
    fig = Figure(figsize=(10, 8) if component is None else (7, 5), layout="constrained")
    for index, name in enumerate(choices, start=1):
        if name == "RCS":
            angles = np.linspace(0, 2*np.pi, 361)
            ax = fig.add_subplot(*layout, index, projection="polar")
            directions = np.column_stack((np.cos(angles), np.sin(angles)))
            ax.plot(angles, result.radar_cross_section(directions))
            ax.set_ylim(bottom=0)
            ax.set_title("Bistatic scattering width (m)")
        else:
            ax = fig.add_subplot(*layout, index)
            field, symbol = name.split()
            unit = "V/m" if symbol.startswith("E") else "A/m"
            title = f"{name} ({unit}) · {quantity}"
            if symbol.endswith("t"):
                mesh = result.mesh_data
                centers = mesh.coordinates[mesh.elements].mean(axis=1)
                vectors = getattr(result, f"{field.lower()}_transverse")
                draw_vectors(ax, centers, vectors.T, mesh.coordinates, vertical_axis="y",
                             quantity=quantity, label=f"|{symbol}| ({unit})")
                ax.set_title(title)
                continue
            transform = {"real": np.real, "imag": np.imag, "magnitude": np.abs}[quantity]
            values = transform(getattr(result, field.lower()))
            limit = np.max(np.abs(values)) or 1.
            x, y = result.mesh_data.coordinates.T
            artist = ax.tripcolor(x, y, result.mesh_data.elements, values, shading="gouraud",
                                  vmin=0 if quantity == "magnitude" else -limit, vmax=limit,
                                  cmap="viridis" if quantity == "magnitude" else "coolwarm")
            ax.set(xlabel="x (m)", ylabel="y (m)", aspect="equal",
                   title=title)
            fig.colorbar(artist, ax=ax)
    return fig
