"""Data-only electrostatic archives and Matplotlib inspection."""
from dataclasses import fields
from pathlib import Path
import h5py
import numpy as np

from fem.common import MeshSnapshot
from fem.common.contracts import mesh_snapshot
from fem.common.errors import PersistenceError, ConfigurationError
from fem.common.persistence import atomic_h5, write_envelope, validate_envelope, write_value, read_value
from .results import ElectrostaticResult


def save_result(result, path):
    mesh = result.mesh_data
    payload = {item.name: getattr(result, item.name) for item in fields(result) if item.name != "mesh"}
    with atomic_h5(path) as handle:
        write_envelope(handle, family="electrostatics", kind="potential", dimension=len(mesh.axes),
                       representation="nodal-potential; nodal-and-cell-fields", static=True)
        write_value(handle, "mesh", mesh)
        write_value(handle, "result", payload)
    return Path(path)


def load_result(path):
    try:
        with h5py.File(path, "r") as handle:
            validate_envelope(handle, family="electrostatics", static=True)
            registry = {"MeshSnapshot": MeshSnapshot}
            mesh = read_value(handle["mesh"], registry)
            return ElectrostaticResult(mesh=mesh, **read_value(handle["result"], registry))
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise PersistenceError(f"Cannot load electrostatic result from {path}: {exc}") from exc


def _draw(ax, result, component, quantity, *, colorbar_axes=None):
    import matplotlib.tri as mtri
    points, cells = result.coordinates, result.elements
    ax.set_xlim(points[:, 0].min(), points[:, 0].max())
    ax.set_xlabel("x (m)")
    if points.shape[1] == 2:
        ax.set_ylim(points[:, 1].min(), points[:, 1].max())
        ax.set_aspect("equal", adjustable="box")
        ax.set_ylabel("y (m)")
    if component == "potential":
        value, unit, location = result.potential, "V", "node"
        if getattr(result, "potential_coefficients", None) is not None:
            from skfem import Basis, ElementLineP2, ElementTriP2, MeshLine, MeshTri
            line = points.shape[1] == 1
            native = (MeshLine if line else MeshTri)(points.T, cells.T)
            basis = Basis(native, ElementLineP2() if line else ElementTriP2())
            sampled, value = basis.refinterp(result.potential_coefficients, nrefs=1)
            points, cells = sampled.p.T, sampled.t.T
    elif component == "mesh":
        if points.shape[1] == 1:
            ax.plot(points[:, 0], np.zeros(len(points)), "|-", color="black")
        else:
            ax.triplot(points[:, 0], points[:, 1], cells, color="black", linewidth=.5)
        ax.set_title("FEM mesh")
        return
    elif component in ("Ex", "Ey", "E", "Et", "Dx", "Dy", "D", "Dt"):
        field = result.element_electric_field if component.startswith("E") else result.element_displacement_field
        if field is None:
            raise ConfigurationError("The result has no cell fields.")
        index = {"x": 0, "y": 1}.get(component[-1])
        if index is not None and index >= points.shape[1]:
            raise ConfigurationError(f"{component} is not defined in one dimension.")
        unit, location = ("V/m" if component.startswith("E") else "C/m²"), "cell"
        if component in ("Et", "Dt"):
            from fem.common._vector_plot import draw_vectors
            if quantity not in ("real", "magnitude", "abs"):
                raise ConfigurationError("Static vector fields support real or magnitude.")
            vectors = field.T if field.shape[1] == 2 else np.vstack((field[:, 0], np.zeros(len(field))))
            draw_vectors(ax, points[cells].mean(axis=1), vectors, points, vertical_axis="y",
                         quantity=quantity, label=f"|{component}| ({unit})", colorbar_axes=colorbar_axes)
            if points.shape[1] == 2:
                ax.set_xlim(points[:, 0].min(), points[:, 0].max())
                ax.set_ylim(points[:, 1].min(), points[:, 1].max())
            ax.set_title(component)
            return
        value = np.linalg.norm(field, axis=1) if index is None else field[:, index]
    else:
        raise ConfigurationError("component must be potential, mesh, E, Ex, Ey, Et, D, Dx, Dy, or Dt.")
    if quantity in ("magnitude", "abs"): value = np.abs(value)
    elif quantity != "real": raise ConfigurationError("Static fields support quantity='real' or 'magnitude'.")
    if points.shape[1] == 1:
        x = points[:, 0] if location == "node" else points[cells, 0].mean(axis=1)
        order = np.argsort(x, kind="stable")
        ax.plot(x[order], value[order])
        ax.set_ylabel(unit)
    else:
        triangle = mtri.Triangulation(points[:, 0], points[:, 1], cells)
        artist = ax.tripcolor(triangle, value, shading="gouraud") if location == "node" else ax.tripcolor(triangle, facecolors=value, shading="flat")
        ax.figure.colorbar(artist, ax=ax, cax=colorbar_axes, label=unit)
        ax.set_aspect("equal")
        ax.set_ylabel("y (m)")
    ax.set_xlabel("x (m)")
    ax.set_title(component)


def draw_3d(plotter, result, component, quantity, *, mesh_overlay=False, material=False):
    """Display field arrows at tetrahedron centres or scalars on three midplanes."""
    from fem.common.viewer import _tetra_grid, add_3d_overlays, add_arrows

    if quantity not in ("real", "magnitude", "abs"):
        raise ConfigurationError("Static fields support real or magnitude.")
    mesh = result.mesh_data
    if component == "mesh":
        return add_3d_overlays(plotter, mesh, mesh_overlay=True, material=material)
    grid = _tetra_grid(mesh)
    if component == "potential":
        values, label = result.potential, "Potential (V)"
        coefficients = result.potential_coefficients
        if coefficients is not None:
            from skfem import Basis, ElementTetP2, MeshTet
            import pyvista as pv
            native = MeshTet(mesh.coordinates.T, mesh.elements.T)
            basis = Basis(native, ElementTetP2())
            # Both element conventions order the six edge DOFs as 01, 12, 20, 03, 13, 23.
            cells = np.column_stack((np.full(native.nelements, 10), basis.element_dofs.T))
            grid = pv.UnstructuredGrid(cells.ravel(), np.full(native.nelements, pv.CellType.QUADRATIC_TETRA),
                                       basis.doflocs.T)
            values = coefficients
        grid.point_data[label] = np.abs(values) if quantity != "real" else values
    elif component in ("E", "D", "Ex", "Ey", "Ez", "Dx", "Dy", "Dz", "|E|", "|D|"):
        symbol = component.strip("|")[0]
        field, unit = ((result.element_electric_field, "V/m") if symbol == "E"
                       else (result.element_displacement_field, "C/m^2"))
        if field is None:
            raise ConfigurationError(f"This result has no cell {symbol} field.")
        label = f"{component} ({unit})"
        if component == symbol:
            centers = mesh.coordinates[mesh.elements].mean(axis=1)
            add_arrows(plotter, centers, field.T, quantity, label=f"|{symbol}| ({unit})")
        else:
            values = (np.linalg.norm(field, axis=1) if component.startswith("|")
                      else field[:, "xyz".index(component[-1])])
            grid.cell_data[label] = np.abs(values) if quantity != "real" else values
    else:
        raise ConfigurationError(f"Unknown 3D electrostatic field {component!r}.")
    if component not in ("E", "D"):
        slices = grid.slice_orthogonal()
        plotter.add_mesh(slices.combine(), scalars=label, cmap="viridis")
    plotter.add_mesh(grid.outline(), color="grey")
    return add_3d_overlays(plotter, mesh, mesh_overlay=mesh_overlay, material=material)


def plot_result(result, *, component=None, quantity="real"):
    if result.coordinates.shape[1] == 3:
        import pyvista as pv
        plotter = pv.Plotter()
        try:
            draw_3d(plotter, result, component or "potential", quantity)
        except Exception:
            plotter.close()
            raise
        return plotter
    from matplotlib.figure import Figure
    figure = Figure()
    ax = figure.subplots()
    _draw(ax, result, component or "potential", quantity)
    figure.tight_layout()
    return figure


def show_result(result, *, block=True):
    if result.coordinates.shape[1] == 3:
        from fem.common.viewer import show
        return show(result, block=block)
    import matplotlib.pyplot as plt
    from matplotlib.widgets import RadioButtons
    figure = plt.figure(figsize=(9, 6))
    selector_axes = figure.add_axes((.02, .25, .15, .5))
    choices = ("potential", "Ex", *(("Ey",) if result.coordinates.shape[1] == 2 else ()),
               "E", "D", "mesh", "Et", "Dt")
    buttons = RadioButtons(selector_axes, choices)
    def redraw(component):
        for axis in list(figure.axes):
            if axis is not selector_axes:
                figure.delaxes(axis)
        is_2d = result.coordinates.shape[1] == 2
        axis = figure.add_axes((.25, .15, .57 if is_2d else .65, .75))
        colorbar_axes = (figure.add_axes((.85, .15, .025, .75))
                         if is_2d and component != "mesh" else None)
        _draw(axis, result, component, "real", colorbar_axes=colorbar_axes)
        figure.canvas.draw_idle()
    buttons.on_clicked(redraw)
    figure._fem_selector = buttons
    redraw("potential")
    plt.show(block=block)
    return figure
