"""Interactive Qt/PyVista result viewer shared by every solver family.

One window shows any result: choose the case, mode, field component and
quantity, and toggle the mesh and the material regions. Et/Ht draw transverse
arrows. One- and two-dimensional fields use Matplotlib. Three-dimensional
fields use PyVista vector glyphs, and scattering RCS uses a radial surface.
The same controls operate on results loaded from HDF5.
"""
from __future__ import annotations

import os
from typing import Any

import numpy as np

from ._vector_plot import aligned_components, draw_mode_vectors

_WAVE_COMPONENTS = ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "E", "H")
_QUANTITIES = ("real", "imag", "magnitude", "phase")
_STATIC_QUANTITIES = ("real", "magnitude")
_MAX_ARROWS = 2500
_TRANSVERSE_COMPONENTS = ('Et', 'Ht')


def _relative_permittivity(description: Any) -> complex:
    """Return the leading relative permittivity of a recorded material description."""
    if not isinstance(description, dict):
        return complex("nan")
    for key in ("eps_r", "epsilon"):
        if key in description:
            value = np.asarray(description[key], dtype=complex).ravel()
            return complex(value[0]) if value.size else complex("nan")
    if "permittivity" in description:
        return _relative_permittivity(description["permittivity"])
    if "tensor" in description:
        return complex(np.asarray(description["tensor"], dtype=complex).ravel()[0])
    if "material" in description:
        return _relative_permittivity(description["material"])
    return complex("nan")


def cell_materials(mesh_data: Any) -> tuple[np.ndarray, list[str]] | None:
    """Return a region index per mesh cell and one label per region.

    Uses the cell tags, region names and configured materials that every solver
    records in ``mesh_data``. Returns None when the mesh carries no tags.
    """
    if mesh_data is None:
        return None
    metadata = mesh_data.metadata
    tags = metadata.get("element_tags")
    if tags is None or len(tags) != len(mesh_data.elements):
        return None
    tags = np.asarray(tags)
    names = dict(metadata.get("physical_names") or {})
    context = metadata.get("context") or {}
    background = context.get("background_material", context.get("background"))
    materials = {"background": background, "exterior": background}
    for key in ("regions", "materials"):
        for region in context.get(key) or ():
            if isinstance(region, dict) and "name" in region:
                materials[region["name"]] = region.get("material", region.get("permittivity", region))
    unique = sorted(set(tags.tolist()))
    labels = []
    for tag in unique:
        name = names.get(tag, f"region {tag}")
        epsilon = _relative_permittivity(materials.get(name))
        labels.append(name if np.isnan(epsilon.real) else f"{name} (εr={epsilon.real:.4g})")
    index = np.searchsorted(np.asarray(unique), tags)
    return index, labels


class _Source:
    """Uniform view of a result: named cases, each optionally with modes."""

    def __init__(self, result: Any):
        self.result = result
        name = type(result).__name__
        if name in ("PeriodicSweepResult", "AntennaSweep"):
            self.cases = list(result)
        else:
            self.cases = [result]
        first = self.cases[0]
        self.has_modes = hasattr(first, "mode") and hasattr(first, "__len__")
        self.family = type(first).__module__.rsplit(".", 1)[0]
        self.static = self.family == "fem.electrostatics"
        self.scattering = self.family == "fem.scattering"
        self.antenna = self.family == "fem.antenna"
        self.is_3d = getattr(first, "dimension", None) == 3

    def case_labels(self) -> list[str]:
        labels = []
        frequencies = getattr(self.result, "frequencies", None)
        for index in range(len(self.cases)):
            if frequencies is not None:
                labels.append(f"Case {index}: {float(frequencies[index]) / 1e9:.6g} GHz")
            else:
                labels.append(f"Case {index}")
        return labels if len(labels) > 1 else ["Result"]

    def components(self, case: Any) -> list[str]:
        if self.antenna:
            from fem.antenna.visualization import COMPONENTS
            components = list(COMPONENTS)
            if case.magnetic_field is None:
                components.remove("Total H")
            return [*components, *(["Frequency sweep"] if hasattr(self.result, "cases") else [])]
        if self.scattering:
            from fem.scattering.visualization import field_views
            return field_views(case)
        if self.static:
            two_d = case.coordinates.shape[1] == 2
            return ["potential", "Ex", *(["Ey"] if two_d else []), "E",
                    "Et", "Dx", *(["Dy"] if two_d else []), "D", "Dt"]
        fields = case.mode(0).fields if self.has_modes else getattr(case, "fields", None)
        if fields is not None and hasattr(fields, "components"):
            names = list(fields.components)
            vectors = ('E', 'H', *_TRANSVERSE_COMPONENTS)
            return names + [n for n in vectors if n not in names]
        return [*_WAVE_COMPONENTS, *_TRANSVERSE_COMPONENTS]

    def mesh_data(self, case: Any) -> Any:
        mesh = getattr(case, "mesh_data", None)
        return mesh[0] if isinstance(mesh, tuple) else mesh

    def info(self, case: Any, mode: int) -> str:
        if self.antenna:
            return (f"{case.frequency/1e9:.6g} GHz   ·   incident {case.incident_power:.4g} W   ·   "
                    f"accepted {case.accepted_power:.4g} W   ·   radiated {case.radiated_power:.4g} W")
        if self.scattering:
            return f"{case.frequency / 1e9:.6g} GHz   ·   {case.polarization}   ·   {case.solve_info['method']}"
        if self.static:
            return f"energy = {case.energy:.6g} J   (static field: E and D; no H)"
        if self.has_modes or hasattr(case, "neff"):
            selected = case.mode(mode) if self.has_modes else case
            neff = complex(selected.neff)
            return f"neff = {neff.real:.6g} {neff.imag:+.3g}j"
        return ""


def _application():
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance()
    return app if app is not None else QApplication([])


def _display_points(source: _Source, points: np.ndarray) -> np.ndarray:
    """Return coordinates as a (nodes, dimension) array."""
    points = np.asarray(points, dtype=float)
    if points.ndim == 1:
        points = points[:, None]
    return points


def _draw_mesh(ax: Any, source: _Source, mesh: Any) -> None:
    points = _display_points(source, mesh.coordinates)
    cells = np.asarray(mesh.elements)
    if points.shape[1] == 1:
        ax.vlines(points[:, 0], 0, 1, transform=ax.get_xaxis_transform(), color="black",
                  linewidth=0.4, alpha=0.5)
    else:
        ax.triplot(points[:, 0], points[:, 1], cells[:, :3], color="black", linewidth=0.3, alpha=0.6)


def _material_figure(source: _Source, mesh: Any, regions: tuple[np.ndarray, list[str]], case: Any = None):
    from matplotlib.colors import ListedColormap
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch

    index, labels = regions
    colours = _region_colours(labels)
    figure = Figure()
    ax = figure.subplots()
    points = _display_points(source, mesh.coordinates)
    cells = np.asarray(mesh.elements)
    if points.shape[1] == 1:
        for cell, region in zip(cells, index):
            ax.axvspan(points[cell, 0].min(), points[cell, 0].max(), color=colours[region], linewidth=0)
        ax.set_yticks([])
    else:
        ax.tripcolor(points[:, 0], points[:, 1], cells[:, :3], facecolors=index.astype(float),
                     cmap=ListedColormap(colours), vmin=-0.5, vmax=len(labels) - 0.5)
        ax.set_aspect("equal", adjustable="box")
    axes = [f"{name} (m)" for name in mesh.axes]
    ax.set_xlabel(axes[0])
    if len(axes) > 1:
        ax.set_ylabel(axes[1])
    handles = [Patch(color=colours[i], label=label) for i, label in enumerate(labels)]
    ax.legend(handles=handles, loc="upper right", fontsize="small", framealpha=0.85)
    ax.set_title("Material regions")
    return figure


class ResultViewer:
    """Window with case/mode/component/quantity selectors and mesh/material toggles."""

    def __init__(self, result: Any, *, title: str | None = None):
        from PySide6.QtCore import QTimer
        from PySide6.QtWidgets import (QPushButton, QCheckBox, QComboBox, QHBoxLayout, QLabel, QMainWindow,
                                       QSpinBox, QVBoxLayout, QWidget)

        self._app = _application()
        self._source = _Source(result)
        self._window = QMainWindow()
        self._window.setWindowTitle(title or f"FEM viewer — {type(result).__name__}")
        self._window.resize(1100, 760)
        self._canvas: Any = None
        self._plotter: Any = None
        self._update_phase = None

        central = QWidget()
        self._layout = QVBoxLayout(central)
        bar = QHBoxLayout()
        self._case = QComboBox()
        self._case.addItems(self._source.case_labels())
        self._mode = QSpinBox()
        self._component = QComboBox()
        self._quantity = QComboBox()
        self._mesh = QCheckBox("Mesh")
        self._material = QCheckBox("Material")
        self._info = QLabel()
        for label, widget in (("Case", self._case), ("Mode", self._mode),
                              ("Component", self._component), ("Quantity", self._quantity)):
            bar.addWidget(QLabel(label))
            bar.addWidget(widget)
        bar.addWidget(self._mesh)
        bar.addWidget(self._material)
        bar.addStretch(1)
        self._layout.addLayout(bar)
        if self._source.antenna:
            controls = QHBoxLayout()
            self._phase = QSpinBox()
            self._phase.setRange(0, 359)
            self._phase.setSuffix("°")
            self._phase.setWrapping(True)
            self._phase.setPrefix("Phase ")
            self._play = QPushButton("Animate phase")
            self._play.setCheckable(True)
            self._timer = QTimer(self._window)
            self._timer.setInterval(80)
            self._timer.timeout.connect(self._advance_phase)
            self._play.toggled.connect(self._toggle_animation)
            self._phase.valueChanged.connect(self._draw_phase)
            controls.addWidget(self._phase)
            controls.addWidget(self._play)
            self._port = QComboBox()
            self._inspect = QPushButton("View port modes")
            self._port_windows = []
            self._inspect.clicked.connect(self._open_port)
            controls.addWidget(self._port)
            controls.addWidget(self._inspect)
            controls.addStretch()
            self._layout.addLayout(controls)
        self._body = QVBoxLayout()
        self._layout.addLayout(self._body, 1)
        self._layout.addWidget(self._info)
        self._window.setCentralWidget(central)
        self._status = self._window.statusBar()

        self._case.currentIndexChanged.connect(self._case_changed)
        self._mode.valueChanged.connect(self._render)
        self._quantity.currentIndexChanged.connect(self._render)
        self._component.currentIndexChanged.connect(self._component_changed)
        self._mesh.toggled.connect(self._render)
        self._material.toggled.connect(self._render)
        self._case_changed()

    # -- state -------------------------------------------------------------
    def _current_case(self) -> Any:
        return self._source.cases[max(self._case.currentIndex(), 0)]

    def _quantities(self) -> list[str]:
        if self._source.antenna:
            view = self._component.currentText()
            from fem.antenna.visualization import VECTOR_FIELDS
            if view in VECTOR_FIELDS:
                return ["real", "imag"]
            return ["dB", "phase"] if view in ("S-parameters", "Frequency sweep") else ["dBi"]
        if self._source.scattering:
            if self._component.currentText() == "RCS":
                return ["m²" if self._source.is_3d else "m"]
            return ["real", "imag", "magnitude"]
        if not self._source.static and self._component.currentText() in _TRANSVERSE_COMPONENTS:
            return ['real', 'imag', 'magnitude']
        return list(_STATIC_QUANTITIES if self._source.static else _QUANTITIES)

    def _set_quantities(self) -> None:
        want = self._quantities()
        if [self._quantity.itemText(i) for i in range(self._quantity.count())] != want:
            previous = self._quantity.currentText()
            self._quantity.blockSignals(True)
            self._quantity.clear()
            self._quantity.addItems(want)
            if previous in want:
                self._quantity.setCurrentText(previous)
            self._quantity.blockSignals(False)

    def _case_changed(self) -> None:
        for widget in (self._mode, self._component):
            widget.blockSignals(True)
        try:
            previous = self._component.currentText()
            self._component.clear()
            case = self._current_case()
            self._component.addItems(self._source.components(case))
            offset = int("port_mesh" in getattr(case, "metadata", {}))
            self._mode.setRange(offset, len(case) - 1 + offset if self._source.has_modes else 0)
            if self._source.antenna:
                previous_port = self._port.currentText()
                self._port.clear()
                self._port.addItems(list(case.port_modes))
                self._port.setCurrentText(previous_port)
                self._inspect.setVisible(bool(case.port_modes))
                self._port.setVisible(bool(case.port_modes))
            default = self._component.findText("potential" if self._source.static else "E")
            if self._source.scattering:
                default = next(i for i in range(self._component.count())
                               if self._component.itemText(i).lower().startswith("total"))
            if self._source.antenna:
                default = self._component.findText(previous or "S-parameters")
            if default >= 0:
                self._component.setCurrentIndex(default)
        finally:
            for widget in (self._mode, self._component):
                widget.blockSignals(False)
        self._component_changed()

    def _component_changed(self) -> None:
        self._set_quantities()
        self._mode.setEnabled(self._source.has_modes)
        spatial = self._component.currentText() != "RCS"
        if self._source.antenna:
            from fem.antenna.visualization import VECTOR_FIELDS
            spatial = self._component.currentText() in VECTOR_FIELDS
            self._play.setVisible(spatial)
            self._phase.setVisible(spatial)
            if not spatial:
                self._play.setChecked(False)
        self._mesh.setEnabled(spatial)
        self._material.setEnabled(spatial)
        self._quantity.setEnabled(self._quantity.count() > 1)
        self._render()

    # -- drawing -----------------------------------------------------------
    def _toggle_animation(self, active):
        if active:
            self._timer.start()
        else:
            self._timer.stop()

    def _advance_phase(self):
        if not self._window.isVisible():
            self._play.setChecked(False)
            return
        self._phase.setValue((self._phase.value() + 10) % 360)

    def _draw_phase(self):
        if self._update_phase is not None:
            self._update_phase(self._phase.value())
            self._plotter.render()
            self._plotter.interactor.update()

    def _open_port(self):
        self._port_windows = [w for w in self._port_windows if w._window.isVisible()]
        self._port_windows.append(self._current_case().show_port_modes(port=self._port.currentText(),
                                                                      block=False))

    def _clear_body(self) -> None:
        self._update_phase = None
        if self._plotter is not None:
            self._plotter.close()
            self._plotter = None
        while self._body.count():
            item = self._body.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()
        self._canvas = None

    def _render(self) -> None:
        component = self._component.currentText()
        quantity = self._quantity.currentText()
        mode = self._mode.value() - self._mode.minimum()
        self._clear_body()
        message = ""
        try:
            case = self._current_case()
            if self._source.antenna and component in ("S-parameters", "Frequency sweep"):
                result = case if component == "S-parameters" else self._source.result
                widget = self._matplotlib(result.plot(quantity=quantity))
            elif self._source.is_3d:
                widget, message = self._pyvista(case, mode, component, quantity)
            else:
                widget, message = self._field_2d(case, mode, component, quantity)
            self._info.setText(self._source.info(case, mode))
            if message:
                self._status.showMessage(message)
            else:
                self._status.clearMessage()
        except Exception as exc:  # keep the window usable when a combination is invalid
            from PySide6.QtWidgets import QLabel
            widget = QLabel(f"Cannot display this selection: {exc}")
            self._status.showMessage(str(exc))
        self._body.addWidget(widget)

    def _field_figure(self, case: Any, mode: int, component: str, quantity: str):
        family = self._source.family
        if family == "fem.waveguide":
            from fem.waveguide.visualization import visualize

            figure, _ = visualize(case, mode=mode, component=component, quantity=quantity,
                                  material=False, show=False)
            return figure
        if family == "fem.periodic":
            from matplotlib.figure import Figure

            from fem.periodic.visualization import visualize

            figure = Figure()
            visualize(case.mode(mode), component=component, quantity=quantity,
                      ax=figure.subplots(), show=False)
            return figure
        return case.plot(component=component, quantity=quantity)

    def _field_2d(self, case: Any, mode: int, component: str, quantity: str):
        mesh = self._source.mesh_data(case)
        if self._source.scattering:
            figure = case.plot(component=component, quantity=quantity)
            message = ""
            if component != "RCS":
                message = _transverse_overlays(figure.axes[0], self._source, mesh,
                                               self._mesh.isChecked(), self._material.isChecked())
                if self._material.isChecked():
                    from matplotlib.collections import LineCollection
                    boundaries = mesh.metadata.get("boundary_facets", {})
                    for kind, edges in boundaries.items():
                        figure.axes[0].add_collection(LineCollection(mesh.coordinates[edges],
                            colors="#b87333" if kind == "sibc" else "black", linewidths=1.8))
            return self._matplotlib(figure), message
        if component in (*_TRANSVERSE_COMPONENTS, "Dt"):
            if self._source.static:
                figure = case.plot(component=component, quantity=quantity)
            else:
                from matplotlib.figure import Figure
                figure = Figure()
                selected = case.mode(mode) if self._source.has_modes else case
                axes = "xz" if self._source.family == "fem.periodic" else "xy"
                draw_mode_vectors(figure.subplots(), selected, component, quantity, axes=axes)
            message = _transverse_overlays(figure.axes[0], self._source, mesh,
                                          self._mesh.isChecked(), self._material.isChecked())
            return self._matplotlib(figure), message
        message = ""
        figure = None
        if self._material.isChecked():
            regions = cell_materials(mesh)
            if regions is None:
                message = "This result records no material regions."
            else:
                figure = _material_figure(self._source, mesh, regions, case)
        if figure is None:
            figure = self._field_figure(case, mode, component, quantity)
        if self._mesh.isChecked():
            if mesh is None:
                message = "This result records no mesh."
            else:
                _draw_mesh(figure.axes[0], self._source, mesh)
        return self._matplotlib(figure), message

    def _matplotlib(self, figure: Any):
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
        from PySide6.QtWidgets import QVBoxLayout, QWidget

        canvas = FigureCanvasQTAgg(figure)
        holder = QWidget()
        layout = QVBoxLayout(holder)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(canvas, holder))
        layout.addWidget(canvas)
        self._canvas = canvas
        return holder

    def _pyvista(self, case: Any, mode: Any, component: str, quantity: str):
        from pyvistaqt import QtInteractor

        plotter = QtInteractor(self._window)
        overlays = dict(mesh_overlay=self._mesh.isChecked(), material=self._material.isChecked())
        if self._source.antenna:
            from fem.antenna.visualization import draw_3d
            message, self._update_phase = draw_3d(plotter, case, component, quantity,
                                                  phase_deg=self._phase.value(), **overlays)
        elif self._source.scattering:
            from fem.scattering.visualization import draw_3d
            message = draw_3d(plotter, case, component, quantity, **overlays)
        else:
            message = populate_3d(plotter, case.mode(mode), self._source.mesh_data(case), component,
                                  quantity, boundary_facets=case.metadata.get("boundary_facets"),
                                  **overlays)
        self._plotter = plotter
        return plotter.interactor, message

    # -- public ------------------------------------------------------------
    def show(self, *, block: bool = True) -> ResultViewer:
        self._window.show()
        if block:
            self._app.exec()
        return self


def _transverse_overlays(ax, source, mesh, mesh_overlay, material):
    """Keep transverse arrows visible when mesh or material is enabled."""
    if not (mesh_overlay or material):
        return ''
    if mesh is None:
        return 'This result records no mesh.'
    points = _display_points(source, mesh.coordinates)
    cells = np.asarray(mesh.elements)
    regions = cell_materials(mesh) if material else None
    message = 'This result records no material regions.' if material and regions is None else ''
    colours = _region_colours(regions[1]) if regions is not None else None
    if points.shape[0] and cells.shape[1] == 2:
        for i, cell in enumerate(cells):
            if regions is not None:
                ax.axvspan(points[cell, 0].min(), points[cell, 0].max(),
                           color=colours[regions[0][i]], alpha=.25, zorder=0)
        if mesh_overlay:
            ax.vlines(points[:, 0], *ax.get_ylim(), color='black', linewidth=.3, alpha=.3)
    else:
        if regions is not None:
            from matplotlib.colors import ListedColormap
            ax.tripcolor(points[:, 0], points[:, 1], cells[:, :3],
                         facecolors=regions[0].astype(float), cmap=ListedColormap(colours),
                         vmin=-.5, vmax=len(colours)-.5, alpha=.25, zorder=0)
        if mesh_overlay:
            ax.triplot(points[:, 0], points[:, 1], cells[:, :3], color='black', linewidth=.3, alpha=.4)
    if regions is not None:
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=colour, label=label, alpha=.5)
                           for colour, label in zip(colours, regions[1])], fontsize='small')
    return message


def populate_3d(plotter: Any, mode: Any, mesh: Any, component: str, quantity: str, *,
                mesh_overlay: bool = False, material: bool = False,
                boundary_facets: Any = None) -> str:
    """Draw one 3D periodic mode into a PyVista plotter; return a status message.

    ``E``/``H`` and transverse ``Et``/``Ht`` are arrows; Cartesian components are coloured points.
    The mesh is a wireframe; Material shows translucent regions and boundary surfaces.
    Boundary facets contain triangle node indices, grouped by boundary kind.
    """
    import pyvista as pv

    fields = mode.fields
    points = np.asarray(fields.coordinates, dtype=float)
    if component in ('E', 'H', *_TRANSVERSE_COMPONENTS):
        axes = 'xy' if component in _TRANSVERSE_COMPONENTS else 'xyz'
        components = [fields.component(f"{component[0]}{axis}")
                      if f'{component[0]}{axis}' in fields.components else np.zeros(len(points))
                      for axis in axes]
        if component in _TRANSVERSE_COMPONENTS:
            components.append(np.zeros_like(components[0]))
            if quantity == 'magnitude':
                components = aligned_components(components)
        add_arrows(plotter, points, components,
                   quantity, label=f"|{component}|")
    else:
        cloud = pv.PolyData(points)
        label = f"{quantity}({component})"
        cloud[label] = np.asarray(fields.quantity(component, quantity))
        plotter.add_mesh(cloud, scalars=label, cmap="RdBu_r" if quantity in ("real", "imag") else "viridis",
                         render_points_as_spheres=True, point_size=6)
    return add_3d_overlays(plotter, mesh, mesh_overlay=mesh_overlay, material=material,
                           boundary_facets=boundary_facets)


def add_3d_overlays(plotter, mesh, *, mesh_overlay, material, boundary_facets=None):
    import pyvista as pv

    message = ""
    if mesh is not None and (mesh_overlay or material):
        grid = _tetra_grid(mesh)
        if material:
            legend = []
            regions = cell_materials(mesh)
            if regions is None:
                message = "This result records no material regions."
            else:
                index, labels = regions
                colours = _region_colours(labels)
                labels = [label.replace("εr", "eps_r") for label in labels]  # VTK fonts lack Greek
                for region, (label, colour) in enumerate(zip(labels, colours)):
                    if label.split(" (")[0] in ("background", "exterior"):
                        legend.append([label, "#8c8c8c"])
                        continue  # leave the background open so the field stays visible
                    legend.append([label, _hex(colour)])
                    cells = grid.extract_cells(np.flatnonzero(index == region))
                    plotter.add_mesh(cells, color=_hex(colour), opacity=0.35, show_edges=False)
                plotter.add_mesh(grid.extract_surface(), style="wireframe", color="grey", opacity=0.3)
            boundaries = (mesh.metadata.get("boundary_facets", {})
                          if boundary_facets is None else boundary_facets)
            for kind, colour in (("pec", "#666666"), ("pmc", "#9467bd"), ("sibc", "#b87333")):
                groups = [np.asarray(boundaries[name], dtype=np.int64).reshape(-1, 3)
                          for name in (kind, f"outer_{kind}")
                          if name in boundaries and np.asarray(boundaries[name]).size]
                if not groups:
                    continue
                triangles = np.concatenate(groups)
                faces = np.column_stack((np.full(len(triangles), 3), triangles))
                surface = pv.PolyData(np.asarray(mesh.coordinates, dtype=float), faces.ravel())
                plotter.add_mesh(surface, name=f"boundary_{kind}", color=colour,
                                 opacity=0.45, show_edges=False)
                legend.append([kind.upper(), colour])
            for name, triangles in boundaries.items():
                if name.startswith("port:"):
                    faces = np.column_stack((np.full(len(triangles), 3), triangles))
                    surface = pv.PolyData(mesh.coordinates, faces.ravel())
                    plotter.add_mesh(surface, name=name, color="teal", opacity=.7)
                    legend.append([name.removeprefix("port:"), "teal"])
            if legend:
                plotter.add_legend(legend, bcolor="white", face="rectangle", loc="upper left",
                                   size=(0.45, 0.05 * len(legend) + 0.02))
        if mesh_overlay:
            plotter.add_mesh(grid, style="wireframe", color="black", opacity=0.15)
    elif mesh is None and (mesh_overlay or material):
        message = "This result records no mesh."
    plotter.add_axes()
    return message


def _hex(colour: Any) -> str:
    from matplotlib.colors import to_hex

    return to_hex(colour)


def _region_colours(labels: list[str]) -> list[Any]:
    """Light grey for the background, distinct qualitative colours for other regions."""
    from matplotlib import colormaps

    palette = colormaps["tab10"].colors
    colours, others = [], 0
    for label in labels:
        if label.split(" (")[0] in ("background", "exterior"):
            colours.append((0.88, 0.88, 0.88))
        else:
            colours.append(palette[others % len(palette)])
            others += 1
    return colours


def _tetra_grid(mesh: Any):
    import pyvista as pv

    cells = np.asarray(mesh.elements)[:, :4]
    connectivity = np.hstack([np.full((len(cells), 1), 4), cells]).ravel()
    return pv.UnstructuredGrid(connectivity, np.full(len(cells), pv.CellType.TETRA),
                               np.asarray(mesh.coordinates, dtype=float))


def add_arrows(plotter: Any, points: np.ndarray, components: list[np.ndarray], quantity: str,
               *, label: str, max_arrows: int = _MAX_ARROWS, actor_name=None) -> Any:
    """Draw a complex 3D vector field as arrows coloured by its magnitude; return the arrows.

    ``real``/``imag`` draw that part of the phasor; ``magnitude``/``phase`` draw
    the real part (the field at ``omega*t = 0``) coloured by the magnitude.
    """
    import pyvista as pv

    vectors = np.stack([np.asarray(c) for c in components], axis=1)
    shown = np.imag(vectors) if quantity == "imag" else np.real(vectors)
    magnitude = np.linalg.norm(np.abs(vectors), axis=1)
    step = max(1, int(np.ceil(len(points) / max_arrows)))
    cloud = pv.PolyData(np.asarray(points, dtype=float)[::step])
    cloud["vectors"] = shown[::step]
    cloud[label] = magnitude[::step]
    span = float(np.max(np.ptp(points, axis=0))) or 1.0
    largest = float(magnitude.max() if actor_name else np.max(np.linalg.norm(shown, axis=1))) or 1.0
    arrows = cloud.glyph(orient="vectors", scale="vectors", factor=0.08 * span / largest)
    if actor_name and actor_name in plotter.actors:
        # Update the glyphs on phase changes without reframing the camera.
        arrows.set_active_scalars(label)
        mapper = plotter.actors[actor_name].mapper
        mapper.dataset.copy_from(arrows, deep=False)
        mapper.update()
    else:
        plotter.add_mesh(arrows, scalars=label, cmap="viridis", name=actor_name)
    return arrows


def show(result: Any, *, block: bool = True) -> ResultViewer:
    """Open the interactive viewer for a result object or HDF5 archive path."""
    if isinstance(result, (str, os.PathLike)):
        from .persistence import load

        result = load(result)
    if os.environ.get("QT_QPA_PLATFORM") == "offscreen":
        block = False
    return ResultViewer(result).show(block=block)


__all__ = ["ResultViewer", "add_arrows", "cell_materials", "populate_3d", "show"]
