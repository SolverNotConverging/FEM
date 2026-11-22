"""In-plane field arrows on the native sampling mesh."""

import numpy as np


def aligned_components(components):
    values = np.asarray(components, dtype=complex)
    pivot = values.ravel()[np.argmax(abs(values))]
    return values * np.exp(-1j * np.angle(pivot))


def draw_vectors(ax, points, values, domain, *, vertical_axis, quantity, label,
                 horizontal_axis="x", colorbar=True, colorbar_axes=None, cmap="viridis"):
    """Plot signed arrows at cell samples, coloured by the full in-plane norm."""
    step = max(1, int(np.ceil(len(points) / 600)))
    points, values = points[::step], values[:, ::step]
    magnitude = np.linalg.norm(values, axis=0)
    if quantity in ("magnitude", "abs") and np.iscomplexobj(values):
        values = aligned_components(values)
    elif quantity not in ("real", "imag", "magnitude", "abs"):
        raise ValueError("Vector arrows support real, imag, or magnitude.")
    vectors = values.imag if quantity == "imag" else values.real
    span = float(np.ptp(domain, axis=0).max()) or 1.0
    largest = float(np.linalg.norm(vectors, axis=0).max()) or 1.0
    y = points[:, 1] if points.shape[1] == 2 else np.zeros(len(points))
    arrows = ax.quiver(points[:, 0], y, vectors[0], vectors[1], magnitude,
                       angles="xy", scale_units="xy", scale=largest / (.06 * span),
                       pivot="mid", cmap=cmap)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"{horizontal_axis} (m)")
    ax.set_ylabel(f"{vertical_axis} (m)")
    ax.set_xlim(domain[:, 0].min() - .08 * span, domain[:, 0].max() + .08 * span)
    if domain.shape[1] == 1:
        ax.set_ylim(-.15 * span, .15 * span)
    else:
        ax.set_ylim(domain[:, 1].min() - .08 * span, domain[:, 1].max() + .08 * span)
    if colorbar:
        ax.figure.colorbar(arrows, ax=ax, cax=colorbar_axes, label=label)


def draw_mode_vectors(ax, mode, component, quantity, *, axes, colorbar=True, cmap="viridis"):
    fields = mode.fields
    points = fields.coordinates
    if isinstance(points, tuple):
        points = np.column_stack([np.asarray(axis).ravel() for axis in points])
    points = np.asarray(points).reshape(-1, fields.dimension)
    domain = fields.mesh_points if fields.mesh_points is not None else points
    values = np.zeros((2, len(points)), complex)
    for i, axis in enumerate(axes):
        name = component[0] + axis
        if name in fields.components:
            values[i] = fields.component(name).ravel()
    owners = getattr(fields, "sample_element_indices", None)
    if owners is None:
        owners = fields.metadata.get("sample_element_indices")
    if owners is not None:
        # Average each cell's own quadrature samples without crossing interfaces.
        _, owners = np.unique(np.asarray(owners).ravel(), return_inverse=True)
        counts = np.bincount(owners)
        points = np.column_stack([np.bincount(owners, weights=v) / counts for v in points.T])
        values = np.array([(np.bincount(owners, weights=v.real)
                            + 1j * np.bincount(owners, weights=v.imag)) / counts for v in values])
    unit = "V/m" if component == "Et" else "A/m"
    draw_vectors(ax, points, values, domain, vertical_axis=axes[1], quantity=quantity,
                 label=f"|{component}| ({unit})", colorbar=colorbar, cmap=cmap)
    ax.set_title(f"Mode {mode.index}: {quantity}({component}), neff={mode.neff:.6g}")
