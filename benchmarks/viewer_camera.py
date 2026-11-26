"""Measure camera drift and update time through the viewer's field and overlay controls."""

import argparse
import csv
import os
from pathlib import Path
from time import perf_counter

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
from matplotlib.figure import Figure

from analytical.coaxial_waveguide_adaptivity import make_solver


def projected_points(viewer, anchors):
    if viewer._canvas is not None:
        return viewer._canvas.figure.axes[0].transData.transform(anchors[:, :2])
    # Qt offscreen can have a zero-sized native VTK framebuffer. Project via
    # the actual camera into the laid-out widget's logical pixel coordinates.
    rect = viewer._body.geometry()
    width, height = rect.width(), rect.height()
    if min(width, height) <= 1:
        raise RuntimeError("Viewer has no measurable viewport.")
    transform = viewer._plotter.camera.GetCompositeProjectionTransformMatrix(width/height, -1., 1.)
    matrix = np.array([[transform.GetElement(i, j) for j in range(4)] for i in range(4)])
    clip = np.column_stack((anchors, np.ones(len(anchors)))) @ matrix.T
    return (1 + clip[:, :2] / clip[:, 3:]) * np.array([width, height])/2


def measure(result, renderer, fields):
    view = result.show(block=False)
    anchors = np.array([[.001, 0, 0], [-.001, .002, 0], [.002, .001, 0]])
    rows = []
    try:
        view._component.setCurrentText(fields[0])
        view._app.processEvents()
        if renderer == "Matplotlib":
            axes = view._canvas.figure.axes[0]
            axes.set(xlim=(-.0025, .002), ylim=(-.0015, .003))
            view._canvas.draw()
        else:
            camera = view._plotter.camera
            camera.Azimuth(32)
            camera.Elevation(17)
            camera.Zoom(1.6)
        baseline = projected_points(view, anchors)
        for repeat in range(4):
            for action in ("Field", "Mesh", "Material"):
                start = perf_counter()
                if action == "Field":
                    view._component.setCurrentText(fields[(repeat+1) % 2])
                elif action == "Mesh":
                    view._mesh.click()
                else:
                    view._material.click()
                view._app.processEvents()
                if view._status.currentMessage():
                    raise RuntimeError(view._status.currentMessage())
                elapsed = perf_counter()-start
                drift = np.linalg.norm(projected_points(view, anchors)-baseline, axis=1).max()
                rows.append(dict(renderer=renderer, action=action, repeat=repeat, elapsed_ms=1000*elapsed,
                                 max_anchor_shift_px=drift))
        start = perf_counter()
        view._reset_view.click()
        view._app.processEvents()
        elapsed = perf_counter()-start
        drift = np.linalg.norm(projected_points(view, anchors)-baseline, axis=1).max()
        rows.append(dict(renderer=renderer, action="Reset view", repeat=0, elapsed_ms=1000*elapsed,
                         max_anchor_shift_px=drift))
        return rows
    finally:
        view._clear_body()
        view._window.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = make_solver().solve(backend="cpu", num_modes=1, max_refinements=0)
    rows = measure(result, "Matplotlib", ["Ex", "Ey"]) + measure(result, "VTK", ["E vector", "H vector"])
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "samples.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    fig = Figure(figsize=(8.2, 3.7), layout="constrained")
    timing, shift = fig.subplots(1, 2)
    actions = ["Field", "Mesh", "Material", "Reset view"]
    for j, renderer in enumerate(("Matplotlib", "VTK")):
        for ax, key, reduce in [(timing, "elapsed_ms", np.median), (shift, "max_anchor_shift_px", max)]:
            values = [reduce([r[key] for r in rows if r["renderer"] == renderer and r["action"] == a]) for a in actions]
            ax.bar(np.arange(4)+(j-.5)*.36, values, width=.36, label=renderer)
            ax.set_xticks(range(4), actions, rotation=20)
    timing.set_ylabel("Median update time (ms)")
    shift.set_ylabel("Maximum screen shift (pixels)")
    timing.legend()
    fig.savefig(args.output / "comparison.png", dpi=150)
    if any((r["max_anchor_shift_px"] <= 1 if r["action"] == "Reset view" else r["max_anchor_shift_px"] > 1)
           for r in rows):
        raise SystemExit("Expected stable framing for field/overlay controls and a visible reset.")


if __name__ == "__main__":
    main()
