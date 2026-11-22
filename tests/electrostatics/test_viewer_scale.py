"""Mesh selection preserves the physical bounds and GUI drawing area."""
from types import SimpleNamespace
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fem.electrostatics.result_api import show_result


def test_mesh_and_fields_use_same_spatial_scale(monkeypatch):
    result = SimpleNamespace(
        coordinates=np.array([[0., 0.], [4e-3, 0.], [4e-3, 1e-3], [0., 1e-3]]),
        elements=np.array([[0, 1, 2], [0, 2, 3]]),
        potential=np.array([0., 1., 1., 0.]),
        element_electric_field=np.array([[-2., 3.], [4., -5.]]),
        element_displacement_field=np.array([[3e-11, -4e-11], [-1e-11, 2e-11]]),
    )
    monkeypatch.setattr(plt, "show", lambda **kwargs: None)
    figure = show_result(result, block=False)
    try:
        figure.canvas.draw()
        reference = figure.axes[1].get_position().bounds
        from matplotlib.quiver import Quiver
        labels = [label.get_text() for label in figure._fem_selector.labels]
        for name, field in (("Et", result.element_electric_field), ("Dt", result.element_displacement_field)):
            figure._fem_selector.set_active(labels.index(name))
            figure.canvas.draw()
            arrows = next(a for a in figure.axes[1].collections if isinstance(a, Quiver))
            np.testing.assert_allclose(arrows.get_offsets(), result.coordinates[result.elements].mean(axis=1))
            np.testing.assert_allclose([arrows.U, arrows.V], field.T, atol=0)
            np.testing.assert_allclose(arrows.get_array(), np.linalg.norm(field, axis=1), atol=0)
        for index in (5, 1, 0):  # Mesh, Ex, potential.
            figure._fem_selector.set_active(index)
            figure.canvas.draw()
            axis = figure.axes[1]
            np.testing.assert_allclose(axis.get_position().bounds, reference)
            np.testing.assert_allclose(axis.get_xlim(), (0., 4e-3))
            np.testing.assert_allclose(axis.get_ylim(), (0., 1e-3))
            assert axis.get_aspect() == 1
            assert axis.get_xlabel() == "x (m)"
            assert axis.get_ylabel() == "y (m)"
    finally:
        plt.close(figure)
