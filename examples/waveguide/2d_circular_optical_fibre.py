"""Solve the two polarizations of the fundamental mode of a step-index fibre."""

from math import pi, sqrt
from pathlib import Path

import fem

wavelength = 1e-6  # Vacuum wavelength in metres
core_radius = 1.3e-6
n_core = 1.50
n_cladding = 1.44
core = fem.Material(name="core glass", epsilon=n_core ** 2)
cladding = fem.Material(name="cladding glass", epsilon=n_cladding ** 2)

# Truncate the cladding several core radii away, where the guided field has decayed.
half_width = 4 * core_radius
solver = fem.ModeSolver2D(
    frequency=fem.C_0 / wavelength,
    x_range=(-half_width, half_width),
    y_range=(-half_width, half_width),
    background_material=cladding,
    boundary=fem.materials.PEC,
)
solver.add_geometry(
    name="core",
    shape=fem.shapes.Circle(center=(0.0, 0.0), radius=core_radius),
    material=core,
)
modes = solver.solve(num_modes=4, max_refinements=2, backend="auto")

v_number = 2 * pi * core_radius / wavelength * sqrt(n_core ** 2 - n_cladding ** 2)
print(f"V = {v_number:.3f}, guided index range: {n_cladding} < neff < {n_core}")
for mode in modes:
    print(f"Mode {mode.index}: neff={mode.neff:.8g}, residual={mode.residual:.2e}")

output_dir = Path(__file__).parent / "outputs" / Path(__file__).stem
output_dir.mkdir(parents=True, exist_ok=True)
modes.save(output_dir / "results.h5")
modes.show()
