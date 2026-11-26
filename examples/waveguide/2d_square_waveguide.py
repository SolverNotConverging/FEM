"""Solve and display 24 full-vector modes of a square PEC waveguide."""
import fem

solver = fem.ModeSolver2D(
    frequency=60e9,
    x_range=(0., 20e-3),
    y_range=(0., 20e-3),
)
modes = solver.solve(
    num_modes=24,
    track_modes=True,
    max_refinements=2,
    backend="cpu",
)

for mode in modes:
    print(f"Mode {mode.index}: neff={mode.neff:.9g}, {mode.polarization}, "
          f"irrep={mode.metadata['symmetry_irrep']}, "
          f"partner={mode.metadata['symmetry_partner']}")

# Select each of the 24 modes and its field components in the viewer.
modes.show()
