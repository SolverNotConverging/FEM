"""Inset-fed copper patch on lossy FR-4 with a numerical microstrip wave port."""

import numpy as np
import fem

solver = fem.AntennaSolver3D(
    frequency=2.45e9, x_range=(-0.091, 0.091), y_range=(-0.091, 0.091),
    z_range=(-0.069, 0.077),
)
board = (-0.029, 0.029)
height = 0.0016
copper = fem.materials.copper
loss_tangent = 0.02
substrate = fem.Material(name="FR-4", epsilon=4.4 * (1 - 1j * loss_tangent))
solver.add_geometry(
    shape=fem.shapes.Box(bounds=(board, board, (-0.0002, 0.0))),
    material=copper,
)
solver.add_geometry(
    shape=fem.shapes.Box(bounds=(board, board, (0.0, height))),
    material=substrate,
)
# A 3 mm strip enters a 5 mm inset and joins the patch on the same metal layer.
metal_z = (height, height + 0.0002)
port_x = board[0]
patch_length = 0.0276
patch_edge = -patch_length / 2
inset_end = patch_edge + 0.009
patch = fem.shapes.Box(bounds=((patch_edge, -patch_edge), (-0.0185, 0.0185), metal_z))
inset = fem.shapes.Box(bounds=((patch_edge - 0.001, inset_end), (-0.0025, 0.0025), metal_z))
line = fem.shapes.Box(bounds=((port_x, inset_end + 0.001), (-0.0015, 0.0015), metal_z))
solver.add_geometry(
    shape=fem.shapes.Union(shapes=(line, fem.shapes.Difference(shape=patch, tool=inset))),
    material=copper,
)
# A short launch enclosure defines the port perimeter at the board edge.
# Its forward end opens onto the uniform microstrip before the inset.
port_y, port_z = (-0.009, 0.009), (0., 0.008)
launch = fem.shapes.Box(bounds=((port_x - 0.0005, port_x + 0.002),
                                (-0.0095, 0.0095), (0., 0.0085)))
bore = fem.shapes.Box(bounds=((port_x, port_x + 0.003), port_y, port_z))
solver.add_geometry(shape=fem.shapes.Difference(shape=launch, tool=bore), material=copper)
solver.add_waveguide_port(
    name="microstrip", axis="x", position=port_x,
    cross_section=fem.shapes.Rectangle(bounds=(port_y, port_z)), mode=1, num_modes=3,
)
solver.add_pml(thickness="auto")
solver.add_nf2ff(distance="auto")
sweep = solver.sweep(frequencies=np.linspace(2.35e9, 2.55e9, 5), max_refinements=4, backend="auto")
for case in sweep:
    gain = case.far_field([[0, 0, 1]]).realized_gain[0]
    neff = case.port_modes["microstrip"][0].neff
    print(f"{case.frequency / 1e9:.2f} GHz   Feed neff={neff:.4g}"
          f"   S11={20 * np.log10(abs(case.s_parameters[0, 0])):.2f} dB"
          f"   Broadside realized gain={10 * np.log10(gain):.2f} dBi")
sweep.show()
