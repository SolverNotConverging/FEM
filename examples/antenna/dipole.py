"""Drive a round-wire PEC dipole across a cylindrical 50-ohm discrete port."""

import numpy as np
import fem

solver = fem.AntennaSolver3D(
    frequency=fem.C_0, x_range=(-0.75, 0.75), y_range=(-0.75, 0.75), z_range=(-0.75, 0.75),
)
radius = 0.003
for arm in ((-0.235, -0.002), (0.002, 0.235)):
    solver.add_geometry(shape=fem.shapes.Cylinder(center=(0, 0), radius=radius, z_range=arm),
                        material=fem.materials.PEC)
solver.add_discrete_port(
    name="feed",
    axis="z",
    gap=fem.shapes.Cylinder(center=(0, 0), radius=radius, z_range=(-0.002, 0.002)),
    impedance=50.0,
)

solver.add_pml(thickness="auto")
solver.add_nf2ff(distance="auto")
sweep = solver.sweep(frequencies=fem.C_0 * np.linspace(0.9, 1.1, 5), backend="auto")
for case in sweep:
    print(f"{case.frequency / 1e6:.1f} MHz   S11={case.s_parameters[0, 0]:.4f}"
          f"   Feed voltage={case.port_voltage['feed']:.4g} V"
          f"   Feed current={case.port_current['feed']:.4g} A")
sweep.show()
