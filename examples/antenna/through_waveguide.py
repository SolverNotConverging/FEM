"""Two-port network of a straight dielectric-filled rectangular guide."""

from pathlib import Path
import numpy as np
import fem

solver = fem.AntennaSolver3D(
    frequency=8e9, x_range=(-0.036, 0.036), y_range=(-0.030, 0.030), z_range=(-0.040, 0.040),
)
u, v, ends = (-0.010, 0.010), (-0.004, 0.004), (-0.015, 0.015)
bore = fem.shapes.Box(bounds=(u, v, ends))
housing = fem.shapes.Box(bounds=((-0.012, 0.012), (-0.006, 0.006), (-0.017, 0.017)))
solver.add_geometry(
    shape=fem.shapes.Difference(shape=housing, tool=bore), material=fem.materials.PEC
)
solver.add_geometry(shape=bore, material=fem.materials.Material(name="Dielectric", epsilon=2.25))
for name, z, drive in (("input", ends[0], 1.0), ("output", ends[1], 0.0)):
    solver.add_waveguide_port(
        name=name,
        axis="z",
        position=z,
        cross_section=fem.shapes.Rectangle(bounds=(u, v)),
        num_modes=3,
        mode=1,
        amplitude=drive,
    )

solver.add_pml(thickness="auto")
solver.add_nf2ff(distance="auto")
sweep = solver.sweep(frequencies=np.linspace(7e9, 9e9, 5), backend="auto")
for case in sweep:
    print(
        f"{case.frequency / 1e9:g} GHz   S11={case.s_parameters[0, 0]:.4f}   S21={case.s_parameters[1, 0]:.4f}"
    )
output = Path(__file__).parent / "outputs" / "through_waveguide"
sweep.save(output / "network.h5")
sweep.plot().savefig(output / "network.png")
sweep.show()
