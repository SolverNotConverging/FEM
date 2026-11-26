"""H-plane sectoral horn with a straight throat and flared conducting walls."""

import numpy as np
import fem


def flared_section(throat, mouth, bottom, top, half_height):
    outline = fem.shapes.Polygon(
        points=(
            (-throat, bottom),
            (throat, bottom),
            (throat, 0),
            (mouth, top),
            (-mouth, top),
            (-throat, 0),
        )
    )
    solid = fem.shapes.Extrusion(shape=outline, z_range=(-half_height, half_height))
    return solid.rotated(angle=90, axis=(1, 0, 0))


solver = fem.AntennaSolver3D(
    frequency=10e9, x_range=(-0.054, 0.054), y_range=(-0.036, 0.036), z_range=(-0.047, 0.068),
)
outer = flared_section(0.0125, 0.027, -0.018, 0.032, 0.006)
inner = flared_section(0.01143, 0.026, -0.016, 0.033, 0.00508)
solver.add_geometry(
    shape=fem.shapes.Difference(shape=outer, tool=inner), material=fem.materials.PEC
)
throat = fem.shapes.Rectangle(bounds=((-0.01143, 0.01143), (-0.00508, 0.00508)))
solver.add_waveguide_port(
    name="throat", axis="z", position=-0.016, cross_section=throat, num_modes=3, mode=1
)
solver.add_pml(thickness="auto")
solver.add_nf2ff(distance="auto")
sweep = solver.sweep(frequencies=np.linspace(9e9, 11e9, 5), max_refinements=2, backend="auto")
for case in sweep:
    gain = case.far_field([[0, 0, 1]]).realized_gain[0]
    print(f"{case.frequency / 1e9:.2f} GHz   S11={case.s_parameters[0, 0]:.4f}"
          f"   Boresight realized gain={10 * np.log10(gain):.2f} dBi")
sweep.show()
