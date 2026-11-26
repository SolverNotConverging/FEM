"""Complex dispersion of a copper microstrip with surface impedance."""
from pathlib import Path

import numpy as np
import fem

OUTPUT = Path(__file__).parent / "outputs" / Path(__file__).stem
OUTPUT.mkdir(parents=True, exist_ok=True)

substrate = fem.Material(name="microwave laminate", epsilon=3.55 * (1. - 1j * 0.0027))
frequencies = np.linspace(8e9, 12e9, 9)
neff_sweep = []
solver = fem.ModeSolver2D(frequency=frequencies[-1], x_range=(-6e-3, 6e-3),
                      y_range=(-35e-6, 6e-3), boundary=fem.materials.PEC, background_material=fem.materials.air)
solver.add_rectangle(x_range=(-6e-3, 6e-3), y_range=(0., 1.524e-3),
                     name="substrate", material=substrate)
solver.add_geometry(name="copper_ground", material=fem.materials.copper,
                    shape=fem.shapes.Rectangle(bounds=((-6e-3, 6e-3), (-35e-6, 0.))))
solver.add_geometry(name="copper_strip", material=fem.materials.copper,
                    shape=fem.shapes.Rectangle(bounds=((-1.5e-3, 1.5e-3), (1.524e-3, 1.559e-3))))
cases = solver.sweep(frequencies=frequencies, max_refinements=2, backend="auto", num_modes=1)
for case, modes in enumerate(cases, start=1):
    modes.save(OUTPUT / f"case_{case:03d}.h5")
    neff_sweep.append(modes.neff)

rows = [(frequency, mode, value.real, value.imag)
        for frequency, values in zip(frequencies, neff_sweep)
        for mode, value in enumerate(values)]
np.savetxt(OUTPUT / "dispersion.csv", rows, delimiter=",", comments="",
           header="frequency,mode,neff_real,neff_imag",
           fmt=("%.12e", "%d", "%.12e", "%.12e"))
fem.plot_dispersion(frequencies, neff_sweep)
