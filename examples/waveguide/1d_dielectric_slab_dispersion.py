"""Guided-mode dispersion of a silicon slab in silica."""
from pathlib import Path

import numpy as np
import fem

OUTPUT = Path(__file__).parent / "outputs" / Path(__file__).stem
OUTPUT.mkdir(parents=True, exist_ok=True)

cladding = fem.Material(name="silica cladding", epsilon=1.44 ** 2)
silicon = fem.Material(name="silicon core", epsilon=3.45 ** 2)
frequencies = np.linspace(180e12, 210e12, 11)
neff_sweep = []
solver = fem.ModeSolver1D(frequency=frequencies[-1], x_range=(-3e-6, 3e-6),
                      background_material=cladding)
solver.add_layer(x_range=(-250e-9, 250e-9), name="core", material=silicon)
cases = solver.sweep(frequencies=frequencies, max_refinements=2, backend="auto", num_modes=4)
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
