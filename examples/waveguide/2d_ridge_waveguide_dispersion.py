"""Full-vector guided-mode dispersion of a silicon ridge on silica."""
from pathlib import Path

import numpy as np
import fem

OUTPUT = Path(__file__).parent / "outputs" / Path(__file__).stem
OUTPUT.mkdir(parents=True, exist_ok=True)

silicon = fem.Material(name="silicon core", epsilon=3.45**2)
silica = fem.Material(name="silica slab", epsilon=1.44**2)
frequencies = np.linspace(900e12, 1100e12, 7)
neff_sweep = []
solver = fem.ModeSolver2D(frequency=frequencies[-1], x_range=(-2e-6, 2e-6),
    y_range=(-1.5e-6, 1.5e-6), boundary=fem.materials.PEC, background_material=fem.materials.vacuum)
solver.add_rectangle(x_range=(-100e-9, 100e-9), y_range=(-110e-9, 110e-9),
    name="silicon_core", material=silicon)
solver.add_rectangle(x_range=(-1e-6, 1e-6), y_range=(-200e-9, -110e-9),
    name="slab", material=silica)
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
