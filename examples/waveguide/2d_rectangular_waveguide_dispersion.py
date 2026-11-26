"""TE10 cutoff dispersion in a rectangular waveguide."""
from pathlib import Path

import numpy as np
import fem

OUTPUT = Path(__file__).parent / "outputs" / Path(__file__).stem
OUTPUT.mkdir(parents=True, exist_ok=True)

frequencies = np.linspace(200e6, 400e6, 21)
neff_sweep = []
solver = fem.ModeSolver2D(frequency=frequencies[-1], x_range=1e0, y_range=500e-3)
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
