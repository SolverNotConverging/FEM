"""Retest two FEM convergence gates while varying window width and absorber resolution."""

import argparse
import csv
from pathlib import Path

import fem
import numpy as np
from matplotlib.figure import Figure
from scipy.optimize import brentq
from scipy.special import hankel2, jv


def slab(width, resolution):
    k = 2*np.pi*10e9/fem.C_0
    half = .005
    v = k*half*np.sqrt(3.)
    u = brentq(lambda x: x*np.tan(x)-np.sqrt(v*v-x*x), 1e-9, min(v, np.pi/2)-1e-9)
    exact = np.sqrt(4-(u/(k*half))**2)
    solver = fem.ModeSolver1D(frequency=10e9, x_range=width)
    solver.add_layer(x_range=(width/2-half, width/2+half), material=fem.Material(epsilon=4.))
    solver.mesh(resolution=resolution)
    modes = solver.solve(backend="cpu", num_modes=4, neff_guess=exact, max_refinements=0)
    mode = min((m for m in modes if m.polarization == "TE"), key=lambda m: abs(m.neff-exact))
    distance = abs(mode.fields.coordinates[0]-width/2)
    profile = np.where(distance <= half, np.cos(u*distance/half),
                       np.cos(u)*np.exp(-k*np.sqrt(exact**2-1)*(distance-half)))
    field = mode.component("Ey")
    amplitude = np.vdot(profile, field)/np.vdot(profile, profile)
    return np.linalg.norm(field/amplitude-profile)/np.linalg.norm(profile)


def cylinder(layers, resolution):
    frequency, radius = 3e9, .025
    k = 2*np.pi*frequency/fem.C_0
    n = np.arange(-15, 16)
    angles = np.linspace(0, np.pi, 31)
    coefficients = -jv(n, k*radius)/hankel2(n, k*radius)
    exact = 4/k*abs(np.exp(1j*angles[:, None]*n) @ coefficients)**2
    solver = fem.ScatteringSolver2D(frequency=frequency, x_range=(-.18, .18), y_range=(-.18, .18),
                                    polarization="TM", pml_layers=layers)
    solver.add_geometry(shape=fem.shapes.Circle(center=(0, 0), radius=radius), material=fem.materials.PEC)
    solver.mesh(max_element_size=fem.C_0/frequency/resolution)
    result = solver.solve(backend="cpu", max_refinements=0)
    width = result.radar_cross_section(np.column_stack((np.cos(angles), np.sin(angles))))
    return np.linalg.norm(width-exact)/np.linalg.norm(exact)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    fig = Figure(figsize=(8.5, 3.8), layout="constrained")
    for ax, case, configs in zip(fig.subplots(1, 2), ("Slab field", "PEC cylinder width"),
                                 ([(.06, [120, 240]), (.12, [240, 480])], [(5, [24, 40]), (10, [24, 40])])):
        for value, resolutions in configs:
            is_slab = case == "Slab field"
            errors = [(slab if is_slab else cylinder)(value, n) for n in resolutions]
            label = f"{value*1000:g} mm window" if is_slab else f"{value} PML layers"
            for n, error in zip(resolutions, errors):
                rows.append(dict(case=case, control=label, resolution=n, relative_error=error))
            h = [1000*(value/n if is_slab else fem.C_0/(3e9*n)) for n in resolutions]
            ax.semilogy(h, 100*np.array(errors), marker="o", label=label)
        ax.set(title=case, xlabel="Nominal spacing (mm)", ylabel="Relative error (%)")
        ax.invert_xaxis()
        ax.grid(alpha=.2)
        ax.legend()
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    fig.savefig(args.output / "comparison.png", dpi=150)


if __name__ == "__main__":
    main()
