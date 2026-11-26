"""Compare gap-fed PEC dipole directivity and radiation resistance with the short-dipole limit.

Reference: MIT 6.013 Electromagnetics and Applications, lecture 19, p. 3.
The finite arms have a total tip-to-tip length of 0.12 free-space wavelengths.
"""

import argparse
from pathlib import Path
import fem

import numpy as np
from matplotlib.figure import Figure

ROOT = Path(__file__).resolve().parents[2]


def build_dipole():
    solver = fem.AntennaSolver3D(
        frequency=fem.C_0, x_range=(-0.75, 0.75), y_range=(-0.75, 0.75), z_range=(-0.75, 0.75)
    )
    radius = 0.003
    for arm in ((-0.06, -0.002), (0.002, 0.06)):
        solver.add_geometry(
            shape=fem.shapes.Cylinder(center=(0, 0), radius=radius, z_range=arm),
            material=fem.materials.PEC,
        )
    solver.add_discrete_port(
        name="feed",
        axis="z",
        gap=fem.shapes.Cylinder(center=(0, 0), radius=radius, z_range=(-0.002, 0.002)),
        impedance=50.0,
    )
    return solver


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=float, default=[0.082, 0.064])
    parser.add_argument("--backend", choices=("cpu", "gpu", "auto"), default="auto")
    parser.add_argument("--element-order", choices=(1, 2), type=int, default=2)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/benchmarks/analytical/short_dipole"
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    theta = np.linspace(0, np.pi, 181)
    reference = 1.5 * np.sin(theta) ** 2
    angles = np.rad2deg(theta)
    columns = [angles, reference]
    labels = ["theta_deg", "analytical_directivity"]
    # A short, center-fed wire with triangular current has effective length L/2.
    exact_resistance = 20 * np.pi**2 * 0.12**2
    power_rows = []
    fig = Figure(figsize=(7, 4), layout="constrained")
    ax = fig.subplots()
    ax.plot(angles, reference, color="black", label="Short-dipole limit")
    for size in sorted(set(args.sizes), reverse=True):
        solver = build_dipole()
        solver.mesh(max_element_size=size, element_order=args.element_order)
        case = solver.solve(backend=args.backend, max_refinements=0)
        errors = []
        for phi in (0.0, 90.0):
            az = np.deg2rad(phi)
            directions = np.column_stack(
                (np.sin(theta) * np.cos(az), np.sin(theta) * np.sin(az), np.cos(theta))
            )
            directivity = case.far_field(directions).directivity
            errors.append(np.max(np.abs(directivity - reference)) / reference.max())
            columns.append(directivity)
            labels.append(f"directivity_h{size:g}_phi{phi:g}")
            ax.plot(angles, directivity, "--", label=f"h={size:g} m, azimuth {phi:g}°")
        error = max(errors)
        resistance = 2 * case.radiated_power / abs(case.port_current["feed"]) ** 2
        power_rows.append(
            (size, case.accepted_power, case.radiated_power, resistance, exact_resistance)
        )
        print(
            f"h={size:g} m   broadside D={directivity[90]:.5f} (ideal 1.5)   "
            f"max error / ideal peak={error:.2%}   accepted={case.accepted_power:.6g} W   "
            f"radiated={case.radiated_power:.6g} W   Rrad={resistance:.3f} ohm "
            f"(short-wire approximation {exact_resistance:.3f} ohm)"
        )
    ax.set(xlabel="Polar angle (degrees)", ylabel="Directivity (linear)")
    ax.legend()
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.output / "comparison.csv",
        np.column_stack(columns),
        delimiter=",",
        header=",".join(labels),
        comments="",
    )
    fig.savefig(args.output / "comparison.png", dpi=150)
    np.savetxt(
        args.output / "power.csv",
        power_rows,
        delimiter=",",
        comments="",
        header="h_m,accepted_W,radiated_W,radiation_resistance_ohm,short_wire_resistance_ohm",
    )
    if args.check and not (
        error < 0.10 and abs(resistance / exact_resistance - 1) < 0.25 and case.accepted_power > 0
    ):
        raise SystemExit(
            "Short-dipole check failed for directivity, radiation resistance or accepted power."
        )


if __name__ == "__main__":
    main()
