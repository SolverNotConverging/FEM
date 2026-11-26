"""Compare shell-only refinement to the analytical scattering width of a cylinder."""

import argparse
import csv
from pathlib import Path

import fem
import numpy as np
from matplotlib.figure import Figure

from amr_observables import cylinder_width


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    solver = fem.ScatteringSolver2D(frequency=fem.C_0, x_range=(-1, 1), y_range=(-1, 1),
                                    polarization="TM", pml_layers=1)
    solver.add_geometry(shape=fem.shapes.Circle(center=(0, 0), radius=.18), material=fem.Material(epsilon=2.25))
    solver.mesh(max_element_size=.08)
    angles = np.linspace(0., np.pi, 19)
    directions = np.column_stack((np.cos(angles), np.sin(angles)))
    exact = cylinder_width(.18, 1.5, angles)
    result = solver.solve(backend="cpu", adaptive_tolerance=.015, pml_tolerance=.005, max_refinements=7,
                          pml_observable=lambda fields: fields.radar_cross_section(directions))
    rows = []
    for step in result.solve_info["adaptive_history"]:
        error = np.linalg.norm(step["observable_values"] - exact) / np.linalg.norm(exact)
        rows.append(dict(step=step["refinement"], physical_cells=step["elements"], field_dofs=step["field_dofs"],
                         layers=step["pml_layers"], physical_residual=step["residual"],
                         pml_residual=step["pml_residual"], interface_residual=step["pml_interface_residual"],
                         relative_width_error=error, status=step["status"]))
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    fig = Figure(figsize=(7, 3.7), layout="constrained")
    ax = fig.subplots()
    ax.semilogy([r["step"] for r in rows], [100*r["relative_width_error"] for r in rows], "o-")
    ax.set(xlabel="Adaptive step", ylabel="Angular scattering-width error (%)",
           xticks=[r["step"] for r in rows], title="Dielectric cylinder · shell convergence")
    ax.grid(alpha=.25)
    fig.savefig(args.output / "convergence.png", dpi=150)
    if rows[-1]["relative_width_error"] >= .01 or not result.solve_info["adaptive_converged"]:
        raise SystemExit("Final angular error or shell convergence check failed.")


if __name__ == "__main__":
    main()
