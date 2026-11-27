"""One native FEM refinement of the existing wave-port patch with each marking policy."""

import argparse
import csv
from pathlib import Path

import fem
import numpy as np


def inset_patch():
    """Build the fixed one-layer PML reference model used for the marking comparison."""
    solver = fem.AntennaSolver3D(
        frequency=2.45e9, x_range=(-0.091, 0.091), y_range=(-0.091, 0.091),
        z_range=(-0.069, 0.077), pml_layers=1,
    )
    board = (-0.029, 0.029)
    height = 0.0016
    copper = fem.materials.copper
    loss_tangent = 0.02
    substrate = fem.Material(name="FR-4", epsilon=4.4 * (1 - 1j * loss_tangent))
    solver.add_geometry(
        shape=fem.shapes.Box(bounds=(board, board, (-0.0002, 0.0))),
        material=copper,
    )
    solver.add_geometry(
        shape=fem.shapes.Box(bounds=(board, board, (0.0, height))),
        material=substrate,
    )
    # A 3 mm strip enters a 5 mm inset and joins the patch on the same metal layer.
    metal_z = (height, height + 0.0002)
    port_x = board[0]
    patch_length = 0.0276
    patch_edge = -patch_length / 2
    inset_end = patch_edge + 0.009
    patch = fem.shapes.Box(bounds=((patch_edge, -patch_edge), (-0.0185, 0.0185), metal_z))
    inset = fem.shapes.Box(bounds=((patch_edge - 0.001, inset_end), (-0.0025, 0.0025), metal_z))
    line = fem.shapes.Box(bounds=((port_x, inset_end + 0.001), (-0.0015, 0.0015), metal_z))
    solver.add_geometry(
        shape=fem.shapes.Union(shapes=(line, fem.shapes.Difference(shape=patch, tool=inset))),
        material=copper,
    )
    # A short launch enclosure defines the port perimeter at the board edge.
    # Its forward end opens onto the uniform microstrip before the inset.
    port_y, port_z = (-0.009, 0.009), (0., 0.008)
    launch = fem.shapes.Box(bounds=((port_x - 0.0005, port_x + 0.002),
                                    (-0.0095, 0.0095), (0., 0.0085)))
    bore = fem.shapes.Box(bounds=((port_x, port_x + 0.003), port_y, port_z))
    solver.add_geometry(shape=fem.shapes.Difference(shape=launch, tool=bore), material=copper)
    solver.add_waveguide_port(
        name="microstrip", axis="x", position=port_x,
        cross_section=fem.shapes.Rectangle(bounds=(port_y, port_z)), mode=1, num_modes=3,
    )
    solver.add_pml(thickness="auto")
    solver.add_nf2ff(distance="auto")
    return solver


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = []
    for label, weights in (("unit", {}), ("fitted", {"interface": 4.0, "element_pde": 4.0})):
        solver = inset_patch()
        solver.mesh(wavelength_elements=4)
        solver.configure_adaptivity(marking_weights=weights)
        solver.configure_s_parameter_convergence(tolerance=0.02)
        reflections = []

        def record(result):
            reflections.append(result.s_parameters[0, 0])
            return result.s_parameters

        result = solver.solve(backend="cpu", max_refinements=1, pml_observable=record)
        for step, (row, s11) in enumerate(
            zip(result.solve_info["adaptive_history"], reflections, strict=True), 1
        ):
            rows.append(
                [
                    label,
                    step,
                    row["elements"],
                    row["residual"],
                    s11.real,
                    s11.imag,
                    20 * np.log10(abs(s11)),
                ]
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "marking",
                "pass",
                "physical_cells",
                "physical_residual",
                "s11_real",
                "s11_imag",
                "s11_db",
            ]
        )
        w.writerows(rows)


if __name__ == "__main__":
    main()
