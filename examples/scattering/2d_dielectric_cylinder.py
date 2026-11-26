"""Plot TE and TM plane-wave scattering by an infinite dielectric cylinder."""
from pathlib import Path

import fem

output = Path(__file__).parent / "outputs" / Path(__file__).stem
output.mkdir(parents=True, exist_ok=True)
for polarization in ("TE", "TM"):
    solver = fem.ScatteringSolver2D(frequency=3e9, x_range=(-.163, .163), y_range=(-.163, .163),
                                    polarization=polarization, direction=(1, 0))
    solver.add_geometry(shape=fem.shapes.Circle(center=(0, 0), radius=.02),
                        material=fem.Material(name="dielectric", epsilon=2.25))
    solver.add_pml(thickness="auto")
    solver.add_nf2ff(distance="auto")
    result = solver.solve(max_refinements=2, backend="auto")
    print(f"{polarization} backscatter width = {result.radar_cross_section([[-1, 0]])[0]:.5g} m")
    result.plot().savefig(output / f"{polarization}.png", dpi=150)
    solver.show()
