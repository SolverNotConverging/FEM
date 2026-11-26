"""View plane-wave scattering by a copper sphere and export both VTK blocks."""
from pathlib import Path

import fem

output = Path(__file__).parent / "outputs" / Path(__file__).stem
solver = fem.ScatteringSolver3D(frequency=3e9, x_range=(-.093, .093), y_range=(-.093, .093),
                                z_range=(-.093, .093), direction=(1, 0, 0),
                                electric_vector=(0, 1, 0))
solver.add_geometry(shape=fem.shapes.Sphere(center=(0, 0, 0), radius=.015),
                    material=fem.materials.copper)
solver.add_pml(thickness="auto")
solver.add_nf2ff(distance="auto")
result = solver.solve(max_refinements=2, backend="auto")
print(f"Backscatter RCS = {result.radar_cross_section([[-1, 0, 0]])[0]:.5g} m²")
result.export_vtk(output / "sphere.vtm")
solver.show()
