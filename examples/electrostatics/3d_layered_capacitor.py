"""Inspect the 3D electric field and compare plate charge with a series capacitor."""

import fem

width, height, gap = .006, .008, .004
volts = 3.
solver = fem.ElectrostaticSolver(dim=3, x_range=width, y_range=height,
                                 z_range=gap, outer_potential=None)
solver.add_box(x_range=solver.x_range, y_range=solver.y_range, z_range=(gap/2, gap),
               material=fem.Material(epsilon=4.), name="dielectric")
solver.set_potential(geometry="z-", potential=0., name="ground")
solver.set_potential(geometry="z+", potential=volts, name="drive")
solver.mesh(max_element_size=.0015, boundary_refinement=None, interface_refinement=None)
result = solver.solve()
expected = fem.EPSILON_0 * width * height / (gap/2 + gap/8)
print(f"FEM capacitance {result.conductor_charge('drive') / volts * 1e12:.5f} pF")
print(f"Analytical capacitance {expected * 1e12:.5f} pF")
result.show()
