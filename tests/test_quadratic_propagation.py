"""On-mesh accuracy of quadratic scalar and periodic Maxwell spaces."""
import numpy as np
import pytest

import fem


@pytest.mark.parametrize("solver_type,geometry,size,tolerances", [
    (fem.ModeSolver1D, {}, .20, (.005, .00004)),
    (fem.PeriodicModeSolver2D, {"z_range": .25, "polarization": "TE"}, .20, (.004, .000005)),
    (fem.PeriodicModeSolver3D, {"z_range": .25, "y_range": .5}, .12, (.003, .000004)),
])
def test_propagation_against_analytic_cutoff(solver_type, geometry, size, tolerances):
    exact = np.sqrt(1 - .5**2)  # One-metre wavelength and one-metre broad wall.
    for order in (1, 2):
        solver = solver_type(frequency=fem.C_0, x_range=1., **geometry)
        mesh = solver.mesh(max_element_size=size, **({"element_order": 1} if order == 1 else {}))
        mode = solver.solve(backend="cpu", max_refinements=0, num_modes=1, neff_guess=.865)[0]
        assert mesh.info.element_order == order
        assert abs(mode.neff - exact) < tolerances[order - 1]
