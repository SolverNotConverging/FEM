"""Opt-in real-device CPU/CUDA agreement checks."""
import os
import numpy as np
import pytest
from scipy.constants import c

pytestmark = pytest.mark.skipif(os.environ.get("FEM_TEST_CUDA") != "1", reason="Set FEM_TEST_CUDA=1 on a cuDSS host")

@pytest.fixture
def cuda_ritz(monkeypatch):
    """Prove family GPU solves invoke the shared CUDA dense Ritz/SVD core."""
    import cupy as cp
    from fem.common import _gpu_resident_arnoldi as resident
    from fem.common import _gpu_arnoldi as arnoldi
    if not resident.is_available():
        pytest.skip('CUDA geev unavailable')
    original = resident.refined_gpu
    calls = []

    def track(hbar, shift, count):
        assert isinstance(hbar, cp.ndarray)
        result = original(hbar, shift, count)
        assert all(isinstance(array, cp.ndarray) for array in result)
        calls.append(hbar.shape)
        return result

    original_cpu = arnoldi.refined_coefficients
    def cpu_track(hbar, shift, count):
        assert isinstance(hbar, np.ndarray)
        calls.append(hbar.shape)
        return original_cpu(hbar, shift, count)
    monkeypatch.setattr(resident, 'refined_gpu', track)
    monkeypatch.setattr(arnoldi, 'refined_coefficients', cpu_track)
    yield calls
    assert any(shape[1] > 2 for shape in calls), 'The numerical solve did not invoke CUDA Ritz extraction'


@pytest.mark.parametrize("dim", [1, 2])
def test_electrostatic_cuda(dim):
    from fem.electrostatics import ElectrostaticSolver
    solver = ElectrostaticSolver(dim=dim, outer_potential=None)
    solver.set_potential(potential=0., geometry="left")
    solver.set_potential(potential=1., geometry="right")
    solver.mesh(max_element_size=.15)
    cpu = solver.solve(max_refinements=0, backend="cpu")
    gpu = solver.solve(max_refinements=0, backend="gpu")
    np.testing.assert_allclose(gpu.potential, cpu.potential, atol=1e-9)
    np.testing.assert_allclose(gpu.energy, cpu.energy, rtol=1e-8)
    assert gpu.residual_norm < 1e-9

@pytest.mark.parametrize("dim", [1, 2])
@pytest.mark.parametrize("backend", ["gpu", "hybrid", "auto"])
def test_waveguide_cuda(dim, backend, cuda_ritz):
    from fem.waveguide.solver_1d import ModeSolver1D
    from fem.waveguide.solver_2d import ModeSolver2D
    solver = ModeSolver1D(frequency=c, x_range=1.) if dim == 1 else ModeSolver2D(frequency=c, x_range=1., y_range=.5)
    solver.mesh(resolution=24 if dim == 1 else (5, 3))
    cpu = solver.solve(num_modes=1, neff_guess=.86, max_refinements=0, backend="cpu")
    gpu = solver.solve(num_modes=1, neff_guess=.86, max_refinements=0, backend=backend)
    np.testing.assert_allclose(gpu[0].neff, cpu[0].neff, rtol=1e-6, atol=1e-7)
    assert gpu[0].residual < 1e-7


@pytest.mark.parametrize("polarization", ["TE", "TM", "vector"])
def test_plane_wave_cudss(polarization):
    import fem
    outer = (-1/3, 1 + 1/3)
    settings = dict(frequency=c, x_range=outer, y_range=outer)
    if polarization == "vector":
        solver = fem.ScatteringSolver3D(**settings, z_range=outer, pml_layers=2)
        target = fem.shapes.Box(bounds=((.3, .7),)*3)
    else:
        solver = fem.ScatteringSolver2D(**settings, polarization=polarization)
        target = fem.shapes.Rectangle(bounds=((.3, .7),)*2)
    solver.add_geometry(shape=target, material=fem.materials.copper)
    solver.mesh(max_element_size=.23)
    cpu = solver.solve(backend="cpu", max_refinements=0)
    gpu = solver.solve(backend="gpu", max_refinements=0)
    assert gpu.solve_info["method"] == "cudss-direct"
    np.testing.assert_allclose(cpu.scattered, gpu.scattered, rtol=1e-7, atol=1e-8)
    np.testing.assert_allclose(cpu.radar_cross_section(-cpu.direction),
                               gpu.radar_cross_section(-gpu.direction), rtol=1e-7)
    if polarization == "vector":
        np.testing.assert_allclose(cpu.surface_current_vertices, gpu.surface_current_vertices,
                                   rtol=1e-7, atol=1e-9)


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("backend", ["gpu", "hybrid", "auto"])
def test_periodic_cuda(dim, backend, cuda_ritz):
    from fem.common import Material
    from fem.periodic import PeriodicModeSolver2D, PeriodicModeSolver3D

    options = dict(frequency=10e9, x_range=.02, z_range=.005,
                   background_material=Material(epsilon=2.25))
    solver = (PeriodicModeSolver2D(**options, polarization="TE") if dim == 2
              else PeriodicModeSolver3D(**options, y_range=.01))
    solver.mesh(max_element_size=.006, wavelength_elements=8)
    cpu = solver.solve(num_modes=1, neff_guess=1.3, max_refinements=0,
                       eigensolver="sparse" if dim == 3 else "dense", backend="cpu")
    gpu = solver.solve(num_modes=1, neff_guess=1.3, max_refinements=0,
                       eigensolver="sparse", backend=backend)
    np.testing.assert_allclose(gpu.neff, cpu.neff, rtol=1e-6, atol=1e-8)
    assert gpu.metadata["backend"] == ("gpu" if backend == "gpu" else "hybrid")
    assert gpu[0].residual < 1e-8
    # neff is invariant under arbitrary field scaling. Check the independently
    # implemented GPU power normalization using a phase-invariant field norm.
    assert np.linalg.norm(gpu[0].component("Ey")) == pytest.approx(
        np.linalg.norm(cpu[0].component("Ey")), rel=1e-5)


@pytest.mark.parametrize("impedance", [0., .3+.2j])
@pytest.mark.parametrize("order,half,size,layers,tolerance", [
    (1, .85 + 1/3, .10, 5, .12),
    (1, .85, .20, 4, .20),
    (2, .85, .20, 4, .03),
])
def test_metal_sphere_cudss_backscatter_against_mie(impedance, order, half, size, layers, tolerance):
    from scipy.special import spherical_jn, spherical_yn
    from fem.common.constants import ETA_0
    import fem

    radius = .15
    k = 2*np.pi
    n = np.arange(1, 10)
    x = k*radius
    psi = x*spherical_jn(n, x)
    dpsi = spherical_jn(n, x) + x*spherical_jn(n, x, derivative=True)
    xi = psi - 1j*x*spherical_yn(n, x)
    dxi = dpsi - 1j*(spherical_yn(n, x) + x*spherical_yn(n, x, derivative=True))
    a = (dpsi - 1j*impedance*psi) / (dxi - 1j*impedance*xi)
    b = (impedance*dpsi - 1j*psi) / (impedance*dxi - 1j*xi)
    expected = np.pi/k**2 * abs(np.sum((2*n+1)*(-1)**n*(a-b)))**2

    extent = (-half, half)
    solver = fem.ScatteringSolver3D(frequency=c, x_range=extent, y_range=extent, z_range=extent,
                                    pml_layers=layers)
    material = fem.SurfaceImpedance(impedance=impedance*ETA_0) if impedance else fem.materials.PEC
    solver.add_geometry(shape=fem.shapes.Sphere(center=(0,0,0), radius=radius), material=material)
    # Resolve propagation to the vacuum integration surface as well as the sphere.
    solver.mesh(max_element_size=size, element_order=order)
    result = solver.solve(backend="gpu", max_refinements=0)
    assert result.radar_cross_section(-result.direction)[0] == pytest.approx(expected, rel=tolerance)


def test_default_backend_uses_available_gpu(cuda_ritz):
    from fem.common.linalg import resolve_backend
    from fem.waveguide import ModeSolver1D
    assert resolve_backend() == "hybrid"
    solver = ModeSolver1D(frequency=c, x_range=1.)
    solver.mesh(resolution=24)
    result = solver.solve(num_modes=1, neff_guess=.86, max_refinements=0)
    assert result.backend == "hybrid"
    assert all(method == "cudss-hybrid-shift-invert" for method in result.metadata["methods"].values())


@pytest.mark.parametrize('case', ['pec_lossy', 'pec_sibc', 'pec_sibc_lossy'])
def test_waveguide_material_combinations_use_resident_core(case, cuda_ritz):
    from benchmarks.gpu_waveguide_2d import make_solver, solve
    solver, target, _ = make_solver(case, 8, 1)
    reference = solve(solver, target, 'cpu', 1)
    actual = solve(solver, target, 'gpu', 1)
    np.testing.assert_allclose(actual.neff, reference.neff, rtol=1e-7, atol=1e-8)
    assert actual[0].residual < 1e-8
    assert actual[0].divergence_residual < 1e-6


@pytest.mark.parametrize('backend', ['gpu', 'auto'])
def test_periodic_tm_cuda(backend, cuda_ritz):
    from fem import Material, PeriodicModeSolver2D
    solver = PeriodicModeSolver2D(frequency=10e9, x_range=.004, z_range=.005,
        background_material=Material(epsilon=2.25), polarization='TM')
    solver.mesh(max_element_size=.001)
    expected = solver.solve(backend='cpu', num_modes=1, neff_guess=1.5,
                            max_refinements=0, eigensolver='dense')
    result = solver.solve(backend=backend, num_modes=1, neff_guess=1.5, max_refinements=0)
    assert result.metadata['backend'] == ('gpu' if backend == 'gpu' else 'hybrid')
    np.testing.assert_allclose(result.neff, expected.neff, atol=1e-8, rtol=1e-7)
