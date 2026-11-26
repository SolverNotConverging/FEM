"""Time CPU and cuDSS plane-wave solves on the same copper-target mesh."""
from argparse import ArgumentParser
from statistics import median
from time import perf_counter

import numpy as np

import fem


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--size", type=float, default=.015, help="Mesh target size in metres")
    parser.add_argument("--element-order", choices=(1, 2), type=int, default=2)
    settings = parser.parse_args()
    if settings.repeats < 1:
        parser.error("--repeats must be positive")
    for pol in ("TE", "TM", "vector"):
        options = dict(frequency=3e9, x_range=(-.06, .06), y_range=(-.06, .06))
        if pol == "vector":
            solver = fem.ScatteringSolver3D(**options, z_range=(-.06, .06))
            target = fem.shapes.Sphere(center=(0, 0, 0), radius=.015)
        else:
            solver = fem.ScatteringSolver2D(**options, polarization=pol)
            target = fem.shapes.Circle(center=(0, 0), radius=.015)
        solver.add_geometry(shape=target, material=fem.materials.copper)
        solver.mesh(max_element_size=settings.size, element_order=settings.element_order)
        times = {key: [] for key in ("cpu", "gpu")}
        results = {}
        for rep in range(settings.repeats + 1):
            for backend in (("cpu", "gpu") if rep % 2 == 0 else ("gpu", "cpu")):
                start = perf_counter()
                results[backend] = solver.solve(backend=backend, max_refinements=0)
                if rep:
                    times[backend].append(perf_counter() - start)
        cpu, gpu = results["cpu"], results["gpu"]
        assert gpu.solve_info["method"] == "cudss-direct"
        np.testing.assert_allclose(gpu.scattered, cpu.scattered, rtol=1e-7, atol=1e-8)
        np.testing.assert_allclose(gpu.radar_cross_section(-gpu.direction),
                                   cpu.radar_cross_section(-cpu.direction), rtol=1e-7)
        print(f"{pol}: {cpu.metadata['ndofs']} DOFs, CPU {median(times['cpu']):.3f} s, "
              f"cuDSS {median(times['gpu']):.3f} s")


if __name__ == "__main__":
    main()
