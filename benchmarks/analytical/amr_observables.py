"""Compare local AMR and a sequence of uniform size targets against analytical observables."""

import argparse
import csv
import multiprocessing
import sys
from contextlib import ExitStack
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.special import h2vp, hankel2, jv, jvp
import fem


def narrow_gap(h, solve, anisotropic=False):
    solver = fem.ElectrostaticSolver(dim=2, x_range=(-.5, .5), y_range=(-.5, .5))
    circle = fem.shapes.Circle
    r, outer, offset = .17, .38, .18
    ring = fem.shapes.Difference(shape=circle(center=(0., 0.), radius=.43),
                                 tool=circle(center=(0., 0.), radius=outer))
    solver.set_potential(geometry=ring, potential=0., name="shield")
    solver.set_potential(geometry=circle(center=(offset, 0.), radius=r), potential=1., name="wire")
    solver.mesh(max_element_size=h, geometry_tolerance=h**2 / 3)
    solver.configure_adaptivity(anisotropic=anisotropic)
    result = solver.solve(**solve)
    from scipy.constants import epsilon_0
    capacitance = 2*np.pi*epsilon_0 / np.arccosh((outer**2 + r**2 - offset**2) / (2*outer*r))
    return result, "capacitance_F_per_m", result.conductor_charge("wire"), capacitance


def inclusion(h, solve, anisotropic=False):
    radius, index = .14, 3.
    solver = fem.ScatteringSolver2D(frequency=fem.C_0, x_range=(-1., 1.), y_range=(-1., 1.),
                                    pml_layers=6, polarization="TM")
    solver.add_geometry(shape=fem.shapes.Circle(center=(0., 0.), radius=radius),
                        material=fem.Material(epsilon=index**2))
    solver.mesh(max_element_size=h)
    solver.configure_adaptivity(anisotropic=anisotropic)
    result = solver.solve(**solve)
    angles = np.linspace(0., np.pi, 19)
    exact = cylinder_width(radius, index, angles)
    measured = result.radar_cross_section(np.column_stack((np.cos(angles), np.sin(angles))))
    return result, "angular_scattering_width_m", measured, exact


def cylinder_width(radius, index, angles):
    """TM cylindrical-wave series at unit free-space wavelength."""
    k, n = 2*np.pi, np.arange(-18, 19)
    x = k*radius
    alpha = index * jvp(n, index*x) / jv(n, index*x)
    coeff = -(jvp(n, x) - alpha*jv(n, x)) / (h2vp(n, x) - alpha*hankel2(n, x))
    return 4/k * abs(np.exp(1j*angles[:, None]*n) @ coeff)**2


def periodic_defect(h, solve, anisotropic=False):
    period, insert = .2, .04
    solver = fem.PeriodicModeSolver2D(frequency=fem.C_0, x_range=(0., 1.), z_range=(0., period),
                                      polarization="TE")
    solver.add_geometry(shape=fem.shapes.Rectangle(bounds=((0., 1.), (.08, .08+insert))),
                        material=fem.Material(epsilon=4.))
    # TE plate mode in each piece. The unit-cell transfer matrix gives cos(beta * period).
    q, p = np.pi*np.sqrt(3), np.pi*np.sqrt(15)
    a, b = q*(period-insert), p*insert
    trace = np.cos(a)*np.cos(b) - .5*(q/p+p/q)*np.sin(a)*np.sin(b)
    exact = np.arccos(complex(trace)) / (2*np.pi*period)
    solver.mesh(max_element_size=h)
    result = solver.solve(num_modes=1, neff_guess=exact, **solve)
    return result, "bloch_effective_index", result[0].neff, exact


def leaky_mode(h, solve, anisotropic=False):
    solver = fem.PeriodicModeSolver2D(frequency=fem.C_0, x_range=(0., 1.), z_range=(0., .15),
                                      polarization="TE")
    thickness, sigma, order = .1, 3., 3
    solver.add_pml(thickness=thickness, direction="x", sigma_max=sigma, order=order)
    separation = 1 - 2j * thickness * sigma / (order+1)
    exact = np.sqrt(1-(.5/separation)**2)
    solver.mesh(max_element_size=h, pml_layers=max(1, round(.14/h)))
    result = solver.solve(num_modes=1, neff_guess=exact, **solve)
    return result, "complex_effective_index", result[0].neff, exact


def two_port_guide(length=.55, epsilon=1.):
    solver = fem.AntennaSolver3D(frequency=fem.C_0, x_range=(-1., 1.), y_range=(-.85, .85),
                                 z_range=(-.95, .95), pml_layers=1)
    u, v, ends = (-.36, .36), (-.15, .15), (-length/2, length/2)
    bore = fem.shapes.Box(bounds=(u, v, ends))
    outer = fem.shapes.Box(bounds=((-.45, .45), (-.27, .27), (ends[0]-.095, ends[1]+.095)))
    solver.add_geometry(shape=fem.shapes.Difference(shape=outer, tool=bore), material=fem.materials.PEC)
    if epsilon != 1.:
        solver.add_geometry(shape=bore, material=fem.Material(epsilon=epsilon, name="Guide dielectric"))
    for name, z in zip(("input", "output"), ends):
        solver.add_waveguide_port(name=name, axis="z", position=z, num_modes=1,
                                  cross_section=fem.shapes.Rectangle(bounds=(u, v)))
    return solver


def multiport(h, solve, anisotropic=False, element_order=1, wavelength_elements=None):
    solver = two_port_guide()
    solver.mesh(max_element_size=h, element_order=element_order, wavelength_elements=wavelength_elements)
    solver.configure_adaptivity(anisotropic=anisotropic)
    result = solver.solve(**solve)
    beta = np.sqrt((2*np.pi)**2 - (np.pi/.72)**2)
    return result, "complex_S21", result.s_parameters[1, 0], np.exp(-1j*beta*.55)


def charged_plate(h, solve, anisotropic, dim):
    from scipy.constants import epsilon_0
    transverse = dict(y_range=1., z_range=1.) if dim == 3 else {}
    solver = fem.ElectrostaticSolver(dim=dim, x_range=1., outer_potential=None, **transverse)
    solver.set_potential(geometry="left" if dim == 1 else "x-", potential=0., name="left")
    solver.set_potential(geometry="right" if dim == 1 else "x+", potential=0., name="right")
    shape = fem.shapes.Interval(bounds=(0., 1.)) if dim == 1 else fem.shapes.Box(bounds=((0., 1.),)*3)
    solver.add_charge_density(geometry=shape, density=epsilon_0)
    solver.mesh(max_element_size=h, element_order=1, boundary_refinement=None)
    solver.configure_adaptivity(anisotropic=anisotropic)
    # V=x(1-x)/2. Integrate epsilon_0*|grad V|^2/2 over unit measure.
    result = solver.solve(**solve)
    return result, "electrostatic_energy", result.energy, epsilon_0/24


def plate_1d(h, solve, anisotropic=False):
    return charged_plate(h, solve, anisotropic, 1)


def plate_3d(h, solve, anisotropic=False):
    return charged_plate(h, solve, anisotropic, 3)


CASES = {f.__name__: (f, h) for f, h in
         [(narrow_gap, .12), (inclusion, .1), (periodic_defect, .16), (leaky_mode, .14),
          (multiport, .19), (plate_1d, .25), (plate_3d, .30)]}

MMG_CASES = {"narrow_gap", "inclusion", "multiport", "plate_3d"}


def peak_rss_mib():
    if sys.platform != "win32":
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / (1024**2 if sys.platform == "darwin" else 1024)
    import ctypes
    from ctypes import wintypes
    class MemoryCounters(ctypes.Structure):
        _fields_ = [("size", wintypes.DWORD), ("faults", wintypes.DWORD),
                    ("peak_rss", ctypes.c_size_t), ("other_counters", ctypes.c_size_t * 7)]
    counters = MemoryCounters()
    counters.size = ctypes.sizeof(counters)
    query = ctypes.windll.psapi.GetProcessMemoryInfo
    query.argtypes = [wintypes.HANDLE, ctypes.POINTER(MemoryCounters), wintypes.DWORD]
    if not query(wintypes.HANDLE(-1), ctypes.byref(counters), counters.size):
        raise ctypes.WinError()
    return counters.peak_rss / 1024**2


def measure(case, strategy, level):
    build, initial_h = CASES[case]
    h = initial_h / 1.5**level if strategy == "uniform" else initial_h
    started = perf_counter()
    local = strategy in ("local", "anisotropic", "gmsh")
    with ExitStack() as controls:
        if strategy == "gmsh":
            # Select the existing Gmsh candidate builder directly. Forcing an
            # MMG failure would include a discarded candidate in the timing.
            from unittest.mock import patch
            def gmsh_candidate(owner, dim, weights, build, record, **options):
                record["refinement_backend"] = "gmsh"
                return build()
            for module in ("remeshing", "policy"):
                controls.enter_context(patch(f"fem.adaptivity.{module}.adapt_candidate", gmsh_candidate))
        result, observable, measured, exact = build(h, dict(
            backend="cpu", max_refinements=level if local else 0, adaptive_tolerance=1e-10),
            anisotropic=strategy == "anisotropic")
    elapsed = perf_counter() - started
    steps = result.solve_info["adaptive_history"]
    history = steps[-1]
    dofs = history.get("field_dofs")
    if dofs is None:
        # Electrostatics uses scalar nodal elements, including in 3D.
        coefficients = result.potential_coefficients
        dofs = len(result.potential if coefficients is None else coefficients)
    error = float(np.linalg.norm(np.asarray(measured) - exact) / np.linalg.norm(exact))
    return dict(case=case, strategy=strategy, level=level, observable=observable,
                relative_error=error, field_dofs=int(dofs), residual=history["residual"],
                elapsed_s=elapsed, peak_rss_MiB=peak_rss_mib(), elements=history["elements"],
                refinements=len(steps)-1, backends=" / ".join(
                    step["refinement_backend"] for step in steps if "refinement_backend" in step),
                fallback_reason=" / ".join(step["fallback_reason"] for step in steps if "fallback_reason" in step),
                status=history["status"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--levels", nargs="+", type=int, default=[0, 1, 2])
    controls = parser.add_mutually_exclusive_group()
    controls.add_argument("--anisotropic", action="store_true",
                          help="Also measure directional MMG refinement for nonperiodic cases.")
    controls.add_argument("--gmsh-control", action="store_true",
                          help="Run only the Gmsh local control for positive levels of nonperiodic cases.")
    parser.add_argument("--output", type=Path, default=Path("outputs/benchmarks/analytical/amr_observables.csv"))
    args = parser.parse_args()
    if min(args.levels) < 0:
        parser.error("levels must be nonnegative")
    if args.gmsh_control and (not MMG_CASES.intersection(args.cases) or max(args.levels) == 0):
        parser.error("the Gmsh control needs a nonperiodic 2D/3D case and a positive refinement level")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # A fresh process makes peak resident memory comparable between runs.
    with multiprocessing.get_context("spawn").Pool(1, maxtasksperchild=1) as pool, args.output.open("w", newline="") as stream:
        writer = None
        for case in args.cases:
            for level in args.levels:
                strategies = ["local", "uniform"] if level else ["initial"]
                if level and args.anisotropic and case in MMG_CASES:
                    strategies.append("anisotropic")
                if args.gmsh_control:
                    strategies = ["gmsh"] if level and case in MMG_CASES else []
                for strategy in strategies:
                    row = pool.apply(measure, (case, strategy, level))
                    if writer is None:
                        writer = csv.DictWriter(stream, fieldnames=row)
                        writer.writeheader()
                    writer.writerow(row)
                    stream.flush()
                    print(row, flush=True)


if __name__ == "__main__":
    main()
