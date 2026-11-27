"""Measure the two-port mesh study under FEM or the CST 2025 Python environment."""

import argparse
import cmath
import csv
from math import pi, sqrt
from pathlib import Path
import re


CASES = {"air": (.55, 1.), "short": (.275, 1.), "loaded": (.55, 1.44), "short_loaded": (.275, 1.44)}


def fem_passes(density, length, epsilon, fraction, factor, cad_point_size, mesh_dir, metric, backend):
    from two_port_mesh_geometry import fem_bore, write_bore
    from amr_observables import two_port_guide

    solver = two_port_guide(length, epsilon)
    solver.mesh(wavelength_elements=density, element_order=2, cad_point_size=cad_point_size)
    solver.configure_adaptivity(anisotropic=metric == "directional", marking_fraction=fraction, refinement_factor=factor,
                                refinement_backend=backend, marking_weights={})
    solver.configure_s_parameter_convergence(tolerance=.02, minimum_passes=3, consecutive_checks=1)
    counts = []

    def measure(result):
        points, cells = fem_bore(result.mesh_data, length)
        counts.append(len(cells))
        if mesh_dir:
            write_bore(mesh_dir / f'pass_{len(counts)}.vtu', points, cells)
        return result.s_parameters

    result = solver.solve(backend="cpu", max_refinements=7, pml_observable=measure)
    for entry, count in zip(result.solve_info["adaptive_history"], counts, strict=True):
        yield (entry["elements"], entry["field_dofs"], entry["s_parameter_matrix"].ravel(),
               entry["s_parameter_delta"], entry["status"], count)


def configure_cst(model, density, length, epsilon, max_passes):
    add = model.add_to_history
    add("PEC guide domain", '''ChangeSolverType "HF Frequency Domain"
Units.Geometry "mm"
Units.Frequency "GHz"
Solver.FrequencyRange "0.28", "0.299792458"
Boundary.Xmin "electric"
Boundary.Xmax "electric"
Boundary.Ymin "electric"
Boundary.Ymax "electric"
Boundary.Zmin "electric"
Boundary.Zmax "electric"
With Background
.XminSpace "0"
.XmaxSpace "0"
.YminSpace "0"
.YmaxSpace "0"
.ZminSpace "0"
.ZmaxSpace "0"
End With
''')
    material = "Vacuum"
    if epsilon != 1.:
        material = "Guide dielectric"
        add(material, f'''With Material
.Reset
.Name "{material}"
.Type "Normal"
.Epsilon "{epsilon}"
.Mue "1"
.Create
End With''')
    add("Uniform guide filling", f'''With Brick
.Reset
.Name "Guide"
.Component "component1"
.Material "{material}"
.Xrange "-360","360"
.Yrange "-150","150"
.Zrange "{-500*length}","{500*length}"
.Create
End With''')
    for number, face in enumerate(("zmin", "zmax"), 1):
        add(f"Port {number}", f'''With Port
.Reset
.PortNumber "{number}"
.NumberOfModes "1"
.Coordinates "Full"
.Orientation "{face}"
.PortOnBound "True"
.Create
End With''')
    add("Tetrahedral adaptation", f'''Mesh.MeshType "Tetrahedral"
With MeshSettings
.SetMeshType "Tet"
.Set "CellsPerWavelengthPolicy","cellsperwavelength"
.Set "StepsPerWaveNear","{density}"
.Set "StepsPerWaveFar","{density}"
.Set "StepsPerBoxNear","2"
.Set "ModelBoxDescrNear","maxedge"
End With
With FDSolver
.SetMethod "Tetrahedral","Discrete samples only"
.OrderTet "Second"
.Stimulation "All","All"
.MeshAdaptionTet "True"
.ResetSampleIntervals "all"
.AddSampleInterval "0.299792458","0.299792458","1","Single","True"
End With
With MeshAdaption3D
.SetType "HighFrequencyTet"
.MinPasses "{min(3, max_passes)}"
.MaxPasses "{max_passes}"
.MaxDeltaS "0.02"
.NumberOfDeltaSChecks "1"
End With''')


def run_cst(project, density, length, epsilon, max_passes):
    import cst.interface

    project = project.resolve()
    if project.suffix.lower() != ".cst" or project.exists() or project.with_suffix("").exists():
        raise ValueError("Choose a new .cst project filename.")
    project.parent.mkdir(parents=True, exist_ok=True)
    environment = cst.interface.DesignEnvironment.new()
    try:
        environment.set_quiet_mode(True)
        guide = environment.new_mws()
        configure_cst(guide.model3d, density, length, epsilon, max_passes)
        guide.save(str(project))
        guide.model3d.run_solver(timeout=240)
        guide.save()
        guide.close()
    finally:
        environment.close()


def cst_passes(project):
    import cst.results

    module = cst.results.ProjectFile(str(project.resolve())).get_3d()
    prefix, = [key.rsplit("\\", 1)[0] for key in module.get_tree_items()
               if "Adaptive Meshing" in key and key.endswith("\\Meshcells")]

    def series(name):
        item = module.get_result_item(prefix + "\\" + name)
        return dict(zip(item.get_xdata(), item.get_ydata()))

    cells = series("Meshcells")
    matrices = [series(f"S-Parameters\\S{i},{j}") for i in (1, 2) for j in (1, 2)]
    delta = series("S-Parameters\\Delta\\All S-Parameters")
    log = (project.with_suffix("") / "Result" / "Model.log").read_text(errors="replace")
    sections = re.split(r"Adaptive mesh refinement pass \d+", log)[1:]
    if len(sections) != len(cells):
        raise ValueError("Pass counts differ between the solver log and the result tree.")
    for (step, count), section in zip(sorted(cells.items()), sections):
        dofs = int(re.search(r"Degrees of freedom\s*:\s*(\d+)", section).group(1))
        status = "refine" if step != max(cells) else (
            "tolerance" if "desired accuracy limit is reached" in section else "stopped")
        yield int(count), dofs, [s[step] for s in matrices], delta.get(step), status, int(count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--density", type=float, required=True)
    parser.add_argument("--case", choices=CASES, default="air")
    parser.add_argument("--marking-fraction", type=float, default=.5, help="FEM bulk marking fraction.")
    parser.add_argument("--refinement-factor", type=float, default=2., help="FEM marked-cell size divisor.")
    parser.add_argument("--cst-max-passes", type=int, choices=(2, 8), default=8,
                        help="Set to two to save the intermediate CST mesh.")
    parser.add_argument("--cad-point-size", type=float, help="FEM CAD vertex target in metres.")
    parser.add_argument("--refinement-backend", choices=("auto", "gmsh"), default="auto")
    parser.add_argument("--metric", choices=("scalar", "directional"), default="directional")
    parser.add_argument("--mesh-dir", type=Path, help="Export FEM bore tetrahedra at each pass.")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--new-cst-project", type=Path, help="Create and run a new project in licensed CST.")
    source.add_argument("--cst-project", type=Path, help="Read a completed guide run with cst.results.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    project = args.new_cst_project or args.cst_project
    length, epsilon = CASES[args.case]
    if args.new_cst_project:
        run_cst(project, args.density, length, epsilon, args.cst_max_passes)
    rows = list(cst_passes(project) if project else fem_passes(
        args.density, length, epsilon, args.marking_fraction, args.refinement_factor, args.cad_point_size, args.mesh_dir, args.metric, args.refinement_backend))
    exact = cmath.exp(-1j * length * sqrt(epsilon*(2*pi)**2 - (pi/.72)**2))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["cells_per_wavelength", "pass", "physical_cells", "dofs", "max_delta_s",
                         "s11", "s12", "s21", "s22", "relative_s21_error", "status", "bore_cells",
                         "case", "marking_fraction", "refinement_factor", "cad_point_size", "metric", "refinement_backend"])
        previous = None
        for step, (cells, dofs, matrix, delta, status, bore) in enumerate(rows, 1):
            if previous is not None:
                observed = max(abs(a-b) for a, b in zip(matrix, previous))
                if delta is None or abs(observed-delta) > 1e-7:
                    raise ValueError(f"Reported delta {delta} differs from S-matrix change {observed}.")
            writer.writerow([args.density, step, cells, dofs, delta, *matrix,
                             abs(matrix[2]-exact)/abs(exact), status, bore, args.case,
                             "" if project else args.marking_fraction, "" if project else args.refinement_factor,
                             "" if project else args.cad_point_size, "" if project else args.metric,
                             "" if project else args.refinement_backend])
            previous = matrix


if __name__ == "__main__":
    main()
