"""Create small HF tetrahedral models for the shared-mesh refinement benchmark (CST Python)."""

import argparse
import csv
from pathlib import Path

CASES = ("patch", "dipole", "monopole", "iris", "dielectric_step")
FREQUENCY_GHZ = dict(patch=2.45, dipole=2.45, monopole=2.45, iris=10.0, dielectric_step=10.0)
# Coordinates in mm. Solid metal is excluded from the field domain.
PORT_LINE = {
    "patch": ((-28.0, 0.0, 0.0), (-28.0, 0.0, 1.6)),
    "dipole": ((0.0, 0.0, -0.5), (0.0, 0.0, 0.5)),
    "monopole": ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
}


def boxes(case):
    if case == "patch":
        # A 5 mm inset leaves two 1 mm slots beside the 3 mm feed strip.
        return [
            ("ground", "PEC", (-29, 29, -29, 29, -0.2, 0)),
            ("board", "FR4", (-29, 29, -29, 29, 0, 1.6)),
            ("patch", "PEC", (-4.8, 13.8, -18.5, 18.5, 1.6, 1.8)),
            ("wing1", "PEC", (-13.8, -4.8, -18.5, -2.5, 1.6, 1.8)),
            ("wing2", "PEC", (-13.8, -4.8, 2.5, 18.5, 1.6, 1.8)),
            ("feed", "PEC", (-29, -4.8, -1.5, 1.5, 1.6, 1.8)),
        ]
    if case == "dipole":
        return [
            ("arm1", "PEC", (-1, 1, -1, 1, -29, -0.5)),
            ("arm2", "PEC", (-1, 1, -1, 1, 0.5, 29)),
        ]
    if case == "monopole":
        return [
            ("ground", "PEC", (-29, 29, -29, 29, -0.5, 0)),
            ("arm", "PEC", (-1, 1, -1, 1, 1, 30)),
        ]
    guide = [("guide", "Vacuum", (-11.43, 11.43, -5.08, 5.08, -20, 20))]
    if case == "iris":
        return guide + [
            ("iris1", "PEC", (-11.43, -5.5, -5.08, 5.08, -1, 1)),
            ("iris2", "PEC", (5.5, 11.43, -5.08, 5.08, -1, 1)),
        ]
    return guide + [("insert", "Insert", (-11.43, 11.43, -5.08, 5.08, -2, 4))]


def configure(model, case, passes, edge_level, wavelength_prerefine):
    add = model.add_to_history
    antenna = case in PORT_LINE
    frequency = FREQUENCY_GHZ[case]
    add(
        "Frequency domain",
        f'ChangeSolverType "HF Frequency Domain"\nUnits.Geometry "mm"\nUnits.Frequency "GHz"\nSolver.FrequencyRange "{frequency * 0.96}", "{frequency}"',
    )
    for side in ("Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax"):
        bc = "open" if antenna else "electric"
        add(side, f'Boundary.{side} "{bc}"\nBackground.{side}Space "{25 if antenna else 0}"')
    if case in ("patch", "dielectric_step"):
        name, eps, tand = ("FR4", 4.4, 0.02) if antenna else ("Insert", 2.25, 0.0)
        add(
            name,
            f'With Material\n.Reset\n.Name "{name}"\n.Type "Normal"\n.Epsilon "{eps}"\n.Mue "1"\n.TanD "{tand}"\n.TanDFreq "{frequency}"\n.Create\nEnd With',
        )
    for name, mat, (x0, x1, y0, y1, z0, z1) in boxes(case):
        add(
            name,
            f'With Brick\n.Reset\n.Name "{name}"\n.Component "model"\n.Material "{mat}"\n.Xrange "{x0}","{x1}"\n.Yrange "{y0}","{y1}"\n.Zrange "{z0}","{z1}"\n.Create\nEnd With',
        )
    if case == "patch":
        for solid in ("wing1", "wing2", "feed"):
            add("Join patch metal", f'Solid.Add "model:patch", "model:{solid}"')
    if antenna:
        a, b = PORT_LINE[case]
        endpoints = "\n".join(
            f'.SetP{i} "False","{pt[0]}","{pt[1]}","{pt[2]}"' for i, pt in enumerate((a, b), 1)
        )
        add(
            "50 ohm feed",
            f'With DiscretePort\n.Reset\n.PortNumber "1"\n.Type "SParameter"\n.Impedance "50"\n{endpoints}\n.Create\nEnd With',
        )
    else:
        for i, side in enumerate(("zmin", "zmax"), 1):
            add(
                "TE10 port",
                f'With Port\n.Reset\n.PortNumber "{i}"\n.NumberOfModes "1"\n.Coordinates "Full"\n.Orientation "{side}"\n.PortOnBound "True"\n.Create\nEnd With',
            )
    settings = f'''Mesh.MeshType "Tetrahedral"
With FDSolver
.SetMethod "Tetrahedral","Discrete samples only"
.OrderTet "Second"
.Stimulation "1","1"
.MeshAdaptionTet "True"
.ResetSampleIntervals "all"
.AddSampleInterval "{frequency}","{frequency}","1","Single","True"
End With
With MeshAdaption3D
.SetType "HighFrequencyTet"
.MinPasses "{passes}"
.MaxPasses "{passes}"
.MaxDeltaS "0.02"
.NumberOfDeltaSChecks "1"
End With'''
    add("Adaptive frequency sample", settings)
    if not wavelength_prerefine:
        add("Fixed initial tetrahedra", 'MeshAdaption3D.WavelengthBasedRefinement "False"')
    if edge_level is not None:
        add("Singular edge refinement", f'MeshAdaption3D.SingularEdgeRefinement "{edge_level}"')


def measure(project, output):
    import cst.results

    module = cst.results.ProjectFile(str(project)).get_3d()
    tree = module.get_tree_items()
    (prefix,) = [
        key.rsplit("\\", 1)[0]
        for key in tree
        if "Adaptive Meshing" in key and key.endswith("\\Meshcells")
    ]

    def series(key):
        item = module.get_result_item(prefix + "\\" + key)
        return dict(zip(item.get_xdata(), item.get_ydata()))

    cells = series("Meshcells")
    reflection = series("S-Parameters\\S1,1")
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["pass", "tetrahedra", "s11_real", "s11_imag"])
        for step in sorted(cells):
            writer.writerow(
                [int(step), int(cells[step]), reflection[step].real, reflection[step].imag]
            )


def main():
    import cst.interface

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cases", choices=CASES, nargs="+", default=CASES)
    parser.add_argument("--passes", type=int, choices=(2, 3), default=2)
    parser.add_argument("--no-wavelength-prerefine", action="store_true")
    parser.add_argument("--edge-level", type=int, choices=range(4))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    environment = cst.interface.DesignEnvironment.new()
    try:
        environment.set_quiet_mode(True)
        for case in args.cases:
            project = (args.output / (case + ".cst")).resolve()
            if project.exists() or project.with_suffix("").exists():
                raise FileExistsError(project)
            model = environment.new_mws()
            try:
                configure(
                    model.model3d,
                    case,
                    args.passes,
                    args.edge_level,
                    not args.no_wavelength_prerefine,
                )
                model.save(str(project))
                model.model3d.run_solver(timeout=360)
                model.save()
            finally:
                model.close()
            measure(project, project.with_suffix(".csv"))
            print(case, project, flush=True)
    finally:
        environment.close()


if __name__ == "__main__":
    main()
