"""Plot measured error and cost for each AMR observable from one or more CSVs."""

import argparse
import csv
from pathlib import Path

from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator


def plot(inputs, output):
    rows = []
    for path in inputs:
        with path.open(newline="", encoding="utf-8") as stream:
            rows.extend(csv.DictReader(stream))
    output.mkdir(parents=True, exist_ok=True)
    titles = {
        "narrow_gap": "Eccentric capacitor · charge per volt",
        "inclusion": "Dielectric cylinder · angular scattering width",
        "periodic_defect": "Periodic dielectric insert · Bloch index",
        "leaky_mode": "Stretched guide · complex effective index",
        "multiport": "Two-port guide · complex transmission",
        "plate_1d": "1D charged interval · electric energy",
        "plate_3d": "3D charged cube · electric energy",
    }
    styles = {"local": ("Local", "o", "#0072B2"), "uniform": ("Uniform", "s", "#D55E00"),
              "anisotropic": ("Directional MMG", "^", "#009E73"), "gmsh": ("Gmsh control", "D", "#8964A1")}
    for case in sorted({row["case"] for row in rows}):
        case_rows = [r for r in rows if r["case"] == case]
        initial = [r for r in case_rows if r["strategy"] == "initial"]
        fig = Figure(figsize=(8.4, 4.2), layout="constrained")
        error, cost = fig.subplots(1, 2)
        for strategy, (label, marker, color) in styles.items():
            values = sorted((r for r in case_rows if r["strategy"] == strategy), key=lambda r: int(r["level"]))
            if not values:
                continue
            series = initial + values
            x = [int(r["field_dofs"]) for r in series]
            error.semilogy(x, [100*float(r["relative_error"]) for r in series], marker=marker, color=color, label=label)
            cost.plot(x, [float(r["elapsed_s"]) for r in series], marker=marker, color=color)
        for ax in (error, cost):
            ax.set_xlabel("Field degrees of freedom")
            ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
            ax.grid(alpha=.22)
        error.set_ylabel("Relative observable error (%)")
        cost.set_ylabel("Mesh and solve time (s)")
        cost.set_ylim(bottom=0)
        error.legend(fontsize=8)
        fig.suptitle(titles[case])
        fig.savefig(output / f"{case}.png", dpi=150)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    plot(args.csv, args.output)
