"""Plot the four measured two-port refinement histories."""
import argparse
import csv
from pathlib import Path

from matplotlib.figure import Figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    fig = Figure(figsize=(9, 4.4), layout="constrained")
    cells, error = fig.subplots(1, 2)
    for solver, color in (("cst", "#006d77"), ("fem", "#a34714")):
        for density, style in ((4, "o-"), (8, "s--")):
            with (args.data / f"{solver}_{density}.csv").open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            passes = [int(row["pass"]) for row in rows]
            label = f"{solver.upper()}  λ/{density}"
            cells.semilogy(passes, [int(r["bore_cells"]) for r in rows], style, color=color, label=label)
            error.semilogy(passes, [100 * float(r["relative_s21_error"]) for r in rows], style, color=color)
    cells.set_title("Waveguide bore mesh")
    cells.set_ylabel("Tetrahedra inside the bore")
    error.set_title("Transmission accuracy")
    error.set_ylabel("Complex S₂₁ error against TE₁₀ (%)")
    for axis in (cells, error):
        axis.set(xlabel="Physical mesh pass", xticks=[1, 2, 3])
        axis.grid(alpha=.24, which="both")
    cells.legend(loc="upper left")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)


if __name__ == "__main__":
    main()
