"""Compare count and spatial distribution of cells in exported waveguide bores."""
import argparse
import csv
from itertools import combinations
from pathlib import Path

import meshio
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from two_port_mesh_comparison import CASES


def geometry(path):
    mesh = meshio.read(path)
    return mesh.points[mesh.get_cells_type("tetra")]


def distribution(vertices, length):
    centers = vertices.mean(axis=1)
    bounds = [(-.36, .36), (-.15, .15), (-length/2, length/2)]
    histogram, _ = np.histogramdd(centers, bins=(3, 2, 3), range=bounds)
    return histogram / len(centers)


def differences(trial, reference, length):
    count_error = abs(np.log(len(trial)/len(reference)))
    spatial_tv = abs(distribution(trial,length) - distribution(reference,length)).sum() / 2
    return count_error, spatial_tv


def cut_lines(vertices, plane=.001):
    distance = vertices[:, :, 1] - plane
    crosses = (distance.min(axis=1) < 0) & (distance.max(axis=1) > 0)
    loops = []
    for tetra, dy in zip(vertices[crosses], distance[crosses]):
        points = []
        for i, j in combinations(range(4), 2):
            if dy[i] * dy[j] < 0:
                xyz = tetra[i] - dy[i] / (dy[j] - dy[i]) * (tetra[j] - tetra[i])
                points.append(xyz[[0, 2]])
            elif dy[i] == 0:
                points.append(tetra[i, [0, 2]])
            elif dy[j] == 0:
                points.append(tetra[j, [0, 2]])
        polygon = np.unique(points, axis=0)
        offset = polygon - polygon.mean(axis=0)
        ordered = polygon[np.argsort(np.arctan2(offset[:, 1], offset[:, 0]))]
        loops.append(np.vstack([ordered, ordered[0]]))
    return loops


def compare_case(directory, case, length, writer, plot):
    variants = {"cst": "CST", "baseline": "FEM baseline", "tuned": "FEM material fit"}
    if (directory / case / "geometry_fit").is_dir():
        variants["geometry_fit"] = "FEM geometry refit"
    figure = Figure(figsize=(3 * len(variants), 7), layout="constrained") if plot else None
    axes = figure.subplots(3, len(variants), sharex=True, sharey=True) if plot else None
    for step in (1, 2, 3):
        ref = geometry(directory / case / "cst" / f"pass_{step}.vtu")
        for column, (solver, label) in enumerate(variants.items()):
            tet = ref if solver == "cst" else geometry(directory / case / solver / f"pass_{step}.vtu")
            count_error, tv = differences(tet, ref, length)
            centers = tet.mean(axis=1)
            volume = abs(np.linalg.det(tet[:, 1:] - tet[:, :1])) / 6
            writer.writerow([case, step, solver, len(tet), len(ref), count_error, tv,
                             np.mean(abs(centers[:, 2]) > length/3), np.median(volume),
                             *distribution(tet, length).ravel()])
            if plot:
                ax = axes[step-1, column]
                ax.add_collection(LineCollection(cut_lines(tet), linewidths=.53, colors="#225467"))
                ax.set(xlim=(-.36,.36), ylim=(-length/2,length/2), aspect="equal")
                if step == 1:
                    ax.set_title(label)
                ax.set_xlabel(f"{len(tet):,} bore tetrahedra")
                ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
                if column == 0:
                    ax.set_ylabel(f"Pass {step}")
    if plot:
        figure.savefig(plot, dpi=160)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("--plot-dir", type=Path, required=True)
    args = parser.parse_args()
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    with (args.data/"spatial_comparison.csv").open("w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["case","pass","solver","bore_cells","cst_bore_cells","abs_log_count_ratio",
                         "spatial_tv","end_sixths_fraction","median_volume_m3",
                         *(f"fraction_x{x}_y{y}_z{z}" for x in range(3) for y in range(2) for z in range(3))])
        for case,(length,_) in CASES.items():
            plot = args.plot_dir / f"two_port_{case}_mesh.png" if case in ("air","short_loaded") else None
            compare_case(args.data/"meshes", case, length, writer, plot)


if __name__ == "__main__":
    main()
