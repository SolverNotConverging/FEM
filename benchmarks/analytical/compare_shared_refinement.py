"""Compare local residual rankings with measured CST refinement on shared tetrahedra."""

import argparse
import csv
from itertools import combinations, product
from pathlib import Path

import numpy as np
import meshio
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from scipy.stats import spearmanr
from skfem import MeshTet
from fem.adaptivity import bulk_mark

from cst_refinement_cases import CASES, PORT_LINE
from compare_two_port_meshes import cut_lines

TERMS = ("interface", "element_pde", "divergence", "natural_boundary", "port_boundary")
TRAIN = ("patch", "dipole", "iris")


def metal_edges(case, p, t):
    mesh = MeshTet(p.T, t.T, sort_t=True)
    boundary = mesh.boundary_facets()
    triangles = p[mesh.facets[:, boundary].T]
    axes = (0, 1, 2) if case in PORT_LINE else (2,)
    outer = np.zeros(len(boundary), bool)
    for axis in axes:
        for end in (p[:, axis].min(), p[:, axis].max()):
            outer |= np.max(abs(triangles[:, :, axis] - end), axis=1) < 1e-9
    return sharp_edge_cells(mesh, boundary[~outer])


def sharp_edge_cells(mesh, facets):
    vertices = mesh.p[:, mesh.facets[:, facets]].transpose(2, 1, 0)
    normals = np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    edges = mesh.f2e[:, facets].T.ravel()
    order = np.argsort(edges)
    ids = edges[order]
    # Each manifold surface edge has two incident conductor triangles.
    unique, first, count = np.unique(ids, return_index=True, return_counts=True)
    pairs = first[count == 2]
    dots = np.einsum("ij,ij->i", normals[order[pairs] // 3], normals[order[pairs + 1] // 3])
    sharp = unique[count == 2][abs(dots) < np.cos(np.deg2rad(35))]
    return np.isin(mesh.t2e, sharp).any(axis=0)


def regions(case, centers):
    x, y, z = centers.T * 1e3
    labels = np.full(len(x), "surroundings", dtype="U32")
    if case == "patch":
        labels[(abs(x) < 29) & (abs(y) < 29) & (z >= 0) & (z <= 1.6)] = "substrate"
        labels[(abs(x) < 16) & (abs(y) < 21) & (abs(z - 1.6) < 3)] = "patch vicinity"
        labels[(x >= -29) & (x < -4) & (abs(y) < 4) & (abs(z - 1.6) < 3)] = "feed and inset"
    elif case in ("dipole", "monopole"):
        radial = np.hypot(x, y) < 6
        labels[radial & (abs(z) < 30)] = "arm vicinity"
        labels[radial & (abs(abs(z) - 29) < 5)] = "arm tip"
        labels[radial & (abs(z) < 5)] = "feed gap"
    elif case == "iris":
        labels[abs(z) < 5] = "iris vicinity"
        labels[(abs(z) < 5) & (abs(abs(x) - 5.5) < 3)] = "aperture edge"
    else:
        labels[(z >= -2) & (z <= 4)] = "dielectric interior"
        labels[np.minimum(abs(z + 2), abs(z - 4)) < 1] = "dielectric interface"
    return labels


def average_precision(score, actual):
    # Group tied scores before computing a precision/recall step.
    order = np.argsort(-score)
    end = np.r_[np.flatnonzero(np.diff(score[order])), len(order) - 1]
    hits = np.cumsum(actual[order])[end]
    return (np.diff(hits, prepend=0) * hits / (end + 1)).sum() / actual.sum()


def prediction(terms, edge, weights):
    interface, pde, edge_weight = weights
    score = sum(
        value * ({"interface": interface, "element_pde": pde}.get(name, 1.0))
        for name, value in terms.items()
    )
    return np.maximum(score, 0) * np.where(edge, edge_weight, 1.0)


def write_rows(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("plots", type=Path)
    args = parser.parse_args()
    args.plots.mkdir(parents=True, exist_ok=True)
    samples = {}
    regions_out = []
    for case, step in product(CASES, (1, 2)):
        with np.load(args.data / f"{case}_{step}.npz") as data:
            p, t = data["points"], data["cells"]
            tetra = p[t]
            centers = tetra.mean(axis=1)
            longest = np.maximum.reduce(
                [
                    np.linalg.norm(tetra[:, a] - tetra[:, b], axis=1)
                    for a, b in combinations(range(4), 2)
                ]
            )
            edge = metal_edges(case, p, t)
            gain = data["log_gain"]
            actual = gain >= np.log2(1.5)
            terms = {key: data[key] for key in TERMS if key in data}
            samples[case, step] = (terms, edge, actual)
            labels = regions(case, centers)
            for region, selected in [(s, labels == s) for s in np.unique(labels)] + [
                ("sharp edge", edge),
                ("all", np.ones(len(edge), bool)),
            ]:
                if not selected.any():
                    continue
                off = data.get("log_gain_edge_off")
                regions_out.append(
                    dict(
                        case=case,
                        step=step,
                        region=region,
                        cells=selected.sum(),
                        cell_fraction=selected.mean(),
                        median_longest_edge_mm=np.median(longest[selected]) * 1000,
                        refined_fraction=actual[selected].mean(),
                        mean_log2_density_gain=gain[selected].mean(),
                        refined_fraction_edge_off=(off[selected] >= np.log2(1.5)).mean()
                        if off is not None
                        else "",
                    )
                )
    trial_rows = []
    # Equal weight per training geometry, with the simplest weights breaking exact ties.
    for weights in product((0.25, 1.0, 4.0), (0.25, 1.0, 4.0), (1.0, 2.0, 4.0)):
        aps = []
        for case in TRAIN:
            terms, edge, actual = samples[case, 1]
            aps.append(average_precision(prediction(terms, edge, weights), actual))
        trial_rows.append(
            dict(
                interface=weights[0],
                element_pde=weights[1],
                conductor_edge=weights[2],
                mean_train_ap=np.mean(aps),
            )
        )
    best = max(
        trial_rows,
        key=lambda row: (
            row["mean_train_ap"],
            -sum(abs(np.log2(row[key])) for key in ("interface", "element_pde", "conductor_edge")),
        ),
    )
    weights = tuple(best[key] for key in ("interface", "element_pde", "conductor_edge"))
    summary = []
    for (case, step), (terms, edge, actual) in samples.items():
        with np.load(args.data / f"{case}_{step}.npz") as data:
            score = prediction(terms, edge, weights)
            base = prediction(terms, edge, (1, 1, 1))
            marking = {}
            for policy, indicators in (("base", base), ("fitted", score)):
                marked = bulk_mark(indicators)
                marking[policy + "_marked"] = len(marked)
                marking[policy + "_precision"] = actual[marked].mean()
                marking[policy + "_recall"] = actual[marked].sum() / actual.sum()
            if case == "patch":
                meshio.write_points_cells(
                    args.data / f"patch_{step}.vtu",
                    data["points"],
                    [("tetra", data["cells"])],
                    cell_data={
                        "unit_indicator": [base],
                        "fitted_indicator": [score],
                        "observed_log2_density_gain": [data["log_gain"]],
                    },
                )
            summary.append(
                dict(
                    case=case,
                    step=step,
                    split="fit" if case in TRAIN and step == 1 else "validation",
                    prevalence=actual.mean(),
                    base_ap=average_precision(base, actual),
                    **marking,
                    fitted_ap=average_precision(score, actual),
                    base_spearman=spearmanr(base, data["log_gain"]).statistic,
                    fitted_spearman=spearmanr(score, data["log_gain"]).statistic,
                    fem_residual=float(data["residual"]),
                    fem_s11_real=data["fem_s11"].real.item(),
                    fem_s11_imag=data["fem_s11"].imag.item(),
                )
            )
    write_rows(args.data / "calibration.csv", trial_rows)
    write_rows(args.data / "predictions.csv", summary)
    write_rows(args.data / "regions.csv", regions_out)
    fig = Figure(figsize=(10, 4.5), layout="constrained")
    left, right = fig.subplots(1, 2)
    inds = np.arange(len(CASES))
    for j, step in enumerate((1, 2)):
        rows = [r for r in summary if r["step"] == step]
        shift = (j - 0.5) * 0.36
        left.bar(
            inds + shift,
            [r["base_ap"] for r in rows],
            width=0.32,
            color=("0.73", "#57839b")[j],
            label=f"Pass {step} baseline",
        )
        left.scatter(
            inds + shift,
            [r["fitted_ap"] for r in rows],
            marker="D",
            s=30,
            color="#9e4415",
            label="Fitted weights" if not j else None,
        )
        left.scatter(
            inds + shift,
            [r["prevalence"] for r in rows],
            marker="_",
            s=80,
            color="black",
            label="Random ranking" if not j else None,
        )
    rows = [r for r in regions_out if r["step"] == 1 and r["region"] == "sharp edge"]
    right.bar(
        inds - 0.18,
        [r["refined_fraction"] for r in rows],
        0.34,
        label="CST level 2",
        color="#176b71",
    )
    right.bar(
        inds + 0.18,
        [r["refined_fraction_edge_off"] for r in rows],
        0.34,
        label="CST off",
        color="#94bbb5",
    )
    for ax, title, ylabel in [
        (left, "Which starting cells get smaller?", "Average precision"),
        (right, "Cells incident on sharp conductor edges", "Fraction refined"),
    ]:
        ax.set(
            xticks=inds,
            xticklabels=["Patch", "Dipole", "Monopole", "Iris", "Slab"],
            ylim=(0, 1),
            title=title,
            ylabel=ylabel,
        )
        ax.legend(fontsize=8)
    fig.savefig(args.plots / "shared_mesh_residual.png", dpi=160)
    fig = Figure(figsize=(10, 3.5), layout="constrained")
    axes = fig.subplots(1, 3, sharex=True, sharey=True)
    for step, ax in enumerate(axes, 1):
        mesh = meshio.read(args.data / f"patch_{step}.vtu")
        # Swap y and z so the common mesh slicer cuts the patch in the xy plane.
        tetra = mesh.points[mesh.get_cells_type("tetra")][:, :, [0, 2, 1]] * 1000
        lines = cut_lines(tetra, plane=1.7)
        ax.add_collection(LineCollection(lines, colors="#17616c", linewidths=0.4))
        ax.set(
            xlim=(-32, 32), ylim=(-26, 26), aspect="equal", xlabel="x (mm)", title=f"Pass {step}"
        )
    axes[0].set_ylabel("y (mm)")
    fig.savefig(args.plots / "shared_mesh_patch.png", dpi=170)
    print("weights", weights)
    for r in summary:
        print(
            r["case"], r["step"], r["split"], "AP", round(r["base_ap"], 3), round(r["fitted_ap"], 3)
        )


if __name__ == "__main__":
    main()
