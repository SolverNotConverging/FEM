"""Export the two-port guide bore as VTK, with coordinates in metres."""

import argparse
import csv
from pathlib import Path

import meshio
import numpy as np


def validate_bore(points, cells, length, expected_count=None):
    if points.ndim != 2 or points.shape[1] != 3 or cells.ndim != 2 or cells.shape[1] != 4:
        raise ValueError("Expected 3D vertices and tetrahedra.")
    if not cells.size or cells.min() < 0 or cells.max() >= len(points):
        raise ValueError("Invalid tetrahedron connectivity.")
    if expected_count is not None and len(cells) != expected_count:
        raise ValueError(f"Mesh has {len(cells)} cells, but the results report {expected_count}.")
    half = np.array([.36, .15, length/2])
    np.testing.assert_allclose(points.min(axis=0), -half, atol=1e-9)
    np.testing.assert_allclose(points.max(axis=0), half, atol=1e-9)
    tetra = points[cells]
    volumes = abs(np.linalg.det(tetra[:, 1:] - tetra[:, :1])) / 6
    if np.any(volumes <= 0) or not np.isclose(volumes.sum(), .72*.30*length, rtol=1e-8, atol=1e-12):
        raise ValueError("Tetrahedra do not fill the guide bore volume.")


def fem_bore(mesh, length):
    points, cells = mesh.coordinates, mesh.elements
    centers = points[cells].mean(axis=1)
    inside = (abs(centers) < [.36, .15, length/2]).all(axis=1)
    used, inverse = np.unique(cells[inside], return_inverse=True)
    points, cells = points[used], inverse.reshape(-1, 4)
    validate_bore(points, cells, length)
    return points, cells


def read_cst_tetrahedra(path):
    # CST writes element counts, named properties, scalar types and byte order in
    # the mesh header. Only the declared fixed-size records through tetrahedra are needed.
    types = {"float64": "<f8", "float32": "<f4", "int32": "<i4", "uint32": "<u4"}
    groups = []
    with path.open("rb") as stream:
        if stream.readline().strip() != b"SLIM" or not stream.readline().startswith(b"format binary_little_endian "):
            raise ValueError(f"Unsupported mesh format in {path}.")
        for _ in range(512):
            words = stream.readline(512).decode("ascii").split()
            if words == ["end"]:
                break
            if words[:1] == ["element"]:
                groups.append((words[1], int(words[-1]), []))
            elif words[:1] == ["property"]:
                groups[-1][2].append((words[2], types[words[1]]))
        else:
            raise ValueError("Missing mesh header terminator.")
        points = None
        for name, count, properties in groups:
            if not properties or count < 0:
                raise ValueError(f"Unexpected records before tetrahedra: {name}")
            records = np.fromfile(stream, dtype=properties, count=count)
            if len(records) != count:
                raise ValueError(f"Incomplete {name} records.")
            if name == "node":
                # The reference projects use millimetres.
                points = np.column_stack([records[key] for key in ("x","y","z")]) * .001
            elif name == "tetrahedron" and points is not None:
                cells = np.column_stack([records[f"node_index{i}"] for i in range(4)])
                return points, cells
    raise ValueError(f"No tetrahedral mesh in {path}.")


def cst_bore(path, length, expected_count):
    points, cells = read_cst_tetrahedra(path)
    validate_bore(points, cells, length, expected_count)
    return points, cells


def write_bore(path, points, cells):
    path.parent.mkdir(parents=True, exist_ok=True)
    meshio.write_points_cells(path, points, [("tetra", cells)], binary=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project", type=Path)
    parser.add_argument("measurements", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--length", type=float, default=.55, help="Guide length in metres.")
    args = parser.parse_args()
    with args.measurements.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    result = args.project.with_suffix("") / "Result"
    for filename, row in (("3do.slim", rows[0]), ("3d.slim", rows[-1])):
        mesh = cst_bore(result / filename, args.length, int(row["physical_cells"]))
        write_bore(args.output_dir / f'pass_{row["pass"]}.vtu', *mesh)


if __name__ == "__main__":
    main()
