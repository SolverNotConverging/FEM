"""Solve on the exported CST physical tetrahedra and measure FEM residual components."""

import argparse
import csv
from pathlib import Path

import numpy as np
import meshio
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
from skfem import Basis, BilinearForm, FacetBasis, LinearForm, MeshTet, asm
from skfem.helpers import curl, dot

from cst_refinement_cases import CASES, FREQUENCY_GHZ, PORT_LINE, boxes
from fem.adaptivity.driven import driven_residual
from fem.common.constants import C_0, ETA_0
from fem.common.pml import PMLSpec
from fem.common.radiation_pml import extend_box, shell_matrix
from fem.common.tet_hcurl import tetrahedral_element
from two_port_mesh_geometry import read_cst_tetrahedra


def material(case, centers):
    eps = np.ones(len(centers), complex)
    for _, mat, limits in boxes(case):
        bound = np.array(limits).reshape(3, 2) * 0.001
        inside = ((centers > bound[:, 0]) & (centers < bound[:, 1])).all(axis=1)
        if mat == "PEC" and inside.any():
            raise ValueError("A CST field tetrahedron lies in the conductor.")
        if mat in ("FR4", "Insert"):
            eps[inside] = 4.4 * (1 - 0.02j) if mat == "FR4" else 2.25
    return eps


def edge_voltage(basis, endpoints):
    # The integral of the first Whitney function along its oriented edge is
    # one. The second, odd function and all face bubbles integrate to zero.
    p, edges = basis.mesh.p.T, basis.mesh.edges.T
    a, b = np.array(endpoints) * 0.001
    on_line = (
        (abs(p[:, :2] - a[:2]) < 1e-9).all(axis=1)
        & (p[:, 2] >= a[2] - 1e-9)
        & (p[:, 2] <= b[2] + 1e-9)
    )
    chain = np.flatnonzero(on_line[edges].all(axis=1))
    dz = p[edges[chain, 1], 2] - p[edges[chain, 0], 2]
    np.testing.assert_allclose(abs(dz).sum(), b[2] - a[2], rtol=1e-7, atol=1e-12)
    v = np.zeros(basis.N)
    v[basis.dofs.edge_dofs[0, chain]] = np.sign(dz)
    return v


def te10(face):
    a, b = 0.02286, 0.01016

    def shape(x):
        return np.sqrt(2 / (a * b)) * np.cos(np.pi * x[0] / a)

    @LinearForm
    def projection(v, w):
        return v[1] * shape(w.x)

    profile = np.zeros_like(face.global_coordinates())
    profile[1] = shape(face.global_coordinates())
    return asm(projection, face), profile


def solve(case, p, t):
    mesh = MeshTet(p.T, t.T, sort_t=True)
    centers = p[t].mean(axis=1)
    eps = material(case, centers)
    antenna = case in PORT_LINE
    k = 2 * np.pi * FREQUENCY_GHZ[case] * 1e9 / C_0
    bounds = np.column_stack((p.min(axis=0), p.max(axis=0)))
    if antenna:
        profile = PMLSpec(thickness=2 * np.pi / (3 * k))
        profiles = dict.fromkeys(profile.sides("xyz"), profile)
        mesh = extend_box(mesh, bounds, profiles, layers=2)
    basis = Basis(
        mesh, tetrahedral_element(2), intorder=5, elements=np.arange(len(t)), disable_doflocs=True
    )

    @BilinearForm(dtype=complex)
    def maxwell(u, v, w):
        return dot(curl(u), curl(v)) - k * k * eps[:, None] * dot(u, v)

    matrix = asm(maxwell, basis)
    boundary = mesh.boundary_facets()
    traces = []
    if antenna:
        matrix += shell_matrix(basis, len(t), bounds, profiles, k)
        v = edge_voltage(basis, PORT_LINE[case])
        coupling = 1j * k * ETA_0 / 50
        rhs = 2 * coupling * v
        g = csc_matrix(v[:, None])
        matrix += coupling * (g @ g.T)
        reflection = v
    else:
        z = mesh.p[2, mesh.facets[:, boundary]]
        beta = np.sqrt(k * k - (np.pi / 0.02286) ** 2)
        ports = [boundary[np.all(abs(z - z0) < 1e-9, axis=0)] for z0 in bounds[2]]
        boundary = np.setdiff1d(boundary, np.concatenate(ports))
        for i, facets in enumerate(ports):
            face = FacetBasis(mesh, basis.elem, facets=facets, dofs=basis.dofs, intorder=5)
            v, field = te10(face)
            g = csc_matrix(v[:, None])
            matrix += 1j * beta * (g @ g.T)
            if i == 0:
                rhs = 2j * beta * v
                reflection = v
            traces.append((face, v, field))
    fixed = basis.dofs.get_facet_dofs(boundary).all()
    free = np.setdiff1d(np.arange(basis.N), fixed)
    coeff = np.zeros(basis.N, complex)
    coeff[free] = spsolve(matrix[free][:, free].tocsc(), rhs[free])
    if not np.isfinite(coeff).all():
        raise ValueError("Nonfinite shared-mesh solution.")
    port_traces = [
        (face, 1j * beta / k * ((v @ coeff) - 2 * (i == 0)) * field)
        for i, (face, v, field) in enumerate(traces)
    ]
    terms = {}
    _, residual, _ = driven_residual(
        basis,
        coeff,
        eps,
        np.ones(len(t)),
        k,
        5,
        port_traces=port_traces,
        component_indicators=terms,
    )
    return terms, residual, reflection @ coeff - 1


def volumes(tetra):
    return np.abs(np.linalg.det(tetra[:, 1:] - tetra[:, :1])) / 6


def locate(tetra, query):
    tree = cKDTree(tetra.mean(axis=1))
    inv = np.linalg.inv(tetra[:, 1:] - tetra[:, :1])
    result = np.full(len(query), -1)
    # Test barycentric containment, expanding the centroid search for thin cells.
    for start in range(0, len(query), 256):
        ids = np.arange(start, min(start + 256, len(query)))
        n = min(12, len(tetra))
        while len(ids):
            near = tree.query(query[ids], k=n)[1].reshape(len(ids), n)
            local = np.einsum("nki,nkij->nkj", query[ids, None] - tetra[near, 0], inv[near])
            inside = (local.min(axis=-1) >= -1e-8) & (local.sum(axis=-1) <= 1 + 1e-8)
            hit = inside.any(axis=1)
            result[ids[hit]] = near[hit, np.argmax(inside[hit], axis=1)]
            ids = ids[~hit]
            if len(ids) and n == len(tetra):
                raise ValueError("A sample in the source mesh is outside the target mesh.")
            n = min(2 * n, len(tetra))
    return result


def log_density_gain(before, after):
    # Symmetric, interior tetrahedron quadrature positions, away from interfaces.
    bary = np.full((4, 4), 0.1381966011250105)
    np.fill_diagonal(bary, 0.5854101966249685)
    sites = np.einsum("ij,njk->nik", bary, before)
    cells = locate(after, sites.reshape(-1, 3)).reshape(-1, 4)
    return np.log2(volumes(before)) - np.log2(volumes(after)[cells]).mean(axis=1)


def load_case(root, case):
    histories = []
    for folder in ("pass2", "pass3", "edge_off"):
        with (root / folder / (case + ".csv")).open() as f:
            histories.append(list(csv.DictReader(f)))
    short, full, edge_off = histories
    for a, b in zip(short, full[:2], strict=True):
        if a["tetrahedra"] != b["tetrahedra"]:
            raise ValueError("Repeated CST runs gave different pass counts.")
        for key in ("s11_real", "s11_imag"):
            np.testing.assert_allclose(float(a[key]), float(b[key]), atol=1e-7, rtol=0)

    def read(folder, file, count):
        p, t = read_cst_tetrahedra(root / folder / case / "Result" / file)
        if (len(t) != int(count) or not np.isfinite(p).all() or t.min() < 0
                or t.max() >= len(p) or not (volumes(p[t]) > 0).all()):
            raise ValueError("Invalid mesh or adaptive-result count mismatch.")
        return p, t

    p1, t1 = read("pass2", "3do.slim", short[0]["tetrahedra"])
    for folder, start in (("pass3", full[0]), ("edge_off", edge_off[0])):
        p, t = read(folder, "3do.slim", start["tetrahedra"])
        np.testing.assert_allclose(p, p1, rtol=0, atol=1e-12)
        np.testing.assert_array_equal(t, t1)
    p2, t2 = read("pass2", "3d.slim", short[1]["tetrahedra"])
    p3, t3 = read("pass3", "3d.slim", full[2]["tetrahedra"])
    po, to = read("edge_off", "3d.slim", edge_off[1]["tetrahedra"])
    volume = volumes(p1[t1]).sum()
    for tetra in (p2[t2], p3[t3], po[to]):
        np.testing.assert_allclose(volumes(tetra).sum(), volume, rtol=1e-9, atol=1e-15)
    return [(p1, t1), (p2, t2), (p3, t3)], po[to], histories


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("projects", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cases", choices=CASES, nargs="+", default=CASES)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for case in args.cases:
        meshes, off, histories = load_case(args.projects, case)
        p3, t3 = meshes[-1]
        meshio.write_points_cells(args.output / f"{case}_3.vtu", p3, [("tetra", t3)])
        for j, ((p, t), (q, u)) in enumerate(zip(meshes[:-1], meshes[1:], strict=True)):
            terms, residual, s11 = solve(case, p, t)
            gain = log_density_gain(p[t], q[u])
            extra = {"log_gain_edge_off": log_density_gain(p[t], off)} if j == 0 else {}
            np.savez_compressed(
                args.output / f"{case}_{j + 1}.npz",
                points=p,
                cells=t,
                log_gain=gain,
                residual=residual,
                fem_s11=s11,
                **terms,
                **extra,
            )
            print(case, j + 1, "residual", residual, "FEM S11", s11, flush=True)
        with (args.output / f"{case}_passes.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["singular_edge_level", "pass", "tetrahedra", "s11_real", "s11_imag"])
            for level, rows in [(2, histories[1]), (0, histories[2])]:
                for row in rows:
                    w.writerow(
                        [level, row["pass"], row["tetrahedra"], row["s11_real"], row["s11_imag"]]
                    )


if __name__ == "__main__":
    main()
