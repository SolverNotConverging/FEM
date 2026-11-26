"""Conforming Cartesian PML around a triangular or tetrahedral vacuum domain."""

import numpy as np
from skfem import Basis, BilinearForm, asm
from skfem.helpers import curl, dot, grad

def extend_box(mesh, ranges, profiles, layers):
    """Sweep the two faces of each axis in turn, including previous axes' extensions.

    Physical nodes and cells keep their indices. Sorting each base triangle's
    vertex IDs gives adjacent triangular prisms the same diagonal on shared faces.
    """
    tolerance = max(b - a for a, b in ranges) * 1e-9
    for axis, limits in enumerate(ranges):
        surface = mesh.facets[:, mesh.boundary_facets()]
        coords = mesh.p[axis, surface]
        points, cells = [mesh.p], [mesh.t]
        offset = mesh.nvertices
        for suffix, sign, limit in zip("-+", (-1, 1), limits):
            profile = profiles.get("xyz"[axis] + suffix)
            if profile is None:
                continue
            base = np.sort(surface[:, np.all(abs(coords - limit) < tolerance, axis=0)], axis=0)
            nodes, local = np.unique(base, return_inverse=True)
            local = local.reshape(base.shape)
            for layer in range(1, layers + 1):
                plane = mesh.p[:, nodes].copy()
                plane[axis] = limit + sign * profile.thickness * layer / layers
                points.append(plane)
                top = local + offset
                # Split an extruded edge into two triangles or a prism into three tets.
                for split in range(len(base)):
                    cells.append(np.vstack((base[: split + 1], top[split:])))
                base = top
                offset += len(nodes)
        mesh = type(mesh)(np.hstack(points), np.hstack(cells))
    return mesh


def material_diagonal(x, ranges, profiles):
    """Relative epsilon and mu for vacuum stretched along the selected exterior axes."""
    stretch = np.ones((3, *x.shape[1:]), dtype=complex)
    for side, profile in profiles.items():
        axis = "xyz".index(side[0])
        lo, hi = ranges[axis]
        depth = lo - x[axis] if side[1] == "-" else x[axis] - hi
        stretch[axis] *= profile.stretch(depth)
    return np.prod(stretch, axis=0) / stretch**2


def shell_matrix(basis, physical_cells, ranges, profiles, k):
    """Assemble the stretched vacuum on shell cells only."""
    mesh = basis.mesh
    if mesh.dim() == 3:
        @BilinearForm(dtype=complex)
        def form(u, v, w):
            return dot(curl(u) / w.tensor, curl(v)) - k**2 * dot(w.tensor * u, v)
    else:
        @BilinearForm(dtype=complex)
        def form(u, v, w):
            return dot(w.tensor[:2] * grad(u), grad(v)) - k**2 * w.tensor[2] * u * v

    # Bound temporary storage for the edge basis and stretched tensors at
    # every quadrature point, including when the physical mesh is refined.
    matrix = 0
    chunk_size = 8192
    for start in range(physical_cells, mesh.nelements, chunk_size):
        cells = np.arange(start, min(start + chunk_size, mesh.nelements))
        shell = Basis(mesh, basis.elem, elements=cells, dofs=basis.dofs, mapping=basis.mapping,
                      quadrature=basis.quadrature, disable_doflocs=True)
        tensor = material_diagonal(shell.global_coordinates(), ranges, profiles)
        matrix += asm(form, shell, tensor=tensor)
    return matrix
