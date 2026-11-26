"""Flux-jump indicators and conforming local refinement for electrostatics."""

from dataclasses import replace
from itertools import combinations
from math import factorial

import numpy as np

from .meshing import FEMMesh, _mesh_extrema


def flux_indicators(mesh: FEMMesh, displacement: np.ndarray, charge: np.ndarray,
                    prescribed_facets: np.ndarray, corner_displacement=None) -> tuple[np.ndarray, float]:
    """Estimate discretization error from div(D)=rho and normal-flux jumps.

    The relative indicator is normalized by the elementwise L2 norm of D.
    It is a refinement heuristic, not a certified bound on the potential.
    Prescribed-potential faces carry conductor charge and are excluded.
    """
    points = mesh.nodes[mesh.elements]
    native = mesh.mesh
    if mesh.nodes.shape[1] == 1:
        measures = np.abs(points[:, 1, 0] - points[:, 0, 0])
        h = measures
    else:
        measures = np.abs(np.linalg.det(points[:, 1:] - points[:, :1])) / factorial(native.dim())
        h = np.maximum.reduce([
            np.linalg.norm(points[:, i] - points[:, j], axis=1)
            for i, j in combinations(range(native.dim() + 1), 2)
        ])
    if corner_displacement is not None:
        gradients = np.linalg.solve(points[:, 1:] - points[:, :1],
                                    corner_displacement[:, 1:] - corner_displacement[:, :1])
        charge = charge - np.trace(gradients, axis1=1, axis2=2)
    squared = h**2 * measures * charge**2
    adjacent = native.f2t
    first, second = adjacent
    facets = native.facets.T
    active = ~prescribed_facets
    interior = (second >= 0) & active
    boundary = (second < 0) & active
    if mesh.nodes.shape[1] == 1:
        normals = np.ones((len(facets), 1))
        face_measure = np.ones(len(facets))
        face_h = h[first].copy()
        face_h[interior] = 0.5 * (h[first[interior]] + h[second[interior]])
    elif native.dim() == 2:
        edges = mesh.nodes[facets[:, 1]] - mesh.nodes[facets[:, 0]]
        face_measure = np.linalg.norm(edges, axis=1)
        face_h = face_measure
        normals = np.column_stack((-edges[:, 1], edges[:, 0])) / face_measure[:, None]
    else:
        a, b, c = np.moveaxis(mesh.nodes[facets], 1, 0)
        cross = np.cross(b - a, c - a)
        area2 = np.linalg.norm(cross, axis=1)
        normals = cross / area2[:, None]
        face_measure = area2 / 2
        face_h = np.maximum.reduce([np.linalg.norm(v, axis=1) for v in (b-a, c-b, a-c)])

    def normal_flux(cells, selected):
        if corner_displacement is None:
            return np.sum(displacement[cells] * normals[selected], axis=-1)[:, None]
        local = (mesh.elements[cells, None, :] == facets[selected, :, None]).argmax(axis=2)
        trace = corner_displacement[cells[:, None], local]
        return np.sum(trace * normals[selected, None, :], axis=-1)

    left_flux = normal_flux(first, np.arange(len(first)))
    jumps = left_flux.copy()
    jumps[interior] -= normal_flux(second[interior], np.flatnonzero(interior))
    # Integrate the squared linear trace exactly on each simplex facet.
    n = jumps.shape[1]
    squared_jump = (np.sum(jumps**2, axis=1) + jumps.sum(axis=1)**2) / (n * (n + 1))
    contribution = face_h * face_measure * squared_jump
    np.add.at(squared, first[interior], 0.5 * contribution[interior])
    np.add.at(squared, second[interior], 0.5 * contribution[interior])
    np.add.at(squared, first[boundary], contribution[boundary])
    if corner_displacement is None:
        density = np.sum(displacement**2, axis=1)
    else:
        vertices = points.shape[1]
        density = (np.sum(corner_displacement**2, axis=(1, 2))
                   + np.sum(corner_displacement.sum(axis=1)**2, axis=1)) / (vertices * (vertices + 1))
    scale = float(measures @ density)
    relative = np.sqrt(float(np.sum(squared)) / max(scale, np.finfo(float).tiny))
    return squared, float(relative)


def refine_marked(mesh: FEMMesh, marked: np.ndarray) -> FEMMesh:
    """Bisect marked cells and close hanging nodes, preserving parent materials."""
    native = mesh.mesh.refined(marked)
    nodes = np.asarray(native.p.T, dtype=float)
    elements = np.asarray(native.t.T, dtype=np.int64)
    centers = nodes[elements].mean(axis=1)
    parents = mesh.mesh.element_finder()(*centers.T)
    minimum, maximum = _mesh_extrema(nodes, elements)
    return replace(
        mesh, mesh=native, nodes=nodes, elements=elements,
        element_tags=mesh.element_tags[parents],
        charge_density=(None if mesh.charge_density is None else mesh.charge_density[parents]),
        potential_cells=({owner: cells[parents] for owner, cells in mesh.potential_cells.items()}
                         if mesh.potential_cells is not None else None),
        info=replace(mesh.info, nodes=len(nodes), elements=len(elements),
                     minimum_edge=minimum, maximum_edge=maximum),
    )
