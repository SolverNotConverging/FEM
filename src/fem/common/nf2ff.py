"""Identify facets and outward orientation on a conforming closed contour or surface."""

import numpy as np


def interface_facets(mesh, box):
    inside = box.contains(*mesh.p[:, mesh.t].mean(axis=1))
    left, right = mesh.f2t
    facets = np.flatnonzero((right >= 0) & (inside[left] != inside[right]))
    return facets, inside[left[facets]]
