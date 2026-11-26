"""Differentiate cell-local polynomial reconstructions on affine simplices."""

from itertools import combinations_with_replacement

import numpy as np


def polynomial_fit(x):
    dim = len(x)
    terms = [(), *((a,) for a in range(dim)), *combinations_with_replacement(range(dim), 2)]
    values = np.array([np.prod(x[list(term)], axis=0) for term in terms])
    if np.linalg.matrix_rank(values) < len(terms):
        terms = terms[:dim + 1]
        values = values[:dim + 1]
    return terms, np.linalg.pinv(values.T)


def material_from_volume(basis, material):
    """Recover one-sided coefficients when a caller supplies only a point evaluator."""
    coordinates = basis.global_coordinates()
    values = material(*coordinates)
    terms, inverse = polynomial_fit(basis.X)
    fits, diagonals = [], []
    for value in values:
        value = np.asarray(value)
        if value.shape == (3,):
            value = value[:, None, None]
        diagonal = np.broadcast_to(value, (3, *coordinates.shape[1:]))
        diagonals.append(diagonal)
        fits.append(diagonal @ inverse.T)

    def evaluate(cells, *coords):
        if isinstance(cells, slice):
            return tuple(diagonals)
        local = basis.mapping.invF(np.asarray(coords), tind=cells)
        monomials = np.array([np.prod(local[list(term)], axis=0) for term in terms])
        return tuple(np.einsum("vct,tcq->vcq", fit[:, cells], monomials) for fit in fits)

    return evaluate


def derivative_operators(basis):
    """Reference differentiation matrices and inverse affine cell Jacobians.

    Fit quadratic samples when unisolvent, otherwise use a linear fit. This
    differentiates FE polynomials exactly up to the reconstruction degree.
    Stretched coefficients in PMLs use the same local reconstruction.
    """
    cached = getattr(basis, "_fem_derivative_operators", None)
    if cached is not None:
        return cached
    x = basis.X
    dim, count = x.shape
    terms, inverse = polynomial_fit(x)
    matrices = []
    for axis in range(dim):
        slopes = np.zeros((count, len(terms)))
        for j, term in enumerate(terms):
            for i, a in enumerate(term):
                if a == axis:
                    slopes[:, j] += np.prod(x[list(term[:i] + term[i+1:])], axis=0)
        matrices.append(slopes @ inverse)
    jac = basis.mapping.invDF(x[:, :1], tind=basis.tind)[..., 0]
    basis._fem_derivative_operators = np.asarray(matrices), jac
    return basis._fem_derivative_operators


def partial(field, axis, operators):
    diff, inverse = operators
    return sum((field @ matrix.T) * inverse[r, axis, :, None] for r, matrix in enumerate(diff))


def modal_volume(dfield, magnetic, neff, axes, operators, xp=np):
    """Squared curl(H) - i D and div(D) with the Bloch phase derivative."""
    curl = -1j * dfield.copy()
    shift = -1j * neff
    curl[0] -= shift * magnetic[1]
    curl[1] += shift * magnetic[0]
    divergence = shift * dfield[2]
    for local, axis in enumerate(axes):
        b, c = (axis + 1) % 3, (axis + 2) % 3
        curl[c] += partial(magnetic[b], local, operators)
        curl[b] -= partial(magnetic[c], local, operators)
        divergence += partial(dfield[axis], local, operators)
    return xp.sum(abs(curl)**2, axis=0), abs(divergence)**2


def domain_residuals(indicators, residual, pml_cells, xp=np):
    """Partition one normalized squared indicator over PML and physical cells."""
    squared = indicators.sum()
    fraction = xp.clip(indicators[pml_cells].sum() / xp.maximum(squared, np.finfo(float).tiny), 0., 1.)
    return residual * xp.sqrt(1 - fraction), residual * xp.sqrt(fraction)
