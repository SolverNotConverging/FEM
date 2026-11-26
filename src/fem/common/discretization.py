"""Numerical integration defaults shared by the finite-element solvers."""

from .errors import ConfigurationError
from .validation import positive_integer


# Degree five integrates quadratic mass forms and has positive simplex weights.
DEFAULT_QUADRATURE_ORDER = 5


def field_orders(element_order, quadrature_order, supported):
    element = positive_integer(element_order, "element_order")
    if element not in supported:
        raise ConfigurationError(f"Supported element_order values for this solver are {supported}.")
    quadrature = positive_integer(quadrature_order, "quadrature_order", minimum=2 * element)
    return element, quadrature
