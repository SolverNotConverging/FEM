"""Model-dependent defaults for driven meshes, calibrated on the CST comparison cases."""

from dataclasses import dataclass

import numpy as np

from fem.common import materials, shapes


@dataclass(frozen=True)
class DrivenMeshPolicy:
    name: str = "general"
    density: float | None = None
    size_cad_points: bool = False
    local_sizes: tuple = ()
    fitted_weights: bool = False
    fraction: float = 0.5
    factor: float = 2.0
    backend: str = "auto"


def _bounds(shape):
    return np.asarray(shape.bounds).reshape(3, 2)


def _guide_policy(regions, ports, wavelength):
    if len(ports) != 2 or not all(
        isinstance(getattr(p, "cross_section", None), shapes.Rectangle) for p in ports
    ):
        return None
    left, right = ports
    if left.axis != right.axis or left.position == right.position:
        return None
    axis = "xyz".index(left.axis)
    bore = _bounds(left).copy()
    bore[axis] = sorted((left.position, right.position))
    transverse = np.arange(3) != axis
    lengths = bore[:, 1] - bore[:, 0]
    tol = 1e-8 * lengths.max()
    if not np.allclose(bore[transverse], _bounds(right)[transverse], rtol=0, atol=tol):
        return None
    metals = [s for s, m in regions if not isinstance(m, materials.Material)]
    walls = [
        s
        for s, m in regions
        if not isinstance(m, materials.Material)
        and isinstance(s, shapes.Difference)
        and isinstance(s.tool, shapes.Box)
        and isinstance(s.shape, shapes.Box)
        and np.allclose(_bounds(s.tool), bore, atol=tol, rtol=0)
        and np.all(_bounds(s.shape)[transverse, 0] < bore[transverse, 0] - tol)
        and np.all(_bounds(s.shape)[transverse, 1] > bore[transverse, 1] + tol)
        and _bounds(s.shape)[axis, 0] <= bore[axis, 0] + tol
        and _bounds(s.shape)[axis, 1] >= bore[axis, 1] - tol
    ]

    def spans_wall(shape, normal, side):
        b = _bounds(shape)
        along = np.arange(3) != normal
        return (
            isinstance(shape, shapes.Box)
            and abs(b[normal, 1 - side] - bore[normal, side]) <= tol
            and np.all(b[along, 0] <= bore[along, 0] + tol)
            and np.all(b[along, 1] >= bore[along, 1] - tol)
        )

    if not walls and not all(
        any(spans_wall(s, normal, side) for s in metals)
        for normal in np.flatnonzero(transverse)
        for side in (0, 1)
    ):
        return None
    index = 1.0
    for shape, material in regions:
        if any(shape is wall for wall in walls):
            continue
        bounds = _bounds(shape)
        if np.any(
            np.minimum(bounds[:, 1], bore[:, 1]) - np.maximum(bounds[:, 0], bore[:, 0]) <= tol
        ):
            continue
        if (
            isinstance(material, materials.Material)
            and isinstance(shape, shapes.Box)
            and np.allclose(bounds, bore, rtol=0, atol=tol)
        ):
            eps, mu = materials.bulk_values(material, form="scalar")
            calibrated_fill = eps.real > 0 and np.isclose(eps.imag, 0) and np.isclose(mu, 1)
            index = np.sqrt(eps.real) if calibrated_fill else None
        else:
            return DrivenMeshPolicy("guide discontinuity", density=6.0, fitted_weights=True)
    # The coarse uniform-guide fits were measured in the single-mode band.
    narrow, broad = sorted(lengths[transverse])
    if index is None:
        return DrivenMeshPolicy("uniform guide")
    local_wave = wavelength / index
    single_mode = max(broad, 2 * narrow) < local_wave < 2 * broad
    if (
        not single_mode
        or any(port.mode > 1 for port in ports)
        or not 1.0 - 1e-8 <= index <= 1.2 + 1e-8
    ):
        return DrivenMeshPolicy("uniform guide")
    # Nearest of the two measured length / broad-wall ratios, .275/.72 and .55/.72.
    short = lengths[axis] / max(lengths[transverse]) < (0.275 + 0.55) / (2 * 0.72)
    loaded = index > 1.1
    density, fraction, factor, backend = {
        (False, False): (2.5, 0.5, 2.25, "gmsh"),
        (True, False): (2.75, 0.5, 1.7, "gmsh"),
        (False, True): (3.25, 0.12, 1.1, "auto"),
        (True, True): (3.0, 0.5, 2.1, "gmsh"),
    }[short, loaded]
    return DrivenMeshPolicy(
        ("short " if short else "") + ("dielectric guide" if loaded else "air guide"),
        density,
        True,
        fraction=fraction,
        factor=factor,
        backend=backend,
    )


def driven_mesh_policy(regions, ports, dimension, wavelength, order):
    if dimension != 3 or order != 2:
        return DrivenMeshPolicy()
    guide = _guide_policy(regions, ports, wavelength)
    if guide is not None:
        return guide
    metals = [_bounds(s) for s, m in regions if not isinstance(m, materials.Material)]
    boards = []
    printed = False
    for shape, material in regions:
        if not isinstance(shape, shapes.Box) or not isinstance(material, materials.Material):
            continue
        bounds = _bounds(shape)
        spans = np.diff(bounds, axis=1).ravel()
        axis = np.argmin(spans)
        if np.partition(spans, 1)[1] < 4 * spans[axis]:
            continue
        eps, mu = materials.bulk_values(material, form="scalar")
        if eps.real <= 0 or mu.real <= 0:
            continue
        index = np.sqrt(abs(eps * mu))
        if index <= 1:
            continue
        lateral = np.arange(3) != axis
        for metal in metals:
            extent = metal[:, 1] - metal[:, 0]
            overlap = np.minimum(metal[:, 1], bounds[:, 1]) - np.maximum(metal[:, 0], bounds[:, 0])
            if (
                np.min(extent[lateral]) >= 4 * extent[axis]
                and np.all(overlap[lateral] > 0)
                and overlap[axis] >= -1e-8 * spans.max()
            ):
                boards.append((shape, wavelength / (10 * index)))
                printed = True
                break
        else:
            if not ports:
                boards.append((shape, wavelength / (6 * index)))
    if boards:
        return DrivenMeshPolicy(
            "printed conductor" if printed else "dielectric slab",
            local_sizes=tuple(boards),
            fitted_weights=True,
        )
    if len(ports) == 1 and hasattr(ports[0], "gap"):
        axis = "xyz".index(ports[0].axis)
        gap = _bounds(ports[0])
        rods, ground_sides = [], set()
        for metal in metals:
            spans = np.diff(metal, axis=1).ravel()
            transverse = np.arange(3) != axis
            # Count only electrodes overlapping the terminal footprint and touching a gap end.
            overlap = np.minimum(gap[:, 1], metal[:, 1]) - np.maximum(gap[:, 0], metal[:, 0])
            distances = abs(metal[axis, ::-1] - gap[axis])
            side = int(np.argmin(distances))
            if np.any(overlap[transverse] <= 0) or distances[side] > 1e-8 * spans.max():
                continue
            if spans[axis] >= 4 * max(spans[transverse]):
                rods.append((side, shapes.Box(bounds=metal)))
            elif min(spans[transverse]) >= 4 * spans[axis]:
                ground_sides.add(side)
        monopole = len(rods) == 1 and 1 - rods[0][0] in ground_sides
        if monopole or (len(rods) == 2 and rods[0][0] != rods[1][0]):
            return DrivenMeshPolicy(
                "monopole" if monopole else "dipole",
                local_sizes=tuple((s, wavelength / 12) for _, s in rods),
                fitted_weights=monopole,
            )
    return DrivenMeshPolicy()
