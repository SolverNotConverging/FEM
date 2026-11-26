"""Fixed-frequency one/two-strip dielectric image-line pilot.

Lengths in the command line and configuration are in mm; FEM uses metres.
The bottom exterior PEC face is the ground plane. PML covers x-/x+/y+ only.
Widths w1 and w2 are along z; each strip spans the guide in x. The internal
gap is edge-to-edge; the remaining gap closes the cell across its z boundary.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
import gc
import hashlib
from importlib.metadata import PackageNotFoundError, version
from itertools import product
import json
import os
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
from scipy.spatial import cKDTree

MM = 1e-3
C0 = 299792458.0
ETA0 = 376.730313668
DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "outputs" / "image_line_pilot"


@dataclass(frozen=True)
class Config:
    frequency_ghz: float = 20.0
    period_mm: float = 10.0
    guide_width_mm: float = 6.0
    guide_height_mm: float = 2.0
    epsilon_r: float = 6.0
    metal_thickness_mm: float = 0.1
    air_mm: float = 5.0
    pml_mm: float = 3.0
    mesh_mm: float = 1.0
    wavelength_elements: int = 12
    backend: str = "gpu"
    modes: int = 6
    ncv: int = 96
    neff_guess: float = 1.7
    harmonic: int = -1
    min_overlap: float = 0.75
    max_neff_step: float = 0.3
    max_pml_fraction: float = 0.5
    min_core_fraction: float = 0.15

    @property
    def k0(self):
        return 2 * np.pi * self.frequency_ghz * 1e9 / C0

    def validate(self):
        for key in ("frequency_ghz", "period_mm", "guide_width_mm", "guide_height_mm",
                    "epsilon_r", "metal_thickness_mm", "air_mm", "pml_mm", "mesh_mm"):
            if not np.isfinite(getattr(self, key)) or getattr(self, key) <= 0:
                raise ValueError(f"{key} must be finite and positive")
        if self.epsilon_r <= 1:
            raise ValueError("epsilon_r must exceed the air permittivity")
        if self.wavelength_elements < 4 or self.modes < 1:
            raise ValueError("wavelength_elements >= 4 and modes >= 1 are required")
        if self.ncv <= max(4 * self.modes, self.modes + 12) + 1:
            raise ValueError("ncv is too small for the requested candidate modes")


@dataclass(frozen=True)
class Cell:
    strips: int
    w1_mm: float = 0.0
    w2_mm: float = 0.0
    gap_mm: float = 0.0

    @property
    def key(self):
        return f"s{self.strips}_w1_{self.w1_mm:g}_w2_{self.w2_mm:g}_g_{self.gap_mm:g}"

    def intervals(self, period_mm):
        if self.strips == 0:
            return []
        values = (self.w1_mm,) if self.strips == 1 else (self.w1_mm, self.w2_mm, self.gap_mm)
        if self.strips not in (1, 2) or any(not np.isfinite(v) or v <= 0 for v in values):
            raise ValueError("Strip widths and the two-strip internal gap must be positive")
        occupied = self.w1_mm if self.strips == 1 else sum(values)
        if occupied >= period_mm:
            raise ValueError("Strips and internal gap must leave a positive outer gap")
        start = (period_mm - occupied) / 2
        intervals = [(start, start + self.w1_mm)]
        if self.strips == 2:
            start = intervals[0][1] + self.gap_mm
            intervals.append((start, start + self.w2_mm))
        return intervals


def make_solver(config, cell):
    import fem

    intervals = cell.intervals(config.period_mm)
    half_x = config.guide_width_mm / 2 + config.air_mm + config.pml_mm
    top_y = config.guide_height_mm + config.metal_thickness_mm + config.air_mm + config.pml_mm
    solver = fem.PeriodicModeSolver3D(
        frequency=config.frequency_ghz * 1e9, x_range=(-half_x * MM, half_x * MM),
        y_range=(0, top_y * MM), z_range=(0, config.period_mm * MM),
        boundary=fem.materials.PEC)
    solver.add_box(name="image_guide", x_range=(-config.guide_width_mm * MM / 2,
                                               config.guide_width_mm * MM / 2),
                   y_range=(0, config.guide_height_mm * MM),
                   z_range=(0, config.period_mm * MM),
                   material=fem.Material(name="lossless dielectric", epsilon=config.epsilon_r))
    for index, (lower, upper) in enumerate(intervals):
        solver.add_box(name=f"strip_{index + 1}",
                       x_range=(-config.guide_width_mm * MM / 2, config.guide_width_mm * MM / 2),
                       y_range=(config.guide_height_mm * MM,
                                (config.guide_height_mm + config.metal_thickness_mm) * MM),
                       z_range=(lower * MM, upper * MM), material=fem.materials.PEC)
    solver.add_pml(thickness=config.pml_mm * MM, direction="x")
    solver.add_pml(thickness=config.pml_mm * MM, direction="y+")
    solver.mesh(max_element_size=config.mesh_mm * MM,
                wavelength_elements=config.wavelength_elements)
    return solver


def core_mask(points, config):
    return ((np.abs(points[:, 0]) < config.guide_width_mm * MM / 2)
            & (points[:, 1] > 0) & (points[:, 1] < config.guide_height_mm * MM))


def field_signature(mode, config):
    """Approximate overlap on a fixed interior grid, independent of mesh DOFs.

    Nearest barycentre sampling is intentionally restricted to dielectric cells.
    This is a tracking diagnostic, not a rigorous biorthogonal modal product.
    """
    mask = core_mask(mode.fields.coordinates, config)
    points = mode.fields.coordinates[mask]
    if not len(points):
        raise ValueError("No dielectric field samples")
    grid = np.stack(np.meshgrid(
        np.linspace(-0.4, 0.4, 7) * config.guide_width_mm * MM,
        np.linspace(0.15, 0.75, 5) * config.guide_height_mm * MM,
        np.linspace(0.05, 0.95, 11) * config.period_mm * MM, indexing="ij"), axis=-1).reshape(-1, 3)
    _, indices = cKDTree(points).query(grid)
    values = np.stack([mode.component(name)[mask][indices] * scale for name, scale in
                       (("Ex", np.sqrt(config.epsilon_r)), ("Ey", np.sqrt(config.epsilon_r)),
                        ("Ez", np.sqrt(config.epsilon_r)), ("Hx", ETA0), ("Hy", ETA0),
                        ("Hz", ETA0))], axis=1).ravel()
    return values / np.linalg.norm(values)


def overlap(first, second):
    return float(np.clip(abs(np.vdot(first, second)), 0, 1))


def describe_mode(mode, config):
    fields = mode.fields
    tetra = fields.mesh_points[fields.mesh_cells]
    volume = np.abs(np.linalg.det(tetra[:, 1:] - tetra[:, :1])) / 6
    owners = fields.sample_element_indices
    core = core_mask(fields.coordinates, config)
    pml = np.asarray(fields.metadata["cell_pml_fraction"]) > 0
    eps = np.where(core, config.epsilon_r, 1.0)
    density = eps * sum(abs(mode.component(n)) ** 2 for n in ("Ex", "Ey", "Ez"))
    density += ETA0**2 * sum(abs(mode.component(n)) ** 2 for n in ("Hx", "Hy", "Hz"))
    weighted = density * volume[owners]
    fraction = float(weighted[core].sum() / max(weighted[~pml].sum(), 1e-300))
    beta_h = mode.beta.real + 2 * np.pi * config.harmonic / (config.period_mm * MM)
    beta_norm = beta_h / config.k0
    return dict(neff_real=mode.neff.real, neff_imag=mode.neff.imag,
                alpha_per_m=mode.alpha, alpha_over_k0=mode.alpha / config.k0,
                beta_per_m=mode.beta.real, beta_harmonic_per_m=beta_h,
                beta_harmonic_over_k0=beta_norm,
                angle_deg=float(np.degrees(np.arcsin(beta_norm))) if abs(beta_norm) < 1 else None,
                direction=mode.direction, pml_fraction=mode.pml_fraction,
                core_fraction=fraction, residual=mode.residual, gauss_residual=mode.gauss_residual)


def candidate_valid(data, config, bound_reference=False):
    # The nonradiating seed has true alpha=0. Small signed discretization/PML
    # errors must not prevent tracking it, but are retained in the output.
    attenuation_ok = (abs(data["alpha_over_k0"]) < 1e-3 if bound_reference
                      else data["alpha_over_k0"] >= -1e-8)
    return (data["direction"] == "forward" and attenuation_ok
            and data["pml_fraction"] <= config.max_pml_fraction
            and data["core_fraction"] >= config.min_core_fraction)


def _solve_cell(config, cell, reference=None):
    started = perf_counter()
    solver = make_solver(config, cell)
    assembled = perf_counter()
    target = complex(config.neff_guess, -0.002) if reference is None else reference[0]
    modes = solver.solve(backend=config.backend, num_modes=config.modes, direction="all",
                         max_pml_fraction=None, neff_guess=target, max_refinements=0,
                         eigensolver_tolerance=1e-8, ncv=config.ncv)
    solved = perf_counter()
    candidates, signatures = [], []
    for mode in modes:
        data = describe_mode(mode, config)
        signature = field_signature(mode, config)
        data["overlap"] = None if reference is None else overlap(reference[1], signature)
        data["neff_step"] = None if reference is None else abs(mode.neff - reference[0])
        candidates.append(data)
        signatures.append(signature)
    eligible = [i for i, data in enumerate(candidates)
                if candidate_valid(data, config, bound_reference=cell.strips == 0)]
    if reference is None:
        # Seed the bound guide mode; an air/PML root cannot anchor the sweep.
        eligible = [i for i in eligible if 1 < candidates[i]["neff_real"] < np.sqrt(config.epsilon_r)]
        eligible.sort(key=lambda i: candidates[i]["core_fraction"], reverse=True)
    else:
        eligible = [i for i in eligible if candidates[i]["overlap"] >= config.min_overlap
                    and candidates[i]["neff_step"] <= config.max_neff_step]
        eligible.sort(key=lambda i: candidates[i]["overlap"] - 0.1 * candidates[i]["neff_step"],
                      reverse=True)
    record = dict(case=cell.key, **asdict(cell), status="untracked", candidates=candidates,
                  mesh_assembly_seconds=assembled-started, solve_seconds=solved-assembled,
                  dofs=solver.system.ndofs, elements=solver.mesh_data.mesh.nelements,
                  backend=modes.metadata["backend"])
    if not eligible:
        record["error"] = "No candidate passed the physical/branch tracking criteria"
        record["total_seconds"] = perf_counter()-started
        return record, None, modes
    selected = eligible[0]
    data = candidates[selected]
    record.update(data, selected_mode=selected, status="accepted")
    if len(eligible) > 1 and reference is not None:
        margin = data["overlap"] - candidates[eligible[1]]["overlap"]
        record["overlap_margin"] = margin
        if margin < 0.05:
            record.update(status="ambiguous", error="Two physical candidates have similar field overlap")
    record["total_seconds"] = perf_counter()-started
    tracked = (modes[selected].neff, signatures[selected]) if record["status"] == "accepted" else None
    return record, tracked, modes


def solve_cell(config, cell, reference=None):
    """Solve a cell and release resident arrays before constructing another mesh.

    The scene listeners and resident session contain cycles. Python's normal
    object-count GC does not account for their large CUDA allocations.
    Different meshes also leave differently sized blocks in CuPy's pool.
    """
    started = perf_counter()
    try:
        record, tracked, modes = _solve_cell(config, cell, reference)
    finally:
        gc.collect()
        cp = sys.modules.get("cupy")
        if cp is not None and config.backend != "cpu":
            cp.get_default_memory_pool().free_all_blocks()
    record["total_seconds"] = perf_counter()-started
    return record, tracked, modes


def sweep_cells(widths, gaps, single_widths):
    singles = [Cell(1, float(w)) for w in sorted(set(single_widths))]
    pairs = [Cell(2, float(w1), float(w2), float(gap))
             for gap, w2, w1 in product(sorted(set(gaps)), sorted(set(widths)), sorted(set(widths)))]
    cells = singles + pairs
    if len({cell.key for cell in cells}) != len(cells):
        raise ValueError("Sample values produce duplicate case IDs; use more widely spaced widths/gaps")
    return cells


def environment_info():
    packages = {}
    for name in ("electromagnetics-fem", "numpy", "scipy", "gmsh", "scikit-fem",
                 "cupy-cuda12x", "cupy-cuda13x", "nvidia-cudss-cu12", "nvidia-cudss-cu13"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            pass
    info = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=sys.version,
                platform=platform.platform(), packages=packages, command=sys.argv,
                pilot_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    try:
        import cupy as cp
        props = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)
        name = props["name"]
        info["gpu"] = dict(name=name.decode() if isinstance(name, bytes) else name,
                           memory_bytes=props["totalGlobalMem"],
                           cuda_runtime=cp.cuda.runtime.runtimeGetVersion())
    except Exception:
        info["gpu"] = None
    return info


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def write_csv(path, rows, names=None):
    if names is None:
        names = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def cell_from_record(record):
    return Cell(**{key: record[key] for key in ("strips", "w1_mm", "w2_mm", "gap_mm")})


def channel_flags(config, cell):
    outer_gap = config.period_mm - cell.w1_mm - cell.w2_mm - cell.gap_mm
    half_period = (cell.strips == 2 and abs(cell.w1_mm-cell.w2_mm) < 1e-9
                   and abs(cell.gap_mm-outer_gap) < 1e-9)
    return dict(outer_gap_mm=outer_gap, half_period_symmetry=half_period,
                selected_harmonic_forbidden=half_period and config.harmonic % 2 != 0)


def usable(record):
    return (record["status"] == "accepted" and record.get("alpha_over_k0", -1) >= 0
            and abs(record.get("beta_harmonic_over_k0", 2)) < 1
            and not record.get("selected_harmonic_forbidden", False))


def slice_grid(records, widths, gap):
    size = len(widths)
    beta, alpha = np.full((size, size), np.nan), np.full((size, size), np.nan)
    lookup = {(r["w1_mm"], r["w2_mm"]): r for r in records
              if r["strips"] == 2 and r["gap_mm"] == gap}
    for j, w2 in enumerate(widths):
        for i, w1 in enumerate(widths):
            record = lookup.get((w1, w2))
            if record and usable(record):
                beta[j, i] = record["beta_harmonic_over_k0"]
                alpha[j, i] = record["alpha_over_k0"]
    return beta, alpha


def contour_segments(widths, beta, alpha, level):
    """Piecewise-linear beta contours with co-interpolated alpha.

    Drop a complete grid square if any corner is missing/invalid. No Delaunay
    triangulation over holes, smoothing, extrapolation, or beta-band binning.
    Return endpoint rows (w1, w2, alpha/k0), one pair per segment.
    """
    segments = []
    for j in range(len(widths)-1):
        for i in range(len(widths)-1):
            corners = ((j, i), (j, i+1), (j+1, i+1), (j+1, i))
            b = np.array([beta[y, x] for y, x in corners])
            a = np.array([alpha[y, x] for y, x in corners])
            if not (np.isfinite(b).all() and np.isfinite(a).all()):
                continue
            points = np.array([[widths[x], widths[y], a[k]]
                               for k, (y, x) in enumerate(corners)])
            for triangle in ((0, 1, 2), (0, 2, 3)):
                intersections = []
                for u, v in zip(triangle, triangle[1:] + triangle[:1]):
                    bu, bv = b[u]-level, b[v]-level
                    if abs(bu) < 1e-12:
                        intersections.append(points[u])
                    if bu * bv < 0:
                        t = -bu/(bv-bu)
                        intersections.append(points[u] + t*(points[v]-points[u]))
                unique = []
                for point in intersections:
                    if not any(np.allclose(point, old, rtol=0, atol=1e-10) for old in unique):
                        unique.append(point)
                if len(unique) == 2:
                    segments.append(np.array(unique))
                elif len(unique) == 3:
                    # A flat beta triangle: retain its perimeter and full alpha range.
                    segments.extend(np.array([unique[u], unique[v]]) for u, v in ((0, 1), (1, 2), (2, 0)))
    return segments


def single_crossings(records, single_widths, level):
    lookup = {r["w1_mm"]: r for r in records if r["strips"] == 1}
    points = []
    for w1, w2 in zip(single_widths[:-1], single_widths[1:]):
        first, second = lookup.get(w1), lookup.get(w2)
        if not first or not second or not usable(first) or not usable(second):
            continue
        b1, b2 = first["beta_harmonic_over_k0"], second["beta_harmonic_over_k0"]
        if b1 == b2:
            continue
        t = (level-b1)/(b2-b1)
        if 0 <= t <= 1:
            points.append((w1+t*(w2-w1), first["alpha_over_k0"]
                           + t*(second["alpha_over_k0"]-first["alpha_over_k0"])))
    return points


def merge_intervals(intervals):
    merged = []
    for lower, upper in sorted(intervals):
        if merged and lower <= merged[-1][1] + 1e-12:
            merged[-1][1] = max(upper, merged[-1][1])
        else:
            merged.append([lower, upper])
    return merged


def make_plots(output, config, records, widths, gaps, single_widths, levels):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    accepted = [r for r in records if usable(r) and r["strips"]]
    if not accepted:
        raise ValueError("No accepted radiating-branch points to plot; inspect records.json")
    beta_values = [r["beta_harmonic_over_k0"] for r in accepted]
    if not levels:
        levels = np.linspace(min(beta_values), max(beta_values), 7)[1:-1].tolist()
    all_alpha = [r["alpha_over_k0"] for r in accepted]
    norm = Normalize(0, max(max(all_alpha), 1e-8))
    beta_label = rf"Radiating-harmonic phase $\beta_{{{config.harmonic}}}/k_0$"
    alpha_label = r"Modal attenuation $\alpha/k_0$"
    subtitle = (f"{config.frequency_ghz:g} GHz | period {config.period_mm:g} mm | "
                f"PEC, lossless dielectric | mesh {config.mesh_mm:g} mm (pilot)")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), layout="constrained")
    pairs = [r for r in accepted if r["strips"] == 2]
    if pairs:
        scatter = axes[0].scatter([r["beta_harmonic_over_k0"] for r in pairs],
                                 [r["alpha_over_k0"] for r in pairs],
                                 c=[r["gap_mm"] for r in pairs], cmap="viridis", s=32,
                                 label="Two strips")
        fig.colorbar(scatter, ax=axes[0], label="Internal gap [mm]")
    single_lookup = {r["w1_mm"]: r for r in records if r["strips"] == 1}
    sb = [single_lookup[w]["beta_harmonic_over_k0"] if w in single_lookup and usable(single_lookup[w])
          else np.nan for w in single_widths]
    sa = [single_lookup[w]["alpha_over_k0"] if w in single_lookup and usable(single_lookup[w])
          else np.nan for w in single_widths]
    axes[0].plot(sb, sa, "o-", color="black", ms=4, label="One strip")
    axes[0].set(xlabel=beta_label, ylabel=alpha_label, title="Sampled achievable responses")
    axes[0].legend()
    axes[0].grid(alpha=.2)
    endpoints, ranges = [], []
    colors = plt.get_cmap("viridis")(np.linspace(0, 1, len(gaps)))
    for level in levels:
        intervals = []
        for gap, color in zip(gaps, colors):
            beta, alpha = slice_grid(records, widths, gap)
            for segment_index, segment in enumerate(contour_segments(widths, beta, alpha, level)):
                lo, hi = sorted(segment[:, 2])
                intervals.append((float(lo), float(hi)))
                axes[1].plot([level, level], [lo, hi], color=color, lw=4, solid_capstyle="round")
                for endpoint, (w1, w2, a) in enumerate(segment):
                    endpoints.append(dict(beta_harmonic_over_k0=level, gap_mm=gap,
                                          segment=segment_index, endpoint=endpoint,
                                          w1_mm=w1, w2_mm=w2, alpha_over_k0=a,
                                          alpha_per_m=a*config.k0))
        for width, a in single_crossings(records, single_widths, level):
            axes[1].plot(level, a, "D", color="black", ms=5)
        merged = merge_intervals(intervals)
        if merged:
            ranges.append(dict(beta_harmonic_over_k0=level,
                               angle_deg=float(np.degrees(np.arcsin(level))),
                               alpha_min_per_m=merged[0][0]*config.k0,
                               alpha_max_per_m=merged[-1][1]*config.k0,
                               connected_intervals_over_k0=json.dumps(merged),
                               interval_count=len(merged)))
    for gap, color in zip(gaps, colors):
        axes[1].plot([], [], color=color, lw=4, label=f"Two strips, gap {gap:g} mm")
    axes[1].plot([], [], "D", color="black", ms=5, label="One-strip crossing")
    axes[1].set(xlabel=beta_label, ylabel=alpha_label,
                title="Leakage on constant-phase contours\nLinear interpolation inside valid grid cells")
    axes[1].legend(fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.18), ncols=2)
    axes[1].grid(alpha=.2)
    fig.suptitle(subtitle, fontsize=12)
    fig.savefig(output / "achievable_leakage.png", dpi=180)
    plt.close(fig)

    ncols = min(2, len(gaps))
    nrows = int(np.ceil(len(gaps)/ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 4.6*nrows),
                             squeeze=False, layout="constrained")
    mappable = plt.cm.ScalarMappable(norm=norm, cmap="magma")
    for ax, gap in zip(axes.flat, gaps):
        beta, alpha = slice_grid(records, widths, gap)
        if len(widths) >= 2 and np.isfinite(alpha).any():
            ax.contourf(widths, widths, np.ma.masked_invalid(alpha),
                        levels=np.linspace(norm.vmin, norm.vmax, 24), norm=norm,
                        cmap="magma", corner_mask=False)
            for level in levels:
                segments = contour_segments(widths, beta, alpha, level)
                for segment in segments:
                    ax.plot(segment[:, 0], segment[:, 1], color="cyan", lw=1)
                if segments:
                    middle = segments[len(segments)//2].mean(axis=0)
                    span = widths[-1]-widths[0]
                    position = np.clip(middle[:2], widths[0]+.06*span, widths[-1]-.06*span)
                    ax.text(*position, f"{level:.3f}", fontsize=8, color="black", ha="center", va="center",
                            bbox=dict(facecolor="white", edgecolor="none", alpha=.8, pad=1))
        xx, yy = np.meshgrid(widths, widths)
        valid = np.isfinite(beta)
        ax.scatter(xx[valid], yy[valid], s=6, color="white")
        ax.scatter(xx[~valid], yy[~valid], marker="x", s=28, color="red")
        ax.plot(widths, widths, "--", lw=.8, color="white", alpha=.7)
        ax.set(xlabel="First strip axial width w1 [mm]", ylabel="Second strip axial width w2 [mm]",
               title=f"Internal gap {gap:g} mm")
    for ax in list(axes.flat)[len(gaps):]:
        ax.set_visible(False)
    fig.colorbar(mappable, ax=list(axes.flat)[:len(gaps)], label=alpha_label)
    fig.suptitle(subtitle + "\nCyan: constant radiating-harmonic beta/k0; red crosses: excluded samples", fontsize=11)
    fig.savefig(output / "constant_beta_contours.png", dpi=180)
    plt.close(fig)
    write_csv(output / "contour_segments.csv", endpoints,
              ["beta_harmonic_over_k0", "gap_mm", "segment", "endpoint", "w1_mm", "w2_mm",
               "alpha_over_k0", "alpha_per_m"])
    write_csv(output / "constant_beta_ranges.csv", ranges,
              ["beta_harmonic_over_k0", "angle_deg", "alpha_min_per_m", "alpha_max_per_m",
               "connected_intervals_over_k0", "interval_count"])
    return levels, ranges


def run_sweep(output, config, cells, resume=False, save_modes=False):
    import fem

    record_path = output / "records.json"
    records = json.loads(record_path.read_text()) if resume and record_path.exists() else []
    done = {r["case"] for r in records}
    references = []
    reference_file = output / "reference.h5"
    seed_path = output / "reference.json"
    seed = json.loads(seed_path.read_text()) if resume and seed_path.exists() else None
    if seed and seed["status"] == "accepted" and reference_file.exists():
        modes = fem.load(reference_file)
        mode = modes[seed["selected_mode"]]
        reference = (mode.neff, field_signature(mode, config))
    else:
        seed, reference, modes = solve_cell(config, Cell(0))
        write_json(output / "reference.json", seed)
        modes.save(reference_file)
        if reference is None:
            raise RuntimeError("Could not identify a bound reference mode; inspect reference.json")
    print(f"Reference neff={reference[0]:.7g}; alpha/k0={seed['alpha_over_k0']:.3g}", flush=True)
    # Small signature snapshots make resumption independent of large field archives.
    signatures_dir = output / "signatures"
    signatures_dir.mkdir(exist_ok=True)
    for record in records:
        signature_file = signatures_dir / f"{record['case']}.npy"
        if record["status"] == "accepted" and signature_file.exists():
            references.append((cell_from_record(record),
                               (complex(record["neff_real"], record["neff_imag"]),
                                np.load(signature_file)), record["case"]))
    for index, cell in enumerate(cells, 1):
        if cell.key in done:
            continue
        started = perf_counter()
        same_family = [entry for entry in references if entry[0].strips == cell.strips]
        parent_key, parent = "unperturbed", reference
        if same_family:
            nearest = min(same_family, key=lambda entry:
                          (cell.w1_mm-entry[0].w1_mm)**2 + (cell.w2_mm-entry[0].w2_mm)**2
                          + (cell.gap_mm-entry[0].gap_mm)**2)
            parent, parent_key = nearest[1], nearest[2]
        try:
            cell.intervals(config.period_mm)
        except ValueError as exc:
            record = dict(case=cell.key, **asdict(cell), status="invalid_geometry", error=str(exc))
            tracked = None
        else:
            try:
                record, tracked, modes = solve_cell(config, cell, parent)
                record["parent"] = parent_key
                if tracked is not None:
                    record["seed_overlap"] = overlap(reference[1], tracked[1])
                    if record["seed_overlap"] < 0.6:
                        record.update(status="untracked", error="Low overlap with the original guide branch")
                        tracked = None
                if save_modes:
                    modes_dir = output / "modes"
                    modes_dir.mkdir(exist_ok=True)
                    modes.save(modes_dir / f"{cell.key}.h5")
            except Exception as exc:
                record = dict(case=cell.key, **asdict(cell), status="failed",
                              error=f"{type(exc).__name__}: {exc}", total_seconds=perf_counter()-started)
                tracked = None
        if tracked is not None:
            np.save(signatures_dir / f"{cell.key}.npy", tracked[1])
            references.append((cell, tracked, cell.key))
        records.append(record)
        write_json(record_path, records)
        if record["status"] == "accepted":
            detail = f"alpha/k0={record['alpha_over_k0']:.5g}, beta_h/k0={record['beta_harmonic_over_k0']:.5g}"
        else:
            detail = record.get("error", "")
        print(f"[{index}/{len(cells)}] {cell.key}: {record['status']} {detail}", flush=True)
    return records


def validate_samples(output, config, records):
    """Sensitivity checks at the seed plus one-strip, equal- and unequal-pair cells."""
    seed = json.loads((output / "reference.json").read_text())
    selected = [seed]
    groups = ([r for r in records if usable(r) and r["strips"] == 1],
              [r for r in records if usable(r) and r["strips"] == 2 and r["w1_mm"] == r["w2_mm"]],
              [r for r in records if usable(r) and r["strips"] == 2 and r["w1_mm"] != r["w2_mm"]])
    for group in groups:
        if group:
            selected.append(sorted(group, key=lambda r: r["alpha_over_k0"])[len(group)//2])
    import fem
    seed_modes = fem.load(output / "reference.h5")
    seed_mode = seed_modes[seed["selected_mode"]]
    checks_path = output / "validation.json"
    checks = json.loads(checks_path.read_text()) if checks_path.exists() else []
    variants = (("finer_mesh", replace(config, mesh_mm=config.mesh_mm*.7,
                                       wavelength_elements=int(np.ceil(config.wavelength_elements/0.7)))),
                ("larger_exterior", replace(config, air_mm=config.air_mm*1.5, pml_mm=config.pml_mm*1.5)))
    for original in selected:
        cell = cell_from_record(original)
        signature = (field_signature(seed_mode, config) if cell.strips == 0
                     else np.load(output / "signatures" / f"{cell.key}.npy"))
        reference = (complex(original["neff_real"], original["neff_imag"]), signature)
        for name, variant in variants:
            previous = next((entry for entry in checks
                             if entry["case"] == cell.key and entry["variation"] == name), None)
            if previous and previous["status"] == "accepted":
                continue
            print(f"Validating {cell.key} / {name}...", flush=True)
            try:
                checked, tracked, _ = solve_cell(variant, cell, reference)
                entry = dict(case=cell.key, variation=name, status=checked["status"],
                             mesh_mm=variant.mesh_mm, wavelength_elements=variant.wavelength_elements,
                             air_mm=variant.air_mm, pml_mm=variant.pml_mm)
                if tracked is not None:
                    da = checked["alpha_over_k0"]-original["alpha_over_k0"]
                    db = checked["beta_harmonic_over_k0"]-original["beta_harmonic_over_k0"]
                    entry.update(alpha_over_k0=checked["alpha_over_k0"],
                                 beta_harmonic_over_k0=checked["beta_harmonic_over_k0"],
                                 delta_alpha_over_k0=da, delta_beta_over_k0=db,
                                 relative_alpha_change=abs(da)/max(abs(original["alpha_over_k0"]), 1e-8),
                                 dofs=checked["dofs"], total_seconds=checked["total_seconds"])
            except Exception as exc:
                entry = dict(case=cell.key, variation=name, status="failed", error=str(exc))
            checks = [old for old in checks
                      if (old["case"], old["variation"]) != (cell.key, name)]
            checks.append(entry)
            write_json(checks_path, checks)
            print(f"Validation {cell.key} / {name}: {entry}", flush=True)
    write_csv(output / "validation.csv", checks)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--widths-mm", type=float, nargs="+", default=[.5, 1, 1.5, 2])
    parser.add_argument("--gaps-mm", type=float, nargs="+", default=[.5, 1.5, 2.5, 3.5])
    parser.add_argument("--single-widths-mm", type=float, nargs="+", default=[.25, .5, .75, 1, 1.25, 1.5, 1.75, 2])
    parser.add_argument("--beta-levels", type=float, nargs="*", default=[],
                        help="Constant radiating-harmonic beta/k0 levels (automatic if omitted)")
    for name in ("frequency_ghz", "period_mm", "guide_width_mm", "guide_height_mm", "epsilon_r",
                 "metal_thickness_mm", "air_mm", "pml_mm", "mesh_mm", "neff_guess"):
        parser.add_argument("--" + name.replace("_", "-"), type=float, default=getattr(Config(), name))
    parser.add_argument("--backend", choices=["gpu", "hybrid", "cpu", "auto"], default="gpu")
    parser.add_argument("--wavelength-elements", type=int, default=Config().wavelength_elements)
    parser.add_argument("--harmonic", type=int, default=-1)
    parser.add_argument("--modes", type=int, default=6)
    parser.add_argument("--ncv", type=int, default=96)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--plot-only", action="store_true", help="Use saved configuration and results")
    parser.add_argument("--save-modes", action="store_true")
    parser.add_argument("--validate", action="store_true", help="Run mesh/exterior checks on representative cells")
    parser.add_argument("--max-cases", type=int, help="Stop after this many sweep cases; resume later")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if not args.plot_only or args.validate:
        # A read-only home cache can make Windows tempfile retry indefinitely
        # when CuPy compiles a previously unseen kernel. Keep it in this run.
        cache = output / "cuda_cache"
        cache.mkdir(exist_ok=True)
        os.environ.setdefault("CUPY_CACHE_DIR", str(cache))
    config_path = output / "config.json"
    if args.plot_only:
        settings = json.loads(config_path.read_text())
        config = Config(**settings["solver"])
        widths, gaps, singles = (settings[key] for key in ("widths_mm", "gaps_mm", "single_widths_mm"))
        records = json.loads((output / "records.json").read_text())
    else:
        config = Config(**{key: getattr(args, key) for key in asdict(Config()) if hasattr(args, key)})
        config.validate()
        widths, gaps, singles = (sorted(set(values)) for values in
                                 (args.widths_mm, args.gaps_mm, args.single_widths_mm))
        if any(not np.isfinite(v) or v <= 0 for values in (widths, gaps, singles) for v in values):
            raise ValueError("Sweep dimensions must be finite and positive")
        if len(widths) < 2:
            raise ValueError("At least two pair widths are required for a contour plane")
        settings = dict(schema=1, solver=asdict(config), widths_mm=widths, gaps_mm=gaps,
                        single_widths_mm=singles)
        if config_path.exists():
            if not args.resume:
                raise ValueError("Output already contains a run; choose --resume or another --output")
            if json.loads(config_path.read_text()) != settings:
                raise ValueError("Resume settings differ from config.json; use the original arguments")
        else:
            write_json(config_path, settings)
            write_json(output / "environment.json", environment_info())
        cells = sweep_cells(widths, gaps, singles)
        if args.max_cases is not None:
            if args.max_cases < 1:
                raise ValueError("max-cases must be positive")
            cells = cells[:args.max_cases]
        records = run_sweep(output, config, cells, args.resume, args.save_modes)
    for record in records:
        record.update(channel_flags(config, cell_from_record(record)))
    write_json(output / "records.json", records)
    write_csv(output / "samples.csv", [{key: value for key, value in r.items() if key != "candidates"}
                                      for r in records])
    if any(not np.isfinite(v) or abs(v) >= 1 for v in args.beta_levels):
        raise ValueError("beta-levels must be finite and inside the air light cone (-1, 1)")
    levels, ranges = make_plots(output, config, records, widths, gaps, singles, args.beta_levels)
    if args.validate:
        validate_samples(output, config, records)
    summary = dict(samples=len(records), accepted=sum(r["status"] == "accepted" for r in records),
                   radiating_branch_samples=sum(usable(r) for r in records),
                   beta_levels=levels, ranges=ranges,
                   total_sample_seconds=sum(r.get("total_seconds", 0) for r in records),
                   interpretation="Pilot modal attenuation for PEC and lossless dielectric. "
                   "Light-cone membership does not establish harmonic radiation strength or "
                   "open-stopband suppression. Half-period cells suppress odd space harmonics "
                   "and are excluded for an odd selected harmonic. Contours are interpolated within valid grid squares; "
                   "validation.json contains representative discretization/exterior sensitivity checks.")
    write_json(output / "summary.json", summary)
    print(f"Saved {len(records)} samples ({summary['accepted']} accepted) and plots to {output}", flush=True)


if __name__ == "__main__":
    main()
