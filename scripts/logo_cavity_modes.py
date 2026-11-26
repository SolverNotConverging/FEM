"""Solve the gprMax letter cavities with this repository's FEM ModeSolver2D."""
from pathlib import Path
import argparse
import json
import time

import matplotlib
matplotlib.use('Agg')
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.figure import Figure
from PIL import Image, ImageFont
import numpy as np
from scipy.spatial import cKDTree
import fem
from fem import shapes, materials
from fem.common.constants import C_0

ROOT = Path(__file__).resolve().parents[1]
FONT = ROOT.parent / 'gprMin/gprMax/toolboxes/gprMaxLogo/fonts/IBMPlexSans-Bold.ttf'
OUTPUT = ROOT / 'outputs/logo_cavity_modes'
EM_SIZE = 2988 * (0.00025 / 3)


def contours(letter):
    if letter == 'a':
        return raster_contours(letter)
    paths = TextPath((0, 0), letter, size=1000,
                     prop=FontProperties(fname=str(FONT))).to_polygons()
    return [p[:-1] * (EM_SIZE / 1000) for p in paths]


def simplify(points, tolerance):
    """Douglas-Peucker simplification, retaining subpixel geometry detail."""
    if len(points) <= 2:
        return points
    delta = points[-1]-points[0]
    if np.dot(delta, delta) == 0:
        distance = np.linalg.norm(points-points[0], axis=1)
    else:
        t = np.clip((points-points[0])@delta/np.dot(delta, delta), 0, 1)
        distance = np.linalg.norm(points-points[0]-t[:, None]*delta, axis=1)
    index = int(np.argmax(distance))
    if distance[index] <= tolerance:
        return points[[0, -1]]
    return np.vstack((simplify(points[:index+1], tolerance)[:-1],
                      simplify(points[index:], tolerance)))


def raster_contours(letter):
    """Use the authoritative mask where flattened font paths self-intersect."""
    toolbox = FONT.parent.parent
    meta = json.loads((toolbox/'model/gprmax_v4_logo.json').read_text())
    glyph = next(g for g in meta['glyphs'] if g['character'] == letter)
    x0, y0, x1, y1 = glyph['bounds_indices']
    with Image.open(toolbox/'model/gprmax_v4_logo_geometry.png') as image:
        mask = np.asarray(image.crop((x0-1, y0-1, x1+1, y1+1))) > 128
    font = ImageFont.truetype(str(FONT), 2988)
    advances = [font.getlength(c) for c in 'gprMax']
    boxes = [font.getbbox(c, anchor='ls') for c in 'gprMax']
    top, bottom = min(b[1] for b in boxes), max(b[3] for b in boxes)
    baseline = (6000-(bottom-top))/2-top
    origin_x = (12000-sum(advances))/2+sum(advances[:'gprMax'.index(letter)])
    dl = EM_SIZE/2988
    x = (np.arange(mask.shape[1])+x0-1-origin_x)*dl
    y = (baseline-(np.arange(mask.shape[0])+y0-1))*dl
    figure = Figure()
    axis = figure.subplots()
    paths = axis.contour(x, y, mask.astype(float), levels=[.5]).allsegs[0]
    # A 0.06 mm tolerance is smaller than one production FDTD pixel.
    return [simplify(p, .00006)[:-1] for p in paths if len(p) > 4]


def make_solver(letter, h):
    paths = contours(letter)
    areas = [abs(np.sum(p[:, 0] * np.roll(p[:, 1], -1)
                        - p[:, 1] * np.roll(p[:, 0], -1))) for p in paths]
    outer = int(np.argmax(areas))
    p = paths[outer]
    lo, hi = p.min(axis=0) - .003, p.max(axis=0) + .003
    solver = fem.ModeSolver2D(frequency=15e9, x_range=(lo[0], hi[0]),
                             y_range=(lo[1], hi[1]), boundary=materials.PEC)
    frame = shapes.Rectangle(bounds=((lo[0], hi[0]), (lo[1], hi[1])))
    solver.add_geometry(shape=shapes.Difference(shape=frame, tool=shapes.Polygon(points=p)),
                        material=materials.PEC, name='surrounding PEC')
    for i, hole in enumerate(paths):
        if i != outer:
            solver.add_geometry(shape=shapes.Polygon(points=hole), material=materials.PEC,
                                name=f'PEC counter {i}')
    solver.mesh(max_element_size=h, element_order=2, wavelength_elements=4,
                boundary_refinement=None, interface_refinement=None)
    return solver, paths


def solve_letter(letter, h, count):
    start = time.perf_counter()
    solver, paths = make_solver(letter, h)
    print(f'{letter}: meshed {solver.mesh_data.elements.shape}, solving {count} vector modes', flush=True)
    modes = solver.solve(num_modes=count, neff_guess=1.001, backend='cpu', max_refinements=0)
    selected = []
    for i, mode in enumerate(modes):
        fc = float(np.real(C_0 / (2*np.pi) * np.sqrt(solver.k0**2 - mode.beta**2 + 0j)))
        print(f'{letter} {i}: {mode.polarization} fc={fc/1e9:.6f}GHz '
              f'residual={mode.residual:.2e}', flush=True)
        if mode.polarization == 'TM':
            selected.append((fc, mode))
    selected.sort(key=lambda item: item[0])
    OUTPUT.mkdir(parents=True, exist_ok=True)
    modes.save(OUTPUT / f'{letter}_h{h:g}_all_modes.h5')
    refined = solver.mesh_data.mesh.refined(1)
    arrays = {'points': refined.p.T, 'triangles': refined.t.T}
    records = []
    for rank, (fc, mode) in enumerate(selected):
        (_, _), (z, basis) = solver.system.basis.split(mode.metadata['coefficient_vector'])
        # One uniform refinement adds precisely the P2 edge-midpoint DOFs.
        # Match those directly to avoid point-finder roundoff on PEC edges.
        distance, indices = cKDTree(basis.doflocs.T).query(refined.p.T * solver.k0)
        assert distance.max() < 1e-10
        values = z[indices]
        phase = np.angle(values[np.argmax(abs(values))])
        values = np.real(values * np.exp(-1j*phase))
        arrays[f'Ez_{rank}'] = values / np.max(abs(values))
        records.append(dict(rank=rank+1, frequency_hz=fc, residual=mode.residual,
                            divergence_residual=mode.divergence_residual))
    np.savez_compressed(OUTPUT / f'{letter}_h{h:g}.npz', **arrays)
    record = dict(letter=letter, h=h, element_order=2, frequency_hz=15e9,
                  solver='fem.ModeSolver2D', modes=records,
                  seconds=time.perf_counter()-start, contours=[p.tolist() for p in paths])
    (OUTPUT / f'{letter}_h{h:g}.json').write_text(json.dumps(record, indent=2))
    print(f'{letter}: {len(selected)} TM modes, {record["seconds"]:.1f}s', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--letters', default='gprMax')
    parser.add_argument('--h', type=float, default=.004)
    parser.add_argument('--count', type=int, default=36)
    args = parser.parse_args()
    for letter in args.letters:
        solve_letter(letter, args.h, args.count)
