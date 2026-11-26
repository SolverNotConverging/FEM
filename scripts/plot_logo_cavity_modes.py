"""Assemble per-letter FEM cavity fields into five complete wordmarks."""
from pathlib import Path
import argparse
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.font_manager import FontProperties
from PIL import ImageFont
import numpy as np

from logo_cavity_modes import OUTPUT, FONT, EM_SIZE


def plot(h=.004, count=5, magnitude=False):
    font = ImageFont.truetype(str(FONT), 2988)
    advance = np.array([font.getlength(c) * EM_SIZE / 2988 for c in 'gprMax'])
    offsets = np.r_[0, advance.cumsum()[:-1]]
    fig, axes = plt.subplots(count, 1, figsize=(10, 2.75 * count + 1.2))
    fig.patch.set_facecolor('#f9fafc')
    cmap = 'magma' if magnitude else 'RdBu_r'
    norm = Normalize(0 if magnitude else -1, 1)
    for rank, ax in enumerate(axes):
        ax.set_facecolor('#f9fafc')
        for letter, offset, width in zip('gprMax', offsets, advance):
            base = OUTPUT / f'{letter}_h{h:g}'
            record = json.loads(Path(str(base)+'.json').read_text())
            with np.load(Path(str(base)+'.npz')) as data:
                xy = data['points'].copy()
                xy[:, 0] += offset
                tri = mtri.Triangulation(xy[:, 0], xy[:, 1], data['triangles'])
                z = data[f'Ez_{rank}']
                ax.tripcolor(tri, abs(z) if magnitude else z, cmap=cmap, norm=norm,
                             shading='gouraud', rasterized=True)
            for p in record['contours']:
                p = np.array(p)
                p = np.vstack((p, p[0]))
                ax.plot(p[:, 0]+offset, p[:, 1], color='#313b4d', lw=.45)
            fc = record['modes'][rank]['frequency_hz']/1e9
            ax.text(offset + width/2, -.077, f'{fc:.3f} GHz', ha='center', va='center',
                    fontsize=11, color='#475569')
        label = 'Fundamental TM mode' if rank == 0 else f'Higher TM mode {rank}'
        ax.text(.005, .218, f'{rank+1:02d}   {label}', fontsize=12, fontweight='bold', color='#172033')
        ax.set_xlim(-.005, advance.sum()+.005)
        ax.set_ylim(-.098, .244)
        ax.set_aspect('equal')
        ax.axis('off')
    fig.suptitle('gprMax  |  Letter-cavity eigenmodes', x=.07, y=.985, ha='left',
                 fontsize=21, fontweight='bold', color='#172033')
    fig.text(.07, .952, 'FEM ModeSolver2D  /  PEC walls  /  vacuum  /  IBM Plex Sans Bold',
             fontsize=11.5, color='#475569')
    fig.text(.07, .027, 'Each letter is solved independently. Rows show increasing TM cavity frequency within each letter.\n'
             'Fields are normalized per letter and mode; frequencies use the original logo scale (em = 0.249 m).',
             fontsize=9, color='#475569', linespacing=1.6)
    fig.subplots_adjust(left=.06, right=.88, top=.927, bottom=.071, hspace=.12)
    cax = fig.add_axes([.91, .35, .015, .31])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(r'$|E_z|\,/\,\max|E_z|$' if magnitude else r'$E_z\,/\,\max|E_z|$', fontsize=11)
    cb.set_ticks([0, .5, 1] if magnitude else [-1, 0, 1])
    stem = OUTPUT / ('logo_tm_modes_magnitude' if magnitude else 'logo_tm_modes')
    fig.savefig(stem.with_suffix('.png'), dpi=210, facecolor=fig.get_facecolor())
    fig.savefig(stem.with_suffix('.pdf'), dpi=220, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(stem.with_suffix('.png'))


def validate_and_export(h=.004, count=5):
    convergence = json.loads((OUTPUT/'mesh_convergence.json').read_text())
    rows = []
    for letter in 'gprMax':
        record = json.loads((OUTPUT/f'{letter}_h{h:g}.json').read_text())
        selected = record['modes'][:count]
        assert len(selected) == count
        scalar = convergence[letter]['coarse']['frequencies_hz'][:count]
        agreement = np.max(np.abs(np.array([m['frequency_hz'] for m in selected])/scalar-1))
        assert agreement < 1e-6, (letter, agreement)
        with np.load(OUTPUT/f'{letter}_h{h:g}.npz') as data:
            assert np.min(data['Ez_0']) > -1e-5, f'Fundamental sign changes: {letter}'
            for rank, mode in enumerate(selected):
                assert mode['residual'] < 1e-8
                assert np.isfinite(data[f'Ez_{rank}']).all()
                rows.append(dict(letter=letter, mode=rank+1,
                                 frequency_GHz=mode['frequency_hz']/1e9,
                                 eigen_residual=mode['residual'],
                                 mesh_change_percent=100*convergence[letter]['relative_mesh_change'][rank]))
        print(f'{letter}: lowest-five TM check passed; frequency agreement {agreement:.2e}')
    with (OUTPUT/'mode_frequencies.csv').open('w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (OUTPUT/'README.txt').write_text(
        'gprMax letter cavity modes, computed with this FEM repository\n\n'
        'Method: fem.ModeSolver2D; quadratic Nedelec/P2 elements; vacuum; PEC walls.\n'
        'The full-vector propagation solve is evaluated at 15 GHz. For a uniform\n'
        'vacuum guide, fc = c/(2*pi)*sqrt(k0^2-beta^2) gives the TMz cavity\n'
        'resonance at beta=0. The first five TM modes of each letter are retained.\n'
        'Their ordering is independently checked using the longitudinal cutoff\n'
        'block of the same FEM pencil, K Ez = lambda M Ez.\n\n'
        'Geometry: original IBM Plex Sans Bold font, em=0.249 m, no tracking.\n'
        'Smooth flattened outlines are used for g,p,r,M,x. The flattened a\n'
        'outline self-intersects, so a uses contours of the authoritative\n'
        'production FDTD glyph mask, simplified to 0.06 mm (below one pixel).\n'
        'Letter holes remain PEC. Each letter is an independent cavity.\n\n'
        'Mesh: 4 mm maximum target, second-order elements; verified against\n'
        '2.5 mm target. Largest change in the 30 frequencies: '+
        f"{max(row['mesh_change_percent'] for row in rows):.5f} percent.\n"
        'This measures mesh sensitivity, not an absolute frequency-error bound.\n\n'
        'Plots: signed Ez normalized independently for every letter/mode;\n'
        'the magnitude variant uses the same fields. Each row combines modal\n'
        'rank across letters; it is not a common-frequency coupled eigenmode.\n'
        'Overall eigenvector sign is arbitrary. Frequencies scale inversely\n'
        'with uniform changes in logo dimensions.\n\n'
        'Reproduce from repository root with .venv/Scripts/python.exe:\n'
        '  scripts/logo_cavity_modes.py --h .004 --count 36\n'
        '  scripts/check_logo_cavity_modes.py\n'
        '  scripts/plot_logo_cavity_modes.py\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--h', type=float, default=.004)
    parser.add_argument('--count', type=int, default=5)
    args = parser.parse_args()
    validate_and_export(args.h, args.count)
    plot(args.h, args.count)
    plot(args.h, args.count, magnitude=True)
