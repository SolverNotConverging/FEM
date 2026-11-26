"""Transverse electric and magnetic arrows for the aligned hexagonal TM pairs."""
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as PolygonPath
from scipy.interpolate import LinearNDInterpolator

from plot_hexagonal_basis import OUT, solve


def plot():
    data = np.load(OUT/'transverse_fields.npz')
    vertices = data['vertices']
    closed = np.vstack((vertices, vertices[0]))
    axis = np.linspace(-1, 1, 19)
    gx, gy = np.meshgrid(axis, axis)
    points = np.column_stack((gx.ravel(), gy.ravel()))
    inside = PolygonPath(vertices).contains_points(points, radius=-.075)
    points = points[inside]
    dense = np.linspace(-1.01, 1.01, 350)
    xx, yy = np.meshgrid(dense, dense)
    mask = PolygonPath(vertices).contains_points(np.column_stack((xx.ravel(), yy.ravel())))
    labels = [r'$E_1$: $x$-like partner', r'$E_1$: $y$-like partner',
              r'$E_2$: $(x^2-y^2)$-like partner', r'$E_2$: $2xy$-like partner']
    for family, word in [('E', 'Electric'), ('H', 'Magnetic')]:
        fig, axes = plt.subplots(2, 2, figsize=(10.2, 9))
        fig.patch.set_facecolor('#f8fafc')
        for ax, field, title in zip(axes.flat, data[family], labels):
            interp = LinearNDInterpolator(data['points'], field)
            dense_vectors = interp(xx, yy)
            magnitude = np.linalg.norm(dense_vectors, axis=-1)
            scale = np.nanmax(magnitude)
            magnitude = np.ma.masked_where(~mask.reshape(xx.shape) | ~np.isfinite(magnitude), magnitude/scale)
            image = ax.pcolormesh(xx, yy, magnitude, vmin=0, vmax=1,
                                  cmap='Blues', alpha=.6, shading='auto', rasterized=True)
            vectors = interp(points)/scale
            assert np.isfinite(vectors).all()
            ax.quiver(*points.T, *vectors.T, color='#102a43', angles='xy',
                      scale_units='xy', scale=9, width=.005, pivot='mid',
                      headwidth=3.8, headlength=4.5, headaxislength=4)
            ax.plot(*closed.T, color='#0f172a', lw=1.8)
            ax.set_aspect('equal')
            ax.set_xlim(-1.09, 1.09)
            ax.set_ylim(-1, 1)
            ax.set_xticks([-1, 0, 1])
            ax.set_yticks([-.8, 0, .8])
            ax.set_xlabel(r'$x/R$')
            ax.set_ylabel(r'$y/R$')
            ax.set_title(title, fontsize=14, pad=9)
            ax.set_facecolor('#f8fafc')
            for spine in ax.spines.values():
                spine.set_visible(False)
        fig.suptitle(f'{word} transverse field in a hexagonal PEC guide',
                     fontsize=18, fontweight='bold', y=.98)
        fig.text(.5, .935, rf'$({family}_x, {family}_y)$ arrows | '
                 'Lowest TM doublets | Vacuum | 60 GHz',
                 ha='center', fontsize=11, color='#475569')
        fig.subplots_adjust(left=.08, right=.85, top=.875, bottom=.13,
                            wspace=.35, hspace=.42)
        cax = fig.add_axes([.9, .28, .022, .43])
        cb = fig.colorbar(image, cax=cax, ticks=[0, .5, 1])
        cb.set_label(r'$|\mathbf{'+family+r'}_t|\,/\,\max|\mathbf{'+family+r'}_t|$', labelpad=10)
        fig.text(.08, .045,
                 'Arrow length and colour indicate transverse magnitude; each panel is normalized.\n'
                 r'Common phase: arrows show $\mathrm{Re}(i\mathbf{'+family+r'}_t)$ '
                 r'when the previous $E_z$ plots are real. $R=10$ mm.',
                 fontsize=10, color='#475569', linespacing=1.6)
        stem = OUT/f'hexagonal_{family}t_arrows'
        fig.savefig(stem.with_suffix('.png'), dpi=190, facecolor=fig.get_facecolor())
        fig.savefig(stem.with_suffix('.svg'), facecolor=fig.get_facecolor())
        plt.close(fig)
        print(stem.with_suffix('.png'), flush=True)
    # In vacuum TM modes, Ht is parallel to z cross Et, and forward power is positive.
    e, h = data['E'], data['H']
    dot = np.sum(e*h, axis=-1)
    cosine = np.abs(dot)/np.maximum(np.linalg.norm(e, axis=-1)*np.linalg.norm(h, axis=-1), 1e-30)
    power = e[..., 0]*h[..., 1]-e[..., 1]*h[..., 0]
    checks = dict(max_abs_E_H_cosine=float(cosine.max()),
                  minimum_forward_poynting=float(power.min()))
    assert checks['max_abs_E_H_cosine'] < 1e-3
    assert checks['minimum_forward_poynting'] > -1e-7*power.max()
    (OUT/'transverse_checks.json').write_text(json.dumps(checks, indent=2))
    print(checks)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--solve', action='store_true')
    args = parser.parse_args()
    if args.solve or not (OUT/'transverse_fields.npz').exists():
        solve(.00055, export_transverse=True)
    plot()
