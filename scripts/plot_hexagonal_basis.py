"""Plot the lowest E1/E2 TM doublets of a regular hexagonal PEC guide."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from scipy.linalg import polar
from scipy.spatial import cKDTree
from skfem import asm
from skfem.models.poisson import mass

import fem
from fem.common.constants import C_0

OUT = Path(__file__).resolve().parents[1] / 'outputs/hexagonal_basis'


def validate_symmetry(tri, fields):
    """Check closure and C6 characters using independent interior samples."""
    axis = np.linspace(-.8, .8, 35)
    x, y = np.meshgrid(axis, axis)
    keep = x*x+y*y < .8**2
    points = np.column_stack((x[keep], y[keep]))
    theta = np.pi/3
    rotated = points @ np.array([[np.cos(theta), np.sin(theta)],
                                 [-np.sin(theta), np.cos(theta)]])
    result = []
    for start, character in [(0, 1), (2, -1)]:
        interp = [mtri.LinearTriInterpolator(tri, f) for f in fields[start:start+2]]
        original = np.column_stack([np.asarray(f(*points.T)) for f in interp])
        transformed = np.column_stack([np.asarray(f(*rotated.T)) for f in interp])
        representation = np.linalg.lstsq(original, transformed, rcond=None)[0]
        error = np.linalg.norm(original@representation-transformed)/np.linalg.norm(transformed)
        assert error < .005, error
        assert abs(np.trace(representation)-character) < .005
        mirrored = np.column_stack([np.asarray(f(points[:, 0], -points[:, 1])) for f in interp])
        reflection_error = np.linalg.norm(original*np.array([1, -1])-mirrored)/np.linalg.norm(mirrored)
        assert reflection_error < .005, reflection_error
        result.append(dict(C6_character=float(np.trace(representation)),
                           rotation_closure_error=float(error),
                           mirror_parity_error=float(reflection_error)))
    return result


def solve(h, export_transverse=False):
    radius = .01
    angle = np.arange(6)*np.pi/3
    vertices = radius*np.column_stack((np.cos(angle), np.sin(angle)))
    solver = fem.ModeSolver2D(frequency=60e9, x_range=(-.011, .011),
                             y_range=(-.011, .011))
    frame = fem.shapes.Rectangle(bounds=((- .011, .011), (-.011, .011)))
    solver.add_geometry(name='hexagonal PEC wall', material=fem.materials.PEC,
                        shape=fem.shapes.Difference(shape=frame,
                              tool=fem.shapes.Polygon(points=vertices)))
    solver.mesh(max_element_size=h, element_order=2, wavelength_elements=4,
                boundary_refinement=None, interface_refinement=None)
    print(f'Mesh h={h}: {solver.mesh_data.mesh.nelements} triangles', flush=True)
    modes = solver.solve(num_modes=20, neff_guess=1.001,
                         backend='cpu', max_refinements=0)
    tm = []
    for mode in modes:
        kc2 = float(np.real(solver.k0**2-mode.beta**2))
        print(f'{mode.index}: {mode.polarization} kcR={np.sqrt(kc2)*radius:.8f} '
              f'residual={mode.residual:.2e}', flush=True)
        if mode.polarization == 'TM':
            tm.append((kc2, mode))
    tm.sort(key=lambda p: p[0])
    assert len(tm) >= 5
    refined = solver.mesh_data.mesh.refined(1)
    xy = refined.p.T/radius
    triangulation = mtri.Triangulation(*xy.T, refined.t.T)
    records, fields = [], []
    transverse_e, transverse_h = [], []
    for name, selected, kind in [('E1', tm[1:3], 1), ('E2', tm[3:5], 2)]:
        coeffs, phases = [], []
        for _, mode in selected:
            (_, _), (z, basis) = solver.system.basis.split(mode.metadata['coefficient_vector'])
            phase = np.angle(z[np.argmax(abs(z))])
            phases.append(np.exp(-1j*phase))
            coeffs.append(np.real(z*np.exp(-1j*phase)))
            assert mode.residual < 1e-7
        z = np.column_stack(coeffs)
        M = asm(mass, basis)
        # Normalize in the longitudinal L2 inner product, then choose a basis
        # by polar alignment with the two polynomial symmetry templates.
        values, vectors = np.linalg.eigh(z.T @ M @ z)
        whitening = vectors @ np.diag(values**-.5) @ vectors.T
        z = z @ whitening
        x, y = basis.doflocs/(solver.k0*radius)
        templates = np.column_stack((x, y) if kind == 1 else (x*x-y*y, 2*x*y))
        rotation, _ = polar(z.T @ M @ templates)
        z = z @ rotation
        distance, indices = cKDTree(basis.doflocs.T).query(refined.p.T*solver.k0)
        assert distance.max() < 1e-9
        plotted = z[indices]
        scales = np.max(abs(plotted), axis=0)
        plotted /= scales
        fields.extend(plotted.T)
        if export_transverse:
            transform = np.diag(phases) @ whitening @ rotation @ np.diag(1/scales)
            for family, destination in [('E', transverse_e), ('H', transverse_h)]:
                components = []
                for component in [family+'x', family+'y']:
                    raw = np.column_stack([np.ravel(mode.component(component))
                                           for _, mode in selected])
                    aligned = 1j*raw @ transform
                    # Ez is real in the scalar plots. Et and Ht are both -i
                    # times real vectors, so display their common +90deg phase.
                    assert np.linalg.norm(aligned.imag)/np.linalg.norm(aligned.real) < 1e-5
                    components.append(aligned.real)
                destination.extend(np.stack(components, axis=-1).transpose(1, 0, 2))
        k = np.sqrt([p[0] for p in selected])
        split = float(np.ptp(k)/np.mean(k))
        assert split < 1e-3, (name, split)
        records.append(dict(irrep=name, kcR=(k*radius).tolist(),
                            cutoff_GHz=(C_0*k/(2*np.pi*1e9)).tolist(),
                            relative_pair_splitting=split,
                            eigen_residuals=[p[1].residual for p in selected]))
    if export_transverse:
        OUT.mkdir(parents=True, exist_ok=True)
        sampled = selected[0][1].fields
        np.savez_compressed(OUT/'transverse_fields.npz',
                            points=np.column_stack((sampled.x.ravel(), sampled.y.ravel()))/radius,
                            E=np.asarray(transverse_e), H=np.asarray(transverse_h),
                            vertices=vertices/radius)
    return triangulation, np.asarray(fields), vertices/radius, records


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tri, fields, vertices, records = solve(.00055)
    symmetry = validate_symmetry(tri, fields)
    _, _, _, fine = solve(.00038)
    for first, second, check in zip(records, fine, symmetry):
        first.update(check)
        change = abs(np.mean(first['kcR'])/np.mean(second['kcR'])-1)
        first['relative_cutoff_mesh_change'] = float(change)
        assert change < 1e-3
    (OUT/'diagnostics.json').write_text(json.dumps(records, indent=2))
    np.savez_compressed(OUT/'fields.npz', points=np.column_stack((tri.x, tri.y)),
                        triangles=tri.triangles, Ez=fields)
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 9))
    fig.patch.set_facecolor('#f8fafc')
    labels = [r'$E_1$: $x$-like partner', r'$E_1$: $y$-like partner',
              r'$E_2$: $(x^2-y^2)$-like partner', r'$E_2$: $2xy$-like partner']
    closed = np.vstack((vertices, vertices[0]))
    for ax, values, title in zip(axes.flat, fields, labels):
        ax.set_facecolor('#f8fafc')
        image = ax.tripcolor(tri, values, cmap='RdBu_r', vmin=-1, vmax=1,
                             shading='gouraud', rasterized=True)
        ax.tricontour(tri, values, levels=[0], colors='#334155', linewidths=.9,
                      linestyles='dashed')
        ax.plot(*closed.T, color='#0f172a', lw=2)
        ax.set_aspect('equal')
        ax.set_xlim(-1.09, 1.09)
        ax.set_ylim(-1, 1)
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-.8, 0, .8])
        ax.set_xlabel(r'$x/R$')
        ax.set_ylabel(r'$y/R$')
        ax.set_title(title, fontsize=14, pad=10)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle('Hexagonal PEC waveguide: symmetry basis fields',
                 fontsize=19, fontweight='bold', y=.98)
    fig.text(.5, .935, r'Lowest TM doublet in each irrep | Signed $E_z$ | '
             'Quadratic finite elements', ha='center', fontsize=11, color='#475569')
    fig.subplots_adjust(left=.08, right=.85, top=.875, bottom=.13,
                        wspace=.35, hspace=.42)
    cax = fig.add_axes([.9, .28, .022, .43])
    cb = fig.colorbar(image, cax=cax, ticks=[-1, 0, 1])
    cb.set_label(r'$E_z\,/\,\max|E_z|$', labelpad=10)
    fig.text(.08, .045,
             'Black outline: PEC wall ($E_z=0$). Dashed lines: interior nodes.\n'
             'Red / blue: opposite field signs. Each panel is normalized independently; '
             'overall sign is arbitrary.', fontsize=10, color='#475569', linespacing=1.6)
    fig.savefig(OUT/'hexagonal_E1_E2.png', dpi=200, facecolor=fig.get_facecolor())
    fig.savefig(OUT/'hexagonal_E1_E2.svg', facecolor=fig.get_facecolor())
    print(json.dumps(records, indent=2), flush=True)
    print(OUT/'hexagonal_E1_E2.png', flush=True)


if __name__ == '__main__':
    main()
