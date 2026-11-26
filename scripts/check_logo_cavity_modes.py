"""Check the TM cutoff reduction of the FEM pencil and mesh convergence.

At beta=0 in vacuum the longitudinal block decouples: A0_zz=K-M.
The existing FEM Gauss matrix is -i*M. Thus K*z=lambda*M*z gives
f_c=f_reference*sqrt(lambda). All PEC constraints come from ModeSolver2D.
"""
import json
import numpy as np
from fem.common.linalg import generalized_eigenpairs
from logo_cavity_modes import make_solver, OUTPUT


def cutoff(letter, h):
    solver, _ = make_solver(letter, h)
    system = solver.system
    # Map free longitudinal indices into the reduced vector and scalar basis.
    reduced_z = np.flatnonzero(np.isin(system.free_dofs, system.longitudinal_indices))
    scalar_z = np.searchsorted(system.longitudinal_indices, system.free_dofs[reduced_z])
    mass = (1j*system.gauss_longitudinal)[scalar_z][:, scalar_z].tocsr()
    stiffness = system.A0[reduced_z][:, reduced_z] + mass
    values, vectors, method = generalized_eigenpairs(stiffness, mass, target=0., count=5,
                                                   backend='cpu', tolerance=1e-11)
    order = np.argsort(values.real)
    values, vectors = values[order], vectors[:, order]
    residuals = []
    for value, v in zip(values, vectors.T):
        left, right = stiffness@v, value*(mass@v)
        residuals.append(float(np.linalg.norm(left-right)/(np.linalg.norm(left)+np.linalg.norm(right))))
    return dict(h=h, frequencies_hz=(solver.frequency*np.sqrt(values.real)).tolist(),
                residuals=residuals, elements=len(solver.mesh_data.elements), method=method)


if __name__ == '__main__':
    records = {}
    for letter in 'gprMax':
        coarse = cutoff(letter, .004)
        fine = cutoff(letter, .0025)
        difference = np.abs(np.array(coarse['frequencies_hz'])/fine['frequencies_hz']-1)
        records[letter] = dict(coarse=coarse, fine=fine, relative_mesh_change=difference.tolist())
        print(f'{letter}: max mesh change {max(difference):.3%}; '
              f'fine residual {max(fine["residuals"]):.2e}', flush=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT/'mesh_convergence.json').write_text(json.dumps(records, indent=2))
    # The full-vector run must return the actual lowest five TM modes.
    for letter, record in records.items():
        path = OUTPUT/f'{letter}_h0.004.json'
        if path.exists():
            vector = json.loads(path.read_text())['modes'][:5]
            assert len(vector) == 5, f'Missing low TM modes for {letter}'
            error = np.max(np.abs(np.array([v['frequency_hz'] for v in vector]) /
                                  record['coarse']['frequencies_hz'] - 1))
            assert error < 1e-6, (letter, error)
            print(f'{letter}: full-vector / cutoff agreement {error:.2e}', flush=True)
