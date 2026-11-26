# Getting started

FEM is a pure-Python package (`py3-none-any`). Every dependency, including Gmsh,
HDF5, PyVista and Qt for Python (PySide6), installs from PyPI wheels; no compiler,
CMake or system libraries are needed. Python 3.11–3.13 on Windows, Linux and macOS.

## Install

```sh
git clone https://github.com/SolverNotConverging/Electromagnetics-FEM.git
cd Electromagnetics-FEM
python -m pip install -e .
```

For development with the locked toolchain, run `uv sync --locked` instead.

Optional extras: `pip install -e .[cuda12]` (or `cuda13`) for GPU factorization with
NVIDIA cuDSS. See [acceleration](acceleration.md).

## First solve

```python
import fem

solver = fem.ModeSolver2D(frequency=10e9, x_range=(0, 22.86e-3), y_range=(0, 10.16e-3))
modes = solver.solve(num_modes=2, max_refinements=2)
print(modes.neff)
modes.save("modes.h5")
modes.show()
```

Construct the solver with SI keyword arguments, add materials and geometry,
then call `solve()`. It builds the initial mesh automatically. Use `plot()`,
`show()` or `save(path)` to inspect or save the result. `fem.load(path)` reopens any saved result.

## Run examples

Examples live in `examples/<family>/`. Each is a short top-to-bottom script that
uses only the public `fem` API; change the numbers at the top and run it again:

```sh
python examples/waveguide/1d_parallel_plate_waveguide.py
python examples/periodic/2d_uniform_cell.py
```

Each example writes to an ignored `outputs/<example>/` folder next to the script.

## View results

`result.show()` opens Matplotlib for 2D plane-wave scattering and VTK for 3D scattering. Both display the incident, scattered and total field and RCS. In the modal and electrostatic viewers, choose the
case (frequency), mode, field component and quantity, and toggle:

- **Mesh** — the finite-element mesh drawn over the field.
- **Material** — the material regions, labelled with their relative permittivity.

Components cover E and H (Cartesian components and magnitudes); electrostatic
results show the potential, E and D. In 3D periodic cells the vector fields `E`
and `H` are drawn as arrows coloured by magnitude (`real`/`imag` select the
phasor part). Saved archives open from
the command line:

```sh
fem view examples/periodic/outputs/2d_uniform_cell/results.h5
fem inspect result.h5     # print a short summary without opening a window
fem info                  # version and automatic solver backend
```

Headless machines can set `QT_QPA_PLATFORM=offscreen`; `result.plot()` returns a
Matplotlib figure without opening a window.
