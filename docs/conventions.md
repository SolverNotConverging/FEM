# Conventions

## Solver families

| Family | Domain | Materials | Guide |
|---|---|---|---|
| Electrostatics | 1D or 2D | Positive scalar, diagonal, or symmetric positive-definite permittivity | [Electrostatics](solvers/electrostatics.rst) |
| Waveguide modes | 1D or 2D cross-section | Scalar or diagonal permittivity/permeability | [Waveguide modes](solvers/waveguide.rst) |
| Periodic modes | 2D or 3D periodic cell | Scalar or diagonal permittivity/permeability | [Periodic modes](solvers/periodic.rst) |
| Scattering | 2D xy (TE/Hz or TM/Ez), 3D vector E | Scalar permittivity/permeability, PEC and SIBC | [Scattering](solvers/scattering.rst) |

Boundary, PML, element-order and backend support depend on the family and
dimension; see the family guide.

## Workflow

1. Construct a solver with keyword arguments and SI quantities.
2. Define reusable `fem.Material` and `fem.shapes` values and assign them to geometry.
3. Configure boundary conditions and, for scattering, the plane-wave direction and polarization.
4. Call `mesh()` to select the initial discretization.
5. Call `solve()` and inspect `solve_info`, including adaptive stopping status for adaptive solvers.
6. Call `plot()`, `show()` or `save(path)`; `fem.load(path)` reopens a result.

Geometry edits invalidate the mesh and result. Remeshing invalidates the result.
Solvers retain explicit mesh settings for automatic remeshing. Frequency is fixed
on a solver; use a supported sweep or construct another solver.

All solver families accept `element_order` and `quadrature_order` in `mesh()`.
The integration degree defaults to 5 for every solver.

| Solver | Default field element order |
| --- | --- |
| Electrostatics, 1D and 2D | 2 |
| Waveguide modes, 1D | 2 |
| Waveguide modes, 2D | 2 |
| Periodic modes, 2D and 3D | 2 |
| Plane-wave scattering, 2D TE and TM | 2 |
| Plane-wave scattering, 3D | 2 |
| Antenna, 3D | 2 |
| Antenna waveguide port modes | Inherit the antenna order, normally 2 |

Every solver supports `mesh(element_order=1)`.
For a 3D antenna or scatterer, this selects a six-function Nédélec tetrahedron
in place of the 20-function second-order element.

Every `solve()` starts on the current mesh. By default it allows two additional
meshes and stops early when the relative discretization estimator reaches
`adaptive_tolerance=0.05`. Set `max_refinements=0` for a fixed mesh.
Inspect `solve_info["adaptive_history"]` for each pass and
`solve_info["adaptive_converged"]` for the final status. The tolerance applies
to the estimator. Accuracy of a particular observable, such as gain or RCS,
can be checked by comparing it across refinement levels.

For a waveguide or periodic dispersion sweep, call
`solver.sweep(frequencies=frequencies, num_modes=...)` and iterate over the
mode sets in increasing frequency order. The solver adapts at the highest
frequency first and reassembles the operators on that final mesh for each
lower frequency. Antenna sweeps use the same refinement policy and return
an `AntennaSweep` with a frequency selector in the viewer.
The `tqdm` bar counts completed frequencies and displays the active frequency
in Hz. Each adaptive pass prints the element count, relative residual,
tolerance, and decision to refine or stop.

## Units and signs

Lengths are in metres, frequencies in hertz and constitutive values are relative.
Fields use `exp(+i*omega*t)` and guided propagation `exp(-i*beta*z)`. Passive
materials have non-positive imaginary parts and forward attenuation is
`alpha = -Im(beta)`. Periodic results store field envelopes with the Bloch phase
recorded separately. Propagating periodic modes have unit longitudinal power,
averaged over one period of the envelope and integrated over the physical
cross-section. The units are W/m of invariant width in 2D and W in 3D.

## Vocabulary

The same names mean the same thing in every family:

| Name | Meaning |
|---|---|
| `frequency`, `frequencies` | Frequency in hertz (one value, or a sweep) |
| `neff` | Complex effective index |
| `beta`, `alpha` | Propagation constant (rad/m) and attenuation `-Im(beta)` (1/m) |
| `num_modes`, `neff_guess` | Number of modes requested and the shift-invert target |
| `max_element_size`, `element_order` | Mesh size in metres and finite-element order |
| `max_refinements`, `adaptive_tolerance` | Adaptive budget; `0` performs a fixed-mesh solve |
| `backend` | `"auto"`, `"cpu"` or `"gpu"` |
| `solve_info`, `mesh_data` | Solver diagnostics and the mesh snapshot of a result |

Python indices (modes, cases) are zero-based. Scalar controls accept Python and
NumPy numbers, including integral floats for counts; booleans, numeric strings
and arrays are rejected. Errors derive from `fem.FEMError`.
