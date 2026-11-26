Dielectric image-line perturbation pilot
============================================

``scripts/image_line_pilot.py`` sweeps one or two PEC strips on a grounded
dielectric image guide at one fixed frequency and period. It extracts the
complex Bloch propagation constant and plots the available attenuation along
constant-phase contours. It runs without opening a viewer.

Geometry and default experiment
--------------------------------------------

The guide runs along z. Its bottom, y=0, is the PEC ground plane; the outer
x faces and upper y face have PML. The ground has no PML. The dielectric is
lossless with relative permittivity 6. The default guide is 6 mm wide in x and
2 mm high, with a 10 mm period and an operating frequency of 20 GHz.

Each strip spans the full 6 mm guide width and has a 0.1 mm PEC thickness.
The swept widths ``w1`` and ``w2`` are the strip dimensions along z.
``gap`` is the edge-to-edge distance between two strips, not their centre
spacing. The pair is centred in its cell, leaving an outer periodic gap of
``period - w1 - w2 - gap``. Nonpositive outer gaps are recorded as invalid
geometries. A single strip is centred and has no internal-gap parameter.

The default sweep contains eight one-strip widths from 0.25 to 2 mm and
64 two-strip cells: independent widths of 0.5, 1, 1.5 and 2 mm, at internal
gaps of 0.5, 1.5, 2.5 and 3.5 mm. An additional unperturbed image guide
anchors the tracked branch. The two widths are independent; the equal-width
subset remains available for comparison.

The default maximum mesh size is 1 mm with 12 elements per wavelength.
Material-aware meshing also refines the dielectric. A coarser 1.5 mm / 8
elements-per-wavelength pilot showed changes of approximately 20--27 percent
in alpha for some cells under refinement, motivating this finer default.

Run and resume
--------------

From the repository root, in an installed FEM environment:

.. code-block:: sh

   python scripts/image_line_pilot.py --validate

The default ``--backend gpu`` explicitly requests the resident CUDA solver.
Use ``--backend hybrid`` or ``--backend cpu`` to select another implementation.
``--backend auto`` follows the package's automatic hybrid/CPU selection.

Change geometry, discretization, sampling or contour levels explicitly:

.. code-block:: sh

   python scripts/image_line_pilot.py --output outputs/image_line_dense \
       --widths-mm 0.5 0.75 1 1.25 1.5 1.75 2 \
       --gaps-mm 0.5 1.5 2.5 3.5 --mesh-mm 1 \
       --wavelength-elements 12 --beta-levels 0.12 0.16 0.20 0.24 0.28 0.32

Lengths in command-line options are in millimetres, frequency is in GHz,
and the saved physical alpha/beta quantities are in inverse metres.
Run ``--help`` for the other geometry and eigensolver controls.

Results default to ``outputs/image_line_pilot``. Each completed case is
checkpointed. Resident GPU arrays and unused CuPy pool blocks are released
between geometries so a long sweep does not accumulate device memory.
Unless ``CUPY_CACHE_DIR`` is already set, compiled CUDA kernels are cached
under the run's ``cuda_cache`` directory, which is writable with the results.
``--max-cases 3`` runs a small initial subset. Resume with the
same simulation/sampling options and ``--resume``; omit ``--max-cases`` to
finish the full sweep. Resume skips all recorded cases, including recorded
failures. Use a new output directory when changing the mesh or solver
configuration, or when retrying failed cases with new settings.

Replot an existing run without solving:

.. code-block:: sh

   python scripts/image_line_pilot.py --output outputs/image_line_pilot \
       --plot-only --beta-levels 0.12 0.16 0.20 0.24 0.28 0.32

The simulation settings and sample axes are loaded from the saved configuration
in plot-only mode. Add ``--validate`` to run the sensitivity checks against
that saved configuration. ``--save-modes`` retains the candidate-mode HDF5
archive for every solved sample; the unperturbed reference is always saved.

Mode selection and interpretation
--------------------------------------------

The seed is a forward, confined mode of the unperturbed guide with real
effective index between 1 and sqrt(epsilon_r). Subsequent cells use the
nearest accepted geometry of the same family as their reference. Selection
requires forward power, passive attenuation, sufficient dielectric energy,
limited PML energy, continuity of complex effective index, and field overlap.
Ambiguous and untracked cases are retained with their candidate diagnostics
and excluded from the contour interpolation.

Overlap uses electric and magnetic periodic envelopes, sampled by nearest
element barycentre on a fixed grid inside the dielectric. It is independent
of the mesh DOF numbering and arbitrary global modal phase. It is an
approximate branch-tracking diagnostic, not a biorthogonal modal projection.
The original guide signature also guards against cumulative branch drift.

Only for the nonradiating reference, small signed alpha/k0 errors below
0.001 in magnitude are allowed during branch identification. The raw signed
value is retained as a numerical diagnostic, and no attenuation correction
is applied to the strip samples. Growing strip modes remain excluded.

The code convention is exp(+j omega t - j beta z), so alpha = -Im(beta).
The CSV contains both Re(beta) and the selected space-harmonic phase:

.. math::

   \beta_n = \operatorname{Re}\beta + 2\pi n/d.

The default harmonic is n=-1. Contour labels and ``--beta-levels`` are
beta_n/k0, not the raw effective index. The corresponding phase-matching
angle from broadside is asin(beta_n/k0), provided its magnitude is below 1.
At fixed period, constant beta_n and constant Re(beta) are equivalent.
An unperturbed guide can have a folded harmonic inside the light cone with
zero amplitude: light-cone membership alone does not prove radiation.

Two identical strips with equal internal and outer gaps have an actual
period of d/2. Odd harmonics of the nominal d-period cell are then absent
on the tracked branch. The script flags these cells, keeps their raw modal
attenuation, and excludes them from odd-harmonic leakage contours. For the
default grid this occurs at w1=w2=1.5 mm and gap=3.5 mm. Its small numerical
attenuation must not be interpreted as radiation into the n=-1 channel.

With PEC strips and lossless dielectric, attenuation of a converged physical
leaky mode is associated with radiation. These pilot values remain modal
attenuation estimates until mesh/exterior sensitivity and the physical
radiation channel are checked. No numerical seed attenuation is subtracted.
The script does not infer harmonic power fractions, far fields, finite-array
matching, or open-stopband suppression from alpha alone. The chosen default
branch is off broadside; demonstrating stopband suppression requires a
separate dispersion or frequency study across broadside.

Plots and data
--------------

* ``achievable_leakage.png`` compares sampled alpha-beta responses and
  interpolated leakage intervals at specified phase constants. Black diamonds
  mark intersections with the one-strip curve. A one-parameter strip family
  provides curve intersections, rather than a two-dimensional contour.
* ``constant_beta_contours.png`` shows alpha in each fixed-gap width-width
  plane, with constant-beta contours in cyan. The equal-width line is dashed.
  Red crosses mark missing or excluded samples.
* ``samples.csv`` contains selected-mode properties, tracking diagnostics,
  geometry, DOFs and timings. ``records.json`` also contains all returned
  candidate modes and failure explanations.
* ``contour_segments.csv`` records both endpoints of every interpolated
  constant-beta segment, including geometry and alpha. ``constant_beta_ranges.csv``
  retains disconnected attainable intervals as well as their overall extrema.
* ``config.json``, ``environment.json``, ``reference.json``, ``reference.h5``
  and the small ``signatures`` directory make the run inspectable and resumable.
* ``validation.csv/json`` records the optional sensitivity checks.

Contours use piecewise-linear interpolation on two triangles per sampled
grid square. A square with any invalid corner is omitted entirely. There is
no extrapolation or interpolation across failed points. Contour endpoints
are interpolated candidate geometries, not new full-wave evaluations.
Do not treat the extrema as proof that every intermediate leakage is
available; the saved connected intervals preserve gaps in the sampled set.

Validation
----------

``--validate`` checks the unperturbed guide and representative one-strip,
equal-pair and unequal-pair cells. One variation reduces the maximum mesh
size by 30 percent and increases wavelength resolution; a separate variation
increases both air clearance and PML thickness by 50 percent. It records
absolute changes in normalized alpha/beta and relative changes in alpha.
These are sensitivity checks, not an automatic convergence certificate.
Completed successful checks are checkpointed and skipped when validation
is resumed; failed checks are retried.

Run the lightweight geometry/interpolation checks with:

.. code-block:: sh

   python -m pytest tests/test_image_line_pilot.py -q
