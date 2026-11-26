FEM adaptivity policy
=====================

Version 1.0.0 of the shared adaptive policy included in FEM. This library supplies error
estimation and refinement policy to the FEM families. Assembly remains in
the individual solvers.

Users configure adaptation through their solver's ``solve`` method::

    result = solver.solve(max_refinements=2, adaptive_tolerance=0.05)

Zero refinements means one solve on the initial mesh. Results expose stopping
status, element counts, and discretization residuals through ``solve_info``.
Algebraic residuals remain separate. Mesh and numerical failures are visible.

See the `getting started <getting_started.md>`_ for installation and
the `API reference <api/adaptivity.rst>`_ for the supported user surface.

Mesh and estimator controls
---------------------------

Waveguide and periodic modal solvers halve the local size target on cells
carrying half the combined squared indicator. Gmsh uses those targets to
generate a conforming candidate mesh from the CAD model.
Nonperiodic 2D and 3D adaptive passes use MMG2D and MMG3D to remesh each
CAD region's interior. Both neighboring regions keep the same boundary nodes
and facets, including conductor walls, port apertures and far-field monitors.
Initial CAD meshing and periodic adaptation use Gmsh. One-dimensional
adaptation uses native interval bisection. Driven absorber-only passes keep
their layer-doubling algorithm.

Each attempted refinement records ``refinement_backend`` in its starting
history entry. If MMG fails, changes a protected interface, or produces cells
below the configured quality floor, the solver
retries using Gmsh and records ``fallback_reason``. The same element, DOF,
and mesh-quality limits apply. An exceeded resource or growth limit stops
refinement and retains the accepted mesh. Use
``solver.configure_adaptivity(anisotropic=True)`` to request a tensor metric
in MMG. Its shortest direction follows the local gradient of the squared
residual indicator per unit cell measure. The maximum directional size
ratio is two. The default uses scalar size targets.

Periodic nodes share the smaller
target across each seam and the new mesh supplies a fresh constraint map.
Manual refinement and failed-eigensolve recovery scale the local sizes on
the current mesh, including regions refined in previous passes.

For curved 2D meshes and periodic 3D meshes, ``mesh(geometry_tolerance=...)`` sets a physical chord-sagitta
target in metres using sampled CAD curvature. The default curvature setting
requests 24 segments per circle. Electrostatic adaptation regenerates the CAD interfaces in 2D and 3D. In periodic 2D, ``pml_layers`` sets the
number of normal subdivisions per absorbing layer (default three).

Modal and driven solves accept ``max_elements`` and ``max_dofs`` as optional
limits on refinement candidates. DOFs count the field space before constraints
and include the absorbing shell for driven problems. The last accepted solve
is returned when the next mesh exceeds a limit. The history records the
stopping reason, edge extrema, and minimum and fifth-percentile simplex
mean-ratio quality (one for a regular simplex, approaching zero at collapse).

For multidimensional Maxwell problems, ``residual_components`` reports
normalized interface, element-equation, divergence (where applicable), and
natural-boundary contributions for the mode or excitation with the largest
indicator. Interior constitutive values come from the CAD cell tag. The
volume terms differentiate local polynomial reconstructions, including within
a PML. Modal histories with a PML also partition the squared indicator into
physical-domain and absorbing-layer contributions. Use the field error, charge, energy or transmission relevant to the
problem to calibrate an estimator threshold. Each modal pass records the
effective indices. Electrostatic histories record energy and conductor charge.

Driven scattering and antenna solves report separate residuals for the physical
domain, the PML volume and the entrance to the PML. Once the physical residual
meets tolerance, the solver holds the physical mesh fixed and doubles the
number of shell layers. By default it compares scattering width or RCS in
36 directions, or the full complex antenna S-matrix. Pass
``pml_observable=lambda result: ...`` to track another scalar or array.
The relative change is the norm of the difference divided by the larger of
the two observable norms. Set ``pml_tolerance`` to control this check. Its
default is ``adaptive_tolerance``.

Shell doublings count toward ``max_refinements`` and have a separate cap,
``max_pml_refinements`` (default three). An excited shell needs two solutions
on the same physical mesh to establish convergence. ``pml_limit`` indicates
that the physical residual met tolerance and the shell comparison still needs
more refinement. Each history entry records ``physical_converged``,
``pml_converged`` and the sampled observable. This check measures sensitivity to the shell discretization
at the configured thickness and absorption profile.

Use ``solver.configure_adaptivity(minimum_quality=..., max_growth=...)`` to
reject a candidate mesh before assembly. For driven problems the quality floor
applies to the physical mesh and the growth limit includes the shell. By
default three consecutive physical refinements with less than one percent
residual reduction stop the solve with ``stagnation``. Set
``stagnation_window=None`` to disable this check.

Examples
--------

Configure adaptation through a solver's ``solve(max_refinements=2,
adaptive_tolerance=0.05)``. The FEM examples under ``examples/fem/`` demonstrate the
physical solver workflows; ``max_refinements=0`` selects a fixed initial mesh.

Tests
-----

Adaptive convergence, budget exhaustion, failures, and state restoration are
tested across solver families in ``tests/test_adaptivity.py`` at the repository
root. Package-specific residual estimators are additionally exercised by their
own numerical suites.
