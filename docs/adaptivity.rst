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
carrying half the combined squared indicator. A Gmsh background field
carries the targets to the next CAD mesh. Periodic nodes share the smaller
target across each seam and the new mesh supplies a fresh constraint map.
A failed eigensolve can request a global refinement.

For 2D meshes, ``mesh(geometry_tolerance=...)`` sets a physical chord-sagitta
target in metres using sampled CAD curvature. The default curvature setting
requests 24 segments per circle. Electrostatic adaptation remeshes curved CAD
and bisects straight-sided meshes. In periodic 2D, ``pml_layers`` sets the
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
problem to calibrate an estimator threshold.

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
