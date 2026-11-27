fem.periodic user API
=====================

Version 2.1.0. This reference covers the deliberately supported user API.
Implementation helpers are documented in their source modules.

Workflow and units
------------------

All lengths are in metres, frequencies in hertz, and constitutive values are relative.
Construction and configuration use keyword arguments. Call ``mesh()`` to build
the initial mesh, ``solve()`` to obtain and store a typed result, and ``show()``
to inspect it interactively. ``mesh_data`` and ``result`` are initially None.
Geometry edits invalidate both; remeshing invalidates the result. Automatic
meshing reuses the last explicit settings. Calling ``show()`` without a result
raises ``NoResultError``.

Adaptive defaults are ``max_refinements=2`` and ``adaptive_tolerance=0.05``.
Zero refinements performs one solve. ``solve_info`` distinguishes algebraic
residuals from adaptive discretization residuals and records stopping reasons.
Python mode and case indices are zero-based.

``solve()`` and frequency sweeps neither save nor open windows. Call
``result.save(path)`` explicitly. ``load_result(path)`` returns inspection-ready
results without solving. Archives use ``fem-results`` schema ``1.0``;
old or convention-incompatible files raise ``PersistenceError``. Loaded
results can plot, show, and save; they do not restart a solver or restore callbacks.
``plot()`` returns a PyVista plotter for full E/H vector views or a Matplotlib Figure for other views without opening a window.
``show(block=True)`` waits for the viewer; ``block=False`` returns immediately.

Select ``E vector`` or ``H vector`` for a rotatable plot of all three field components.

The time convention is ``exp(+i*omega*t)`` with guided propagation
``exp(-i*beta*z)``. Passive materials have nonpositive imaginary constitutive
values; passive forward attenuation is ``-Im(beta)``.

Materials support scalar or diagonal epsilon and mu. Unsupported off-diagonal
forms raise an explicit configuration error. Periodic archives store periodic
field envelopes; the Bloch phase is recorded separately from these fields.

Supported exports
-----------------

``NotDiscretizedError``, ``PeriodicMode``, ``PeriodicModeSet``, ``PeriodicModeSolver2D``, ``PeriodicModeSolver3D``, ``PeriodicModesError``, ``PeriodicSampledFields``, ``PeriodicSweepResult``, ``StaleDiscretizationError``.

Solver construction and operations
----------------------------------

``PeriodicModeSolver2D``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D(*, frequency: 'float', x_range: 'float | Sequence[float]', z_range: 'float | Sequence[float]', polarization: 'str' = 'both', background_material: 'materials.Material' = Material(name='vacuum', epsilon=(1+0j), mu=(1+0j)), boundary: 'materials.IdealBoundary' = IdealBoundary(name='PEC', kind='pec')) -> 'None'

Solve complex Floquet propagation constants in one ``x-z`` unit cell.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``frequency``
     - ``float``
     - Required
     - ``—``
     - Operating frequency in hertz; finite and positive.
   * - ``x_range``
     - ``float | Sequence[float]``
     - Required
     - ``—``
     - Physical x extent or increasing bounds, in metres.
   * - ``z_range``
     - ``float | Sequence[float]``
     - Required
     - ``—``
     - Physical z extent or increasing bounds, in metres.
   * - ``polarization``
     - ``str``
     - Optional
     - ``'both'``
     - TE, TM, or both for the 2D periodic solver; result filtering selects the named polarization.
   * - ``background_material``
     - ``materials.Material``
     - Optional
     - ``Material(name='vacuum', epsilon=(1+0j), mu=(1+0j))``
     - Predefined bulk Material assigned to unfilled space.
   * - ``boundary``
     - ``materials.IdealBoundary``
     - Optional
     - ``IdealBoundary(name='PEC', kind='pec')``
     - Predefined PEC or PMC exterior-boundary material.

Returns: a configured ``PeriodicModeSolver2D``.

``PeriodicModeSolver2D.add_geometry``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.add_geometry(*, shape, material, name=None, clip=False)

Assign a predefined material to a continuous shape in metres.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``shape``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined fem.common.shapes object in metres.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.
   * - ``name``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.
   * - ``clip``
     - ``bool``
     - Optional
     - ``False``
     - Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.add_rectangle``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.add_rectangle(*, x_range, z_range, material, name=None, clip=False)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``x_range``
     - ``float | tuple[float, float] / m``
     - Required
     - ``—``
     - Physical x extent or increasing bounds, in metres.
   * - ``z_range``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Physical z extent or increasing bounds, in metres.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.
   * - ``name``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.
   * - ``clip``
     - ``bool``
     - Optional
     - ``False``
     - Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.add_circle``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.add_circle(*, center, radius, material, name=None, clip=False)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``center``
     - ``tuple[float, ...] / m``
     - Required
     - ``—``
     - Physical centre coordinates in metres.
   * - ``radius``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Positive radius in metres.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.
   * - ``name``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.
   * - ``clip``
     - ``bool``
     - Optional
     - ``False``
     - Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.add_polygon``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.add_polygon(*, points, material, name=None, clip=False)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``points``
     - ``sequence[tuple[float, ...]] / m``
     - Required
     - ``—``
     - Ordered polygon vertex coordinates in metres.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.
   * - ``name``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.
   * - ``clip``
     - ``bool``
     - Optional
     - ``False``
     - Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.set_material``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.set_material(*, geometry, material)

Reassign a predefined material and invalidate mesh/result.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``geometry``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A handle returned by add_geometry or a geometry convenience method.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.set_shape``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.set_shape(*, geometry, shape)

Replace a shape in metres and invalidate mesh/result.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``geometry``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A handle returned by add_geometry or a geometry convenience method.
   * - ``shape``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined fem.common.shapes object in metres.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.set_boundary``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.set_boundary(*, material)

Set the exterior PEC/PMC wall; internal SIBC uses geometry assignment.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.add_mesh_refinement``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.add_mesh_refinement(*, shape: 'Shape2D', max_element_size: 'float', name: 'str | None' = None) -> 'MeshRefinement'

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``shape``
     - ``Shape2D``
     - Required
     - ``—``
     - A predefined fem.common.shapes object in metres.
   * - ``max_element_size``
     - ``float``
     - Required
     - ``—``
     - Maximum initial element edge length in metres.
   * - ``name``
     - ``str | None``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.add_pml``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.add_pml(*, thickness: "float | Literal['auto']" = 'auto', order: 'int' = 3, sigma_max: 'float' = 5.0, direction: 'str' = 'all') -> 'PMLSpec'

Add a polynomial PML on selected exterior sides. Auto thickness is one third wavelength.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``thickness``
     - ``float | Literal['auto']``
     - Optional
     - ``'auto'``
     - Thickness in metres inside each selected outer face. Auto uses one third of the free-space wavelength.
   * - ``order``
     - ``int``
     - Optional
     - ``3``
     - Polynomial order of the PML profile.
   * - ``sigma_max``
     - ``float``
     - Optional
     - ``5.0``
     - Backend PML-strength magnitude; the outgoing stretch has negative imaginary sign.
   * - ``direction``
     - ``str``
     - Optional
     - ``'all'``
     - Propagation direction for solve; selected coordinate direction for PML.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.remove``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.remove(*, geometry)

Remove an object and invalidate its mesh/result.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``geometry``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A handle returned by add_geometry or a geometry convenience method.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.refine``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.refine(factor: 'float' = 2.0) -> 'FEMPeriodicMesh2D'

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``factor``
     - ``float``
     - Optional
     - ``2.0``
     - Mesh refinement factor greater than one.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver2D.mesh``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.mesh(*, geometry_tolerance: 'float | None' = None, pml_layers: 'int' = 3, max_element_size: 'float | None' = None, resolution: 'tuple[int, int] | None' = None, wavelength_elements: 'int' = 4, element_order: 'int' = 2, quadrature_order: 'int' = 5)

Build the initial mesh; solve() may subsequently refine it.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``geometry_tolerance``
     - ``float | None``
     - Optional
     - ``None``
     - Target chord sagitta on curved CAD boundaries, in metres. None uses curvature-based sizing.
   * - ``pml_layers``
     - ``int``
     - Optional
     - ``3``
     - Number of partitions across each PML thickness.
   * - ``max_element_size``
     - ``float | None``
     - Optional
     - ``None``
     - Maximum initial element edge length in metres.
   * - ``resolution``
     - ``tuple[int, int] | None``
     - Optional
     - ``None``
     - Initial node counts; use instead of a maximum element size.
   * - ``wavelength_elements``
     - ``int``
     - Optional
     - ``4``
     - Minimum number of initial elements per local wavelength.
   * - ``element_order``
     - ``int``
     - Optional
     - ``2``
     - Finite-element polynomial order supported by this backend.
   * - ``quadrature_order``
     - ``int``
     - Optional
     - ``5``
     - Polynomial integration degree for finite-element assembly. Defaults to five.

Returns: the initial mesh stored in ``mesh_data``.

``PeriodicModeSolver2D.configure_adaptivity``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.configure_adaptivity(*, minimum_quality=None, max_growth=None, stagnation_window=3, minimum_reduction=0.01, anisotropic=False)

Set candidate simplex quality, per-pass growth and residual progress limits.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``minimum_quality``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional lower bound on candidate simplex mean-ratio quality.
   * - ``max_growth``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional maximum candidate-to-previous element-count ratio.
   * - ``stagnation_window``
     - ``int``
     - Optional
     - ``3``
     - Consecutive insufficient reductions before stopping. None disables the check.
   * - ``minimum_reduction``
     - ``float``
     - Optional
     - ``0.01``
     - Minimum fractional residual reduction per refinement for the stagnation check.
   * - ``anisotropic``
     - ``bool``
     - Optional
     - ``False``
     - Use a directional MMG metric based on local residual-density gradients in nonperiodic 2D and 3D.

Returns: None.

``PeriodicModeSolver2D.solve``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.solve(*, neff_guess: 'complex | None' = None, num_modes: 'int' = 4, direction: 'Direction' = 'forward', eigensolver: 'str' = 'auto', backend: 'str' = 'auto', eigensolver_tolerance: 'float' = 1e-10, residual_tolerance: 'float' = 1e-08, propagation_ratio_tolerance: 'float' = 0.001, max_pml_fraction: 'float | None' = 0.5, dense_linearization_limit: 'int' = 700, max_refinements: 'int' = 2, adaptive_tolerance: 'float' = 0.05, max_elements: 'int | None' = None, max_dofs: 'int | None' = None)

Solve and return modes with separate algebraic and adaptive controls.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``neff_guess``
     - ``complex | None``
     - Optional
     - ``None``
     - Dimensionless complex effective-index search target.
   * - ``num_modes``
     - ``int``
     - Optional
     - ``4``
     - Number of modes requested; positive integer.
   * - ``direction``
     - ``Direction``
     - Optional
     - ``'forward'``
     - Propagation direction for solve; selected coordinate direction for PML.
   * - ``eigensolver``
     - ``str``
     - Optional
     - ``'auto'``
     - Eigenproblem algorithm: auto, dense, or sparse where supported.
   * - ``backend``
     - ``str``
     - Optional
     - ``'auto'``
     - Execution backend: auto (default) selects hybrid GPU Arnoldi/CPU Ritz when available, otherwise warns and selects SciPy CPU; hybrid, gpu (resident CUDA Ritz), and cpu explicitly select a backend.
   * - ``eigensolver_tolerance``
     - ``float``
     - Optional
     - ``1e-10``
     - Algebraic eigensolver convergence tolerance.
   * - ``residual_tolerance``
     - ``float``
     - Optional
     - ``1e-08``
     - Maximum accepted eigenproblem residual, separate from adaptation.
   * - ``propagation_ratio_tolerance``
     - ``float``
     - Optional
     - ``0.001``
     - Minimum real-power fraction used to classify a mode as propagating.
   * - ``max_pml_fraction``
     - ``float | None``
     - Optional
     - ``0.5``
     - Maximum permitted modal PML fraction; None disables PML-mode filtering.
   * - ``dense_linearization_limit``
     - ``int``
     - Optional
     - ``700``
     - Largest linearized eigenproblem dimension allowed for dense eigensolution.
   * - ``max_refinements``
     - ``int``
     - Optional
     - ``2``
     - Maximum mesh refinements after the initial solve; zero means one solve.
   * - ``adaptive_tolerance``
     - ``float``
     - Optional
     - ``0.05``
     - Relative discretization-residual stopping threshold.
   * - ``max_elements``
     - ``int | None``
     - Optional
     - ``None``
     - Adaptive mesh element budget.
   * - ``max_dofs``
     - ``int | None``
     - Optional
     - ``None``
     - Maximum field-space size on an adaptive candidate mesh before eliminating boundary DOFs.

Returns: the physics-specific result stored in ``result``.

``PeriodicModeSolver2D.sweep``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.sweep(*, frequencies, max_refinements=2, adaptive_tolerance=0.05, **solve_options)

Yield mode sets in increasing frequency order, adapting only at the highest frequency.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``frequencies``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Strictly increasing positive frequencies in hertz.
   * - ``max_refinements``
     - ``int``
     - Optional
     - ``2``
     - Maximum mesh refinements after the initial solve; zero means one solve.
   * - ``adaptive_tolerance``
     - ``float``
     - Optional
     - ``0.05``
     - Relative discretization-residual stopping threshold.
   * - ``solve_options``
     - ``array-like or scalar``
     - Optional
     - ``—``
     - Additional keyword arguments accepted by this solver's solve method.

Returns: an iterator of mode sets in increasing frequency order.

``PeriodicModeSolver2D.show``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver2D.show(*, block: 'bool' = True)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``block``
     - ``bool``
     - Optional
     - ``True``
     - Wait for the interactive viewer to close when true.

Returns: the interactive viewer window.

``PeriodicModeSolver3D``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D(*, frequency: 'float', x_range: 'float | tuple[float, float]', y_range: 'float | tuple[float, float]', z_range: 'float | tuple[float, float]', background_material: 'materials.Material' = Material(name='vacuum', epsilon=(1+0j), mu=(1+0j)), boundary: 'materials.IdealBoundary' = IdealBoundary(name='PEC', kind='pec')) -> 'None'

Solve complex fixed-frequency Bloch propagation in a tetrahedral cell.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``frequency``
     - ``float``
     - Required
     - ``—``
     - Operating frequency in hertz; finite and positive.
   * - ``x_range``
     - ``float | tuple[float, float]``
     - Required
     - ``—``
     - Physical x extent or increasing bounds, in metres.
   * - ``y_range``
     - ``float | tuple[float, float]``
     - Required
     - ``—``
     - Physical y extent or increasing bounds, in metres.
   * - ``z_range``
     - ``float | tuple[float, float]``
     - Required
     - ``—``
     - Physical z extent or increasing bounds, in metres.
   * - ``background_material``
     - ``materials.Material``
     - Optional
     - ``Material(name='vacuum', epsilon=(1+0j), mu=(1+0j))``
     - Predefined bulk Material assigned to unfilled space.
   * - ``boundary``
     - ``materials.IdealBoundary``
     - Optional
     - ``IdealBoundary(name='PEC', kind='pec')``
     - Predefined PEC or PMC exterior-boundary material.

Returns: a configured ``PeriodicModeSolver3D``.

``PeriodicModeSolver3D.add_geometry``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.add_geometry(*, shape, material, name=None, clip=False)

Assign a predefined material to a continuous shape in metres.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``shape``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined fem.common.shapes object in metres.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.
   * - ``name``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.
   * - ``clip``
     - ``bool``
     - Optional
     - ``False``
     - Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.add_box``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.add_box(*, x_range, y_range, z_range, material, name=None, clip=False)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``x_range``
     - ``float | tuple[float, float] / m``
     - Required
     - ``—``
     - Physical x extent or increasing bounds, in metres.
   * - ``y_range``
     - ``float | tuple[float, float] / m``
     - Required
     - ``—``
     - Physical y extent or increasing bounds, in metres.
   * - ``z_range``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Physical z extent or increasing bounds, in metres.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.
   * - ``name``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.
   * - ``clip``
     - ``bool``
     - Optional
     - ``False``
     - Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.add_sphere``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.add_sphere(*, center, radius, material, name=None, clip=False)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``center``
     - ``tuple[float, ...] / m``
     - Required
     - ``—``
     - Physical centre coordinates in metres.
   * - ``radius``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Positive radius in metres.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.
   * - ``name``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.
   * - ``clip``
     - ``bool``
     - Optional
     - ``False``
     - Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.add_cylinder``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.add_cylinder(*, center, radius, z_range, material, name=None, clip=False)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``center``
     - ``tuple[float, ...] / m``
     - Required
     - ``—``
     - Physical centre coordinates in metres.
   * - ``radius``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Positive radius in metres.
   * - ``z_range``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Physical z extent or increasing bounds, in metres.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.
   * - ``name``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.
   * - ``clip``
     - ``bool``
     - Optional
     - ``False``
     - Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.set_material``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.set_material(*, geometry, material)

Reassign a predefined material and invalidate mesh/result.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``geometry``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A handle returned by add_geometry or a geometry convenience method.
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.set_shape``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.set_shape(*, geometry, shape)

Replace a shape in metres and invalidate mesh/result.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``geometry``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A handle returned by add_geometry or a geometry convenience method.
   * - ``shape``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined fem.common.shapes object in metres.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.set_boundary``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.set_boundary(*, material)

Set the exterior PEC/PMC wall; internal SIBC uses geometry assignment.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``material``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A predefined bulk, ideal-boundary, or supported SIBC material.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.add_mesh_refinement``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.add_mesh_refinement(*, shape: 'Shape3D', max_element_size: 'float', name: 'str | None' = None) -> 'object'

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``shape``
     - ``Shape3D``
     - Required
     - ``—``
     - A predefined fem.common.shapes object in metres.
   * - ``max_element_size``
     - ``float``
     - Required
     - ``—``
     - Maximum initial element edge length in metres.
   * - ``name``
     - ``str | None``
     - Optional
     - ``None``
     - Optional name used for later identification and diagnostics.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.add_pml``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.add_pml(*, thickness: "float | Literal['auto']" = 'auto', order: 'int' = 3, sigma_max: 'float' = 5.0, direction: 'str' = 'all') -> 'PMLSpec'

Add a polynomial PML on selected exterior sides. Auto thickness is one third wavelength.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``thickness``
     - ``float | Literal['auto']``
     - Optional
     - ``'auto'``
     - Thickness in metres inside each selected outer face. Auto uses one third of the free-space wavelength.
   * - ``order``
     - ``int``
     - Optional
     - ``3``
     - Polynomial order of the PML profile.
   * - ``sigma_max``
     - ``float``
     - Optional
     - ``5.0``
     - Backend PML-strength magnitude; the outgoing stretch has negative imaginary sign.
   * - ``direction``
     - ``str``
     - Optional
     - ``'all'``
     - Propagation direction for solve; selected coordinate direction for PML.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.remove``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.remove(*, geometry)

Remove an object and invalidate its mesh/result.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``geometry``
     - ``array-like or scalar``
     - Required
     - ``—``
     - A handle returned by add_geometry or a geometry convenience method.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.refine``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.refine(factor: 'float' = 2.0) -> 'PeriodicMesh3D'

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``factor``
     - ``float``
     - Optional
     - ``2.0``
     - Mesh refinement factor greater than one.

Returns: the configured geometry/excitation handle, or None for in-place configuration.

``PeriodicModeSolver3D.mesh``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.mesh(*, geometry_tolerance: 'float | None' = None, max_element_size: 'float | None' = None, wavelength_elements: 'int' = 4, material_aware: 'bool' = True, element_order: 'int' = 2, quadrature_order: 'int' = 5)

Build the initial mesh; solve() may subsequently refine it.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``geometry_tolerance``
     - ``float | None``
     - Optional
     - ``None``
     - Target chord sagitta on curved CAD boundaries, in metres. None uses curvature-based sizing.
   * - ``max_element_size``
     - ``float | None``
     - Optional
     - ``None``
     - Maximum initial element edge length in metres.
   * - ``wavelength_elements``
     - ``int``
     - Optional
     - ``4``
     - Minimum number of initial elements per local wavelength.
   * - ``material_aware``
     - ``bool``
     - Optional
     - ``True``
     - Use material-dependent initial mesh sizing.
   * - ``element_order``
     - ``int``
     - Optional
     - ``2``
     - Finite-element polynomial order supported by this backend.
   * - ``quadrature_order``
     - ``int``
     - Optional
     - ``5``
     - Polynomial integration degree for finite-element assembly. Defaults to five.

Returns: the initial mesh stored in ``mesh_data``.

``PeriodicModeSolver3D.configure_adaptivity``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.configure_adaptivity(*, minimum_quality=None, max_growth=None, stagnation_window=3, minimum_reduction=0.01, anisotropic=False)

Set candidate simplex quality, per-pass growth and residual progress limits.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``minimum_quality``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional lower bound on candidate simplex mean-ratio quality.
   * - ``max_growth``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional maximum candidate-to-previous element-count ratio.
   * - ``stagnation_window``
     - ``int``
     - Optional
     - ``3``
     - Consecutive insufficient reductions before stopping. None disables the check.
   * - ``minimum_reduction``
     - ``float``
     - Optional
     - ``0.01``
     - Minimum fractional residual reduction per refinement for the stagnation check.
   * - ``anisotropic``
     - ``bool``
     - Optional
     - ``False``
     - Use a directional MMG metric based on local residual-density gradients in nonperiodic 2D and 3D.

Returns: None.

``PeriodicModeSolver3D.solve``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.solve(*, neff_guess: 'complex | None' = None, num_modes: 'int' = 4, direction: 'Direction' = 'forward', eigensolver: 'str' = 'auto', backend: 'str' = 'auto', eigensolver_tolerance: 'float' = 1e-10, residual_tolerance: 'float' = 1e-08, divergence_tolerance: 'float' = 1e-06, propagation_ratio_tolerance: 'float' = 0.001, max_pml_fraction: 'float | None' = 0.5, dense_linearization_limit: 'int' = 700, ncv: 'int | None' = None, max_iterations: 'int | None' = None, random_seed: 'int' = 0, max_refinements: 'int' = 2, adaptive_tolerance: 'float' = 0.05, max_elements: 'int | None' = None, max_dofs: 'int | None' = None)

Solve and return modes with separate algebraic and adaptive controls.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``neff_guess``
     - ``complex | None``
     - Optional
     - ``None``
     - Dimensionless complex effective-index search target.
   * - ``num_modes``
     - ``int``
     - Optional
     - ``4``
     - Number of modes requested; positive integer.
   * - ``direction``
     - ``Direction``
     - Optional
     - ``'forward'``
     - Propagation direction for solve; selected coordinate direction for PML.
   * - ``eigensolver``
     - ``str``
     - Optional
     - ``'auto'``
     - Eigenproblem algorithm: auto, dense, or sparse where supported.
   * - ``backend``
     - ``str``
     - Optional
     - ``'auto'``
     - Execution backend: auto (default) selects hybrid GPU Arnoldi/CPU Ritz when available, otherwise warns and selects SciPy CPU; hybrid, gpu (resident CUDA Ritz), and cpu explicitly select a backend.
   * - ``eigensolver_tolerance``
     - ``float``
     - Optional
     - ``1e-10``
     - Algebraic eigensolver convergence tolerance.
   * - ``residual_tolerance``
     - ``float``
     - Optional
     - ``1e-08``
     - Maximum accepted eigenproblem residual, separate from adaptation.
   * - ``divergence_tolerance``
     - ``float``
     - Optional
     - ``1e-06``
     - Maximum accepted discrete Gauss-law residual.
   * - ``propagation_ratio_tolerance``
     - ``float``
     - Optional
     - ``0.001``
     - Minimum real-power fraction used to classify a mode as propagating.
   * - ``max_pml_fraction``
     - ``float | None``
     - Optional
     - ``0.5``
     - Maximum permitted modal PML fraction; None disables PML-mode filtering.
   * - ``dense_linearization_limit``
     - ``int``
     - Optional
     - ``700``
     - Largest linearized eigenproblem dimension allowed for dense eigensolution.
   * - ``ncv``
     - ``int | None``
     - Optional
     - ``None``
     - Arnoldi subspace dimension; None lets the eigensolver choose.
   * - ``max_iterations``
     - ``int | None``
     - Optional
     - ``None``
     - CPU: maximum ARPACK iterations. GPU: maximum Arnoldi operator applications; None allows thirteen subspace cycles.
   * - ``random_seed``
     - ``int``
     - Optional
     - ``0``
     - Seed for reproducible eigensolver starting vectors.
   * - ``max_refinements``
     - ``int``
     - Optional
     - ``2``
     - Maximum mesh refinements after the initial solve; zero means one solve.
   * - ``adaptive_tolerance``
     - ``float``
     - Optional
     - ``0.05``
     - Relative discretization-residual stopping threshold.
   * - ``max_elements``
     - ``int | None``
     - Optional
     - ``None``
     - Adaptive mesh element budget.
   * - ``max_dofs``
     - ``int | None``
     - Optional
     - ``None``
     - Maximum field-space size on an adaptive candidate mesh before eliminating boundary DOFs.

Returns: the physics-specific result stored in ``result``.

``PeriodicModeSolver3D.sweep``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.sweep(*, frequencies, max_refinements=2, adaptive_tolerance=0.05, **solve_options)

Yield mode sets in increasing frequency order, adapting only at the highest frequency.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``frequencies``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Strictly increasing positive frequencies in hertz.
   * - ``max_refinements``
     - ``int``
     - Optional
     - ``2``
     - Maximum mesh refinements after the initial solve; zero means one solve.
   * - ``adaptive_tolerance``
     - ``float``
     - Optional
     - ``0.05``
     - Relative discretization-residual stopping threshold.
   * - ``solve_options``
     - ``array-like or scalar``
     - Optional
     - ``—``
     - Additional keyword arguments accepted by this solver's solve method.

Returns: an iterator of mode sets in increasing frequency order.

``PeriodicModeSolver3D.show``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSolver3D.show(*, block: 'bool' = True)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``block``
     - ``bool``
     - Optional
     - ``True``
     - Wait for the interactive viewer to close when true.

Returns: the interactive viewer window.

Returned results
----------------

Result objects are returned by solving or loading. Their constructors are
implementation details. Inspect ``mesh_data``, ``metadata``, ``solve_info``, and
``frequency`` where applicable. Modal sets support iteration, indexing, ``neff``
and ``beta``; each mode provides sampled fields and residual diagnostics.

``PeriodicModeSet.plot``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSet.plot(*, component: 'str | None' = None, quantity: 'str' = 'real', mode: 'int' = 0)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``component``
     - ``str | None``
     - Optional
     - ``None``
     - Field component to display, such as Ey; electrostatics also accepts potential or mesh.
   * - ``quantity``
     - ``str``
     - Optional
     - ``'real'``
     - Displayed field quantity: real, imag, magnitude/abs, or phase; static fields support real or magnitude.
   * - ``mode``
     - ``int``
     - Optional
     - ``0``
     - Zero-based mode index.

Returns: a PyVista plotter for full E/H vector views or a Matplotlib Figure for other views.

``PeriodicModeSet.show``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSet.show(*, block: 'bool' = True)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``block``
     - ``bool``
     - Optional
     - ``True``
     - Wait for the interactive viewer to close when true.

Returns: the interactive viewer window.

``PeriodicModeSet.save``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSet.save(path)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``path``
     - ``str | PathLike``
     - Required
     - ``—``
     - Destination/source HDF5 path. Saving is atomic; loading does not run a solver.

Returns: the written Path.

``PeriodicModeSet.mode``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSet.mode(number: 'int') -> 'PeriodicMode'

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``number``
     - ``int``
     - Required
     - ``—``
     - Zero-based mode index.

Returns: the selected mode.

``PeriodicModeSet.by_polarization``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicModeSet.by_polarization(polarization: 'str') -> 'tuple[PeriodicMode, ...]'

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``polarization``
     - ``str``
     - Required
     - ``—``
     - TE, TM, or both for the 2D periodic solver; result filtering selects the named polarization.

Returns: the selected data or diagnostic report.

``PeriodicSweepResult.from_results``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicSweepResult.from_results(results)

Combine solved periodic mode sets into a frequency or parameter sweep.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``results``
     - ``Sequence[PeriodicModeSet]``
     - Required
     - ``—``
     - Nonempty sequence of completed periodic mode sets, in sweep order.

Returns: the selected data or diagnostic report.

``PeriodicSweepResult.plot``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicSweepResult.plot(*, case=0, component='Ey', quantity='real', mode=0)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``case``
     - ``int``
     - Optional
     - ``0``
     - Zero-based sweep case index.
   * - ``component``
     - ``str | None``
     - Optional
     - ``'Ey'``
     - Field component to display, such as Ey; electrostatics also accepts potential or mesh.
   * - ``quantity``
     - ``str``
     - Optional
     - ``'real'``
     - Displayed field quantity: real, imag, magnitude/abs, or phase; static fields support real or magnitude.
   * - ``mode``
     - ``int``
     - Optional
     - ``0``
     - Zero-based mode index.

Returns: a PyVista plotter for full E/H vector views or a Matplotlib Figure for other views.

``PeriodicSweepResult.show``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicSweepResult.show(*, block=True)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``block``
     - ``bool``
     - Optional
     - ``True``
     - Wait for the interactive viewer to close when true.

Returns: the interactive viewer window.

``PeriodicSweepResult.save``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    PeriodicSweepResult.save(path)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``path``
     - ``str | PathLike``
     - Required
     - ``—``
     - Destination/source HDF5 path. Saving is atomic; loading does not run a solver.

Returns: the written Path.

Transverse vector views
-----------------------

In the interactive viewer opened by ``show()``, select ``Et`` or ``Ht``
for electric or magnetic vector arrows. Waveguides display x/y components.
Two-dimensional periodic results display x/z components in the mesh plane,
including in ``plot()``. Three-dimensional periodic results display x/y components.

``real`` and ``imag`` draw the selected phasor part. ``magnitude`` applies one
global phase rotation to make the largest transverse component real, then draws
that real vector field. All arrows are coloured by the complex transverse norm;
their lengths are scaled together within each view. Phase is not offered as a
vector quantity. Mesh and material overlays leave the transverse arrows visible.

One-dimensional waveguide samples lie on y=0. Three-dimensional periodic modes
use PyVista arrows with zero z component.
The same selections work for loaded archives and individual sweep cases.

Result data and diagnostics
---------------------------

``mesh_data.coordinates`` stores physical nodes in metres; ``elements`` stores
zero-based connectivity. ``axes`` identifies physical coordinate order.
``mesh_data.metadata['context']`` records material and boundary configuration.
The result is an inspection snapshot; editing it cannot restart a solver.

Each selected mode exposes ``neff`` (dimensionless), ``beta`` (rad/m),
``fields`` (sampled E/H in V/m and A/m), ``coefficients`` (FEM expansion),
and ``residual`` (algebraic eigenproblem residual). Field component names
include Ex, Ey, Ez, Hx, Hy, and Hz. Coordinate and cell ownership arrays
preserve field locations; they must not be interpreted as interchangeable
nodal or cell values. ``solve_info`` records adaptive history separately.

Periodic fields are Bloch envelopes. ``period`` is in m; each mode also
provides ``gamma``, ``bloch_multiplier``, folded propagation quantities,
and Gauss-law/PML filtering diagnostics. Combine solved cases with
``PeriodicSweepResult.from_results(results)`` and call ``save(path)``.
Loaded multi-case archives index cases lazily.

Geometry and material values
----------------------------

Define reusable materials and shapes with ``fem.common`` before assigning them.
Use ``Material(name=..., epsilon=..., mu=...)`` for bulk media,
``materials.PEC`` or ``materials.PMC`` for ideal boundaries, and the
documented ``materials.copper``-style presets where SIBC is supported.
Continuous primitives and Boolean/transformed shapes live in
``fem.common.shapes``. Solver packages also re-export Material, materials, and shapes.

``fem.load``
~~~~~~~~~~~~

.. code-block:: python

    fem.load(path: 'str | os.PathLike[str]')

Load any FEM result archive; the archive records which solver wrote it.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``path``
     - ``str | os.PathLike[str]``
     - Required
     - ``—``
     - Destination/source HDF5 path. Saving is atomic; loading does not run a solver.

Returns: a typed result; multi-case archives provide lazy case access.

Errors
------

Invalid inputs raise ``ConfigurationError`` or ``GeometryError`` where available.
Mesh and numerical failures raise the corresponding ``MeshError`` or
``SolverError``. ``NoResultError`` requires a successful solve first.
``PersistenceError`` identifies an incompatible or unreadable archive.
Viewer errors include the executable path or installation setting needed to
correct a launch failure. Saving and loading do not require an active GUI.

Backend selection
-----------------

``backend="auto"`` tries hybrid GPU Arnoldi/CPU Ritz, then warns before SciPy CPU fallback.
``hybrid`` explicitly selects that path; ``gpu`` explicitly requires resident CUDA Ritz.
Both accelerated choices use cuDSS for direct scattering solves.
