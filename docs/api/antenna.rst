fem.antenna user API
====================

Version 2.0.2. This reference covers the deliberately supported user API.
Implementation helpers are documented in their source modules.

Waveguide-fed antennas
----------------------

See the `solver guide <../solvers/antenna.rst>`_ for field conventions and examples.

``AntennaSolver3D``, ``AntennaResult``, ``AntennaSweep``, ``RadiationPattern``.

``AntennaSolver3D``
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D(*, frequency, x_range, y_range, z_range, pml_thickness=None, pml_layers=5)

Driven antenna in a physical vacuum box surrounded by a Cartesian PML.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``frequency``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Operating frequency in hertz; finite and positive.
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
   * - ``pml_thickness``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Exterior PML thickness per side in metres. Defaults to one third of the free-space wavelength at the constructor frequency.
   * - ``pml_layers``
     - ``int``
     - Optional
     - ``5``
     - Number of mesh layers normal to each face of the physical box.

Returns: the configured solver.

``AntennaSolver3D.add_geometry``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.add_geometry(*, shape, material)

Add an isotropic bulk, PEC, or impedance target inside the radiation box.

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
     - Isotropic bulk material, PEC, surface impedance or good conductor.

Returns: None.

``AntennaSolver3D.add_mesh_refinement``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.add_mesh_refinement(*, shape, max_element_size)

Set a local size in a rectangle, box, circle or sphere in physical coordinates.

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
   * - ``max_element_size``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Maximum initial element edge length in metres.

Returns: None.

``AntennaSolver3D.add_pml``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.add_pml(*, thickness: "float | Literal['auto']" = 'auto', order: 'int' = 3, sigma_max: 'float' = 5.0, direction: 'str' = 'all') -> 'PMLSpec'

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
     - One face, an axis pair, or all exterior faces.

Returns: the PML profile.

``AntennaSolver3D.add_nf2ff``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.add_nf2ff(distance='auto', *, x_range=None, y_range=None, z_range=None)

Place a closed vacuum NF2FF surface by its inset from the PML or by axis ranges.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``distance``
     - ``str``
     - Optional
     - ``'auto'``
     - Inset in metres from each PML interface, or 'auto' for a wavelength-based vacuum clearance.
   * - ``x_range``
     - ``float | tuple[float, float] / m``
     - Optional
     - ``None``
     - Physical x extent or increasing bounds, in metres.
   * - ``y_range``
     - ``float | tuple[float, float] / m``
     - Optional
     - ``None``
     - Physical y extent or increasing bounds, in metres.
   * - ``z_range``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Physical z extent or increasing bounds, in metres.

Returns: None.

``AntennaSolver3D.add_waveguide_port``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.add_waveguide_port(*, name, axis, position, cross_section, mode=1, num_modes=None, amplitude=1.0, phase_deg=0.0)

Assign a numerical wave port. Mode 1 is the highest effective-index forward mode.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``name``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Unique port name.
   * - ``axis``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Waveguide axis (x, y or z).
   * - ``position``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Reference plane coordinate in metres along the guide axis.
   * - ``cross_section``
     - ``array-like or scalar``
     - Required
     - ``—``
     - 2D shape selecting the conforming aperture.
   * - ``mode``
     - ``int``
     - Optional
     - ``1``
     - One-based effective-index rank of the launched mode.
   * - ``num_modes``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Retained mode count. Defaults to the selected mode rank.
   * - ``amplitude``
     - ``float``
     - Optional
     - ``1.0``
     - Incident power-wave magnitude in square roots of watts.
   * - ``phase_deg``
     - ``float``
     - Optional
     - ``0.0``
     - Phase in degrees for the drive or displayed vector phasor.

Returns: None.

``AntennaSolver3D.add_discrete_port``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.add_discrete_port(*, name, axis, gap, impedance=50.0, amplitude=1.0, phase_deg=0.0)

Assign a box or cylindrical feed gap between conductors with a real reference impedance.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``name``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Unique port name.
   * - ``axis``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Waveguide axis (x, y or z).
   * - ``gap``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Box or cylinder spanning the two feed terminals along the port axis.
   * - ``impedance``
     - ``float``
     - Optional
     - ``50.0``
     - Positive real reference impedance in ohms.
   * - ``amplitude``
     - ``float``
     - Optional
     - ``1.0``
     - Incident power-wave magnitude in square roots of watts.
   * - ``phase_deg``
     - ``float``
     - Optional
     - ``0.0``
     - Phase in degrees for the drive or displayed vector phasor.

Returns: None.

``AntennaSolver3D.solve_port_modes``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.solve_port_modes(*, port, frequency=None)

Compute the selected port's 2D ModeSet for inspection before a driven solve.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``port``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Name of the configured waveguide port.
   * - ``frequency``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Operating frequency in hertz; finite and positive.

Returns: a 2D waveguide ModeSet.

``AntennaSolver3D.set_excitation``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.set_excitation(*, port, amplitude, phase_deg=0.0)

Set the incident power-wave magnitude in sqrt(W) and phase in degrees.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``port``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Name of the configured waveguide port.
   * - ``amplitude``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Incident power-wave magnitude in square roots of watts.
   * - ``phase_deg``
     - ``float``
     - Optional
     - ``0.0``
     - Phase in degrees for the drive or displayed vector phasor.

Returns: None.

``AntennaSolver3D.mesh``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.mesh(*, max_element_size=None, element_order=2, quadrature_order=5)

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``max_element_size``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Maximum initial element edge length in metres.
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

Returns: the physical mesh.

``AntennaSolver3D.configure_adaptivity``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.configure_adaptivity(*, minimum_quality=None, max_growth=None, stagnation_window=3, minimum_reduction=0.01, anisotropic=False)

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

``AntennaSolver3D.solve``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.solve(*, backend='auto', max_refinements=2, adaptive_tolerance=0.05, max_elements=None, max_dofs=None, pml_tolerance=None, pml_observable=None, max_pml_refinements=3)

Adapt against all independent port excitations, then report their coherent superposition.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``backend``
     - ``str``
     - Optional
     - ``'auto'``
     - SciPy CPU or NVIDIA cuDSS (gpu or hybrid). Auto probes cuDSS.
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
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Adaptive mesh element budget.
   * - ``max_dofs``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Maximum field-space size on an adaptive candidate mesh before eliminating boundary DOFs.
   * - ``pml_tolerance``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Relative tolerance on the observable change after doubling shell layers. Defaults to adaptive_tolerance.
   * - ``pml_observable``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional function returning a finite scalar or array from a result for the shell convergence check.
   * - ``max_pml_refinements``
     - ``int``
     - Optional
     - ``3``
     - Maximum shell doublings within the overall refinement budget.

Returns: an AntennaResult.

``AntennaSolver3D.sweep``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.sweep(*, frequencies, backend='auto', max_refinements=2, adaptive_tolerance=0.05, max_elements=None, max_dofs=None, pml_tolerance=None, pml_observable=None, max_pml_refinements=3)

Adapt at the highest frequency in hertz and reuse that mesh at every lower frequency.

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
   * - ``backend``
     - ``str``
     - Optional
     - ``'auto'``
     - SciPy CPU or NVIDIA cuDSS (gpu or hybrid). Auto probes cuDSS.
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
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Adaptive mesh element budget.
   * - ``max_dofs``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Maximum field-space size on an adaptive candidate mesh before eliminating boundary DOFs.
   * - ``pml_tolerance``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Relative tolerance on the observable change after doubling shell layers. Defaults to adaptive_tolerance.
   * - ``pml_observable``
     - ``array-like or scalar``
     - Optional
     - ``None``
     - Optional function returning a finite scalar or array from a result for the shell convergence check.
   * - ``max_pml_refinements``
     - ``int``
     - Optional
     - ``3``
     - Maximum shell doublings within the overall refinement budget.

Returns: an AntennaSweep.

``AntennaSolver3D.show``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSolver3D.show(*, block: 'bool' = True)

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

Returns: the displayed figure or plotter.

``AntennaResult.far_field``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaResult.far_field(directions)

Coherent far electric amplitude (V), intensity (W/sr), and linear power ratios.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``directions``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Cartesian observation directions, one per row.

Returns: a RadiationPattern.

``AntennaResult.show_port_modes``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaResult.show_port_modes(*, port, block=True)

Inspect the numerical port modes in the 2D waveguide viewer.

.. list-table:: Arguments
   :header-rows: 1
   :widths: 16 20 12 16 36

   * - Argument
     - Type / units
     - Required / optional
     - Default
     - Meaning
   * - ``port``
     - ``array-like or scalar``
     - Required
     - ``—``
     - Name of the configured waveguide port.
   * - ``block``
     - ``bool``
     - Optional
     - ``True``
     - Wait for the interactive viewer to close when true.

Returns: the 2D waveguide mode viewer.

``AntennaResult.plot``
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaResult.plot(*, component='S-parameters', quantity=None, phase_deg=0.0)

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
     - ``'S-parameters'``
     - S-parameters, Total E, Total H, Surface current, Directivity, Gain or Realized gain.
   * - ``quantity``
     - ``str``
     - Optional
     - ``None``
     - Network magnitude (dB) or phase; real or imaginary vector field.
   * - ``phase_deg``
     - ``float``
     - Optional
     - ``0.0``
     - Phase in degrees for the drive or displayed vector phasor.

Returns: a Matplotlib figure or VTK plotter.

``AntennaResult.show``
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaResult.show(*, block: 'bool' = True)

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

Returns: the displayed figure or plotter.

``AntennaResult.save``
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaResult.save(path)

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
     - HDF5 result path.

Returns: the HDF5 path.

``AntennaSweep.plot``
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSweep.plot(*, component='S-parameters', case=0, quantity=None, phase_deg=0.0)

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
     - ``'S-parameters'``
     - S-parameters, Total E, Total H, Surface current, Directivity, Gain or Realized gain.
   * - ``case``
     - ``int``
     - Optional
     - ``0``
     - Zero-based sweep case index.
   * - ``quantity``
     - ``str``
     - Optional
     - ``None``
     - Network magnitude (dB) or phase; real or imaginary vector field.
   * - ``phase_deg``
     - ``float``
     - Optional
     - ``0.0``
     - Phase in degrees for the drive or displayed vector phasor.

Returns: a Matplotlib figure or VTK plotter.

``AntennaSweep.show``
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSweep.show(*, block: 'bool' = True)

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

Returns: the displayed figure or plotter.

``AntennaSweep.save``
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    AntennaSweep.save(path)

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
     - HDF5 result path.

Returns: the HDF5 path.
