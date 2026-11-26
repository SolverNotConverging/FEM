fem.scattering user API
=======================

Version 2.0.2. This reference covers the deliberately supported user API.
Implementation helpers are documented in their source modules.

Plane-wave scattering
---------------------

See the `solver guide <../solvers/scattering.rst>`_ for field conventions and examples.

``ScatteringResult``, ``ScatteringSolver2D``, ``ScatteringSolver3D``.

``ScatteringSolver2D``
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver2D(*, frequency, x_range, y_range, polarization='TM', direction=(1, 0), pml_thickness=None, pml_layers=5)

Plane-wave scattering in the xy plane by infinite z-invariant targets.

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
   * - ``polarization``
     - ``str``
     - Optional
     - ``'TM'``
     - TE (Hz) or TM (Ez) relative to the invariant z axis.
   * - ``direction``
     - ``tuple``
     - Optional
     - ``(1, 0)``
     - Cartesian propagation direction, normalized internally.
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

``ScatteringSolver2D.add_geometry``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver2D.add_geometry(*, shape, material)

Add a bounded inclusion. Later objects take precedence in overlaps.

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

``ScatteringSolver2D.add_mesh_refinement``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver2D.add_mesh_refinement(*, shape, max_element_size)

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

``ScatteringSolver2D.add_pml``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver2D.add_pml(*, thickness: "float | Literal['auto']" = 'auto', order: 'int' = 3, sigma_max: 'float' = 5.0, direction: 'str' = 'all') -> 'PMLSpec'

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

``ScatteringSolver2D.add_nf2ff``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver2D.add_nf2ff(distance='auto', *, x_range=None, y_range=None, z_range=None)

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

``ScatteringSolver2D.mesh``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver2D.mesh(*, max_element_size=None, element_order=2, quadrature_order=5)

Build a conforming triangular or tetrahedral mesh with a target size in metres.

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

``ScatteringSolver2D.solve``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver2D.solve(*, backend='auto', max_refinements=2, adaptive_tolerance=0.05, max_elements=None, max_dofs=None)

Adapt the physical mesh using the relative volume, interface and boundary residual.

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

Returns: a ScatteringResult.

``ScatteringSolver2D.show``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver2D.show(*, block: 'bool' = True)

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

``ScatteringSolver3D``
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver3D(*, frequency, x_range, y_range, z_range, direction=(1, 0, 0), electric_vector=(0, 1, 0), pml_thickness=None, pml_layers=5)

Full-vector Maxwell scattering by bounded three-dimensional targets.

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
   * - ``direction``
     - ``tuple``
     - Optional
     - ``(1, 0, 0)``
     - Cartesian propagation direction, normalized internally.
   * - ``electric_vector``
     - ``tuple``
     - Optional
     - ``(0, 1, 0)``
     - Real incident electric polarization perpendicular to propagation.
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

``ScatteringSolver3D.add_geometry``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver3D.add_geometry(*, shape, material)

Add a bounded inclusion. Later objects take precedence in overlaps.

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

``ScatteringSolver3D.add_mesh_refinement``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver3D.add_mesh_refinement(*, shape, max_element_size)

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

``ScatteringSolver3D.add_pml``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver3D.add_pml(*, thickness: "float | Literal['auto']" = 'auto', order: 'int' = 3, sigma_max: 'float' = 5.0, direction: 'str' = 'all') -> 'PMLSpec'

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

``ScatteringSolver3D.add_nf2ff``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver3D.add_nf2ff(distance='auto', *, x_range=None, y_range=None, z_range=None)

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

``ScatteringSolver3D.mesh``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver3D.mesh(*, max_element_size=None, element_order=2, quadrature_order=5)

Build a conforming triangular or tetrahedral mesh with a target size in metres.

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

``ScatteringSolver3D.solve``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver3D.solve(*, backend='auto', max_refinements=2, adaptive_tolerance=0.05, max_elements=None, max_dofs=None)

Adapt the physical mesh using the relative volume, interface and boundary residual.

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

Returns: a ScatteringResult.

``ScatteringSolver3D.show``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringSolver3D.show(*, block: 'bool' = True)

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

``ScatteringResult.radar_cross_section``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringResult.radar_cross_section(directions)

Bistatic width (m) or area (m²) for each Cartesian observation direction.

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

Returns: an array of widths (m) or areas (m²).

``ScatteringResult.plot``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringResult.plot(*, component=None, quantity='real')

Plot a field such as 'Total Et' (TE) or 'Incident Ez' (TM), or 'RCS'.

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
     - Incident, total or scattered field. Choose Et/Hz for 2D TE, Ez/Ht for 2D TM, or E/H in 3D. Full 2D arrows include 'Total E vector' and 'Total H vector'. RCS in either dimension, or Surface current on a 3D conductor.
   * - ``quantity``
     - ``str``
     - Optional
     - ``'real'``
     - Real part, imaginary part, or magnitude.

Returns: a PyVista plotter for 3D or full vector views and a Matplotlib figure for planar fields and scattering width.

``ScatteringResult.show``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringResult.show(*, block: 'bool' = True)

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

``ScatteringResult.export_vtk``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringResult.export_vtk(path)

Write a .vtm collection with its field and angular RCS blocks for ParaView.

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
     - HDF5 result path for save, or .vtm path for export_vtk.

Returns: the VTM path.

``ScatteringResult.vtk_datasets``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringResult.vtk_datasets()

Return volume field vectors, a radial RCS surface and any conductor currents.

Returns: a VTK multiblock dataset.

``ScatteringResult.save``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    ScatteringResult.save(path)

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
     - HDF5 result path for save, or .vtm path for export_vtk.

Returns: the HDF5 path.
