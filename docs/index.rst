fem
===

Finite-element electromagnetic solvers with one workflow for every problem type:

.. code-block:: python

   import fem

   fr4 = fem.Material(name="FR-4", epsilon=4.4 * (1 - 0.02j))  # lossy dielectric, tan(delta) = 0.02
   copper = fem.materials.copper                                # surface-impedance conductor

   solver = fem.ModeSolver2D(frequency=10e9, x_range=(-6e-3, 6e-3), y_range=(-35e-6, 6e-3),
                             boundary=fem.materials.PEC)
   solver.add_rectangle(x_range=(-6e-3, 6e-3), y_range=(0, 1.6e-3), name="substrate", material=fr4)
   solver.add_rectangle(x_range=(-6e-3, 6e-3), y_range=(-35e-6, 0), name="ground", material=copper)
   solver.add_rectangle(x_range=(-1.5e-3, 1.5e-3), y_range=(1.6e-3, 1.635e-3), name="strip", material=copper)
   solver.mesh(max_element_size=0.5e-3)
   modes = solver.solve(num_modes=1)

   print(f"neff = {modes[0].neff:.4f}, loss = {modes[0].alpha * 8.686:.1f} dB/m")
   modes.save("microstrip.h5")
   modes.show()  # fields, mesh and material regions

Selecting a field or toggling **Mesh** or **Material** preserves pan and zoom
in planar plots and the camera pose, zoom and projection in 3D views. Each
viewer remembers the framing of its planar and rotatable views when switching
between them. Use **Reset view** at the bottom to fit the current selection.

.. list-table:: Solvers
   :header-rows: 1
   :widths: 35 65

   * - Class
     - Problem
   * - ``fem.ElectrostaticSolver``
     - Potentials, fields, charges and energies in 1D, 2D or 3D
   * - ``fem.ModeSolver1D``, ``fem.ModeSolver2D``
     - Guided modes of a z-invariant slab or cross-section
   * - ``fem.PeriodicModeSolver2D``, ``fem.PeriodicModeSolver3D``
     - Complex Bloch modes of a structure periodic along z (leaky-wave antennas)
   * - ``fem.ScatteringSolver2D``, ``fem.ScatteringSolver3D``
     - Plane-wave scattering fields and bistatic RCS for 2D TE/TM and 3D targets
   * - ``fem.AntennaSolver3D``
     - Numerical waveguide and discrete feeds, active S-parameters and radiation patterns

.. toctree::
   :maxdepth: 2
   :caption: Guide

   getting_started
   conventions
   materials_and_shapes
   adaptivity
   acceleration
   benchmarks/README

.. toctree::
   :maxdepth: 2
   :caption: Solvers

   solvers/electrostatics
   solvers/waveguide
   solvers/periodic
   solvers/image_line_pilot
   solvers/scattering
   solvers/antenna

.. toctree::
   :maxdepth: 1
   :caption: Reference

   api/common
   api/electrostatics
   api/waveguide
   api/periodic
   api/scattering
   api/antenna
   api/adaptivity
   changelog
