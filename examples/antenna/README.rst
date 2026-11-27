Antenna examples
================

Run a script from the repository root. Each sweeps five frequencies and opens
the shared result viewer. Select a frequency to inspect the fields and radiation
pattern on the mesh adapted at the highest frequency.
The microstrip patch example allows two refinements after the initial solve.
The other examples use the default budget of eight.
At the highest frequency, physical adaptation stops once the largest change
in a complex S-parameter falls below 0.02 on or after the third physical mesh.
The solver also checks convergence of the absorbing layer within the same budget.
Lower frequencies reuse the resulting mesh. All four examples surround the
physical box with a five-layer PML.

* `Through waveguide <through_waveguide.py>`_ sweeps a dielectric-filled
  two-port line from 7 to 9 GHz and saves the network and a frequency-response plot.
* `Microstrip patch <microstrip_patch.py>`_ launches the fundamental microstrip
  mode at the board edge. A continuous copper strip feeds an inset
  in the copper patch above a finite copper ground plane. The FR-4 substrate
  uses relative permittivity 4.4 and loss tangent 0.02. The example reports
  the port effective index, input reflection and broadside realized gain from 2.35 to 2.55 GHz.
* `Sectoral horn <sectoral_horn.py>`_ launches a mode in a straight
  throat and reports boresight realized gain from 9 to 11 GHz.
* `Dipole <dipole.py>`_ excites two round PEC arms through a cylindrical
  50-ohm gap port over a band of approximately 270 to 330 MHz. Select **Surface current** and **Animate phase** to view
  the RF current along both arms.

The `solver guide <../../docs/solvers/antenna.rst>`_ describes port rank,
drive normalization, phase controls and the near-field to far-field transform.
