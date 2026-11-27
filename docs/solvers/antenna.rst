Waveguide-fed antennas
======================

``fem.AntennaSolver3D`` solves the driven curl-curl Maxwell system with
tetrahedral edge elements. Enclose the complete antenna and its feeds in
a vacuum box. A Cartesian perfectly matched layer (PML) surrounds all six
faces. Bulk regions accept passive isotropic epsilon and mu. PEC and
surface-impedance bodies use boundary conditions on the conductor surface.

The x, y and z ranges include the PML. Its inner faces bound the physical
domain and its outer wall is PEC. The polynomial
coordinate stretch acts along each exterior axis, including both axes at
an edge and all three at a corner. The ports and the closed NF2FF surface
stay inside the physical box. Saved electric fields, magnetic fields and surface currents
cover the physical domain.

By default, each side has five mesh layers spanning one third of the
free-space wavelength at the constructor frequency. Use ``add_pml(thickness="auto")``
for this setting. Supply a thickness in metres, polynomial ``order`` (default 3),
peak ``sigma_max`` (default 5), and ``direction`` to customize a profile.
The first explicit call selects the PML sides, with PEC on the others.
Subsequent calls can assign separate profiles to disjoint sides.
Set ``pml_layers`` at construction to change the mesh resolution of the absorber.
A sweep keeps the mesh and stretch profile fixed. Actual reflections depend on incidence angle and mesh
resolution. For wide sweeps, check thickness at the lowest frequency and
layer resolution at the highest. Increase both thickness and layer count
and compare S-parameters and far fields to check convergence.

For finer physical elements near a feed or aperture, call
``add_mesh_refinement(shape=region, max_element_size=h)`` before ``mesh()``.
A box or sphere specifies the refined region in metres. The wave port
eigensolver uses the refined face triangles wherever the region intersects
its reference plane.

Numerical waveguide ports
-------------------------

Define each port on a flat, mesh-conforming conductor end face of a uniform
guide. The normal axis can be x, y or z. The other two coordinates in
increasing Cartesian order define the aperture in metres. Select it with
a 2D shape such as a polygon, rectangle, circle or annulus. Its edges,
internal conductor holes and dielectric interfaces come from the 3D CAD.

The full-vector 2D waveguide eigensolver uses the face triangles and the
adjacent material tags. The port edges inherit the PEC or surface-impedance condition of the
adjoining guide walls. Tangential electric and magnetic eigenfields define the multimode
Dirichlet-to-Neumann condition. Increase the retained basis to check
convergence near discontinuities. For a configured through-guide::

    aperture = fem.shapes.Rectangle(bounds=((-0.010, 0.010), (-0.004, 0.004)))
    solver.add_waveguide_port(name="input", axis="z", position=-0.015,
                              cross_section=aperture, num_modes=3, mode=1)
    modes = solver.solve_port_modes(port="input")
    print(modes.neff)
    modes.show()

Port ranks start at 1. Passing only ``mode=2`` retains two modes and drives
the second. Supply ``num_modes`` to retain additional modes in the termination.
Propagating modes are ordered by decreasing real
effective index, followed by evanescent modes in order of increasing
attenuation. The selected mode must propagate throughout a sweep. Ranks
are recomputed at each frequency. At crossings and degeneracies, inspect
the fields to identify the polarization. ModeSet sequence indexing in
Python starts at zero.

Choose a port and **View port modes** in the result viewer, or call
``case.show_port_modes(port="input")`` for the same 2D waveguide viewer.

Discrete gap feeds
------------------

Specify a box or cylinder of vacuum or bulk material. Both end faces
normal to the port axis must contact conductors across their entire area.
For a 50-ohm round-wire feed with a 3 mm radius and a 4 mm gap::

    gap = fem.shapes.Cylinder(center=(0, 0), radius=0.003, z_range=(-0.002, 0.002))
    solver.add_discrete_port(name="feed", axis="z", gap=gap, impedance=50.)

Rotate the cylinder to align its length with an x- or y-directed feed.
Voltage is the area-average of electric-field line integrals along the
positive axis. The reciprocal current distribution is uniform across the
gap. Use a gap small compared with the local wavelength. ``case.port_voltage``
and ``case.port_current`` contain peak terminal phasors in V and A. Accepted
terminal power is Re(V I*)/2. Discrete and waveguide ports share a drive
vector and may be combined in the same frequency sweep.

Network and coherent drive
--------------------------

Amplitudes are incident power waves in square roots of watts and phases
are in degrees, with the exp(+iωt) convention. For example,
``solver.set_excitation(port="feed", amplitude=.5, phase_deg=90.)`` sets
a 0.25 W incident wave with phase +90°. Zero amplitude supplies a passive
terminated port. Each waveguide port has one selected launch mode.

Each solve factors the system once and solves a unit-drive column for each
selected port mode. ``result.channels`` gives the order of (port, mode) pairs
in ``result.s_parameters``. Rows receive and columns launch.
``result.modal_s_parameters[(port, mode)]`` gives the outgoing row for any
retained propagating mode. The columns are independent of drive amplitude
and phase. The combined fields use the complex incident vector a.
``result.active_s_parameters`` gives (Sa)ᵢ/aᵢ at each selected port mode.
It is undefined (NaN) at an undriven port. Its outgoing wave remains available
in ``result.outgoing``.

Radiation patterns and frequency sweeps
---------------------------------------

The solver meshes a closed Huygens box in vacuum between the geometry and
outer boundary. It samples both equivalent currents, n × H and −n × E,
from the solved fields. Call ``add_nf2ff(distance="auto")`` for a target inset
of 0.3 free-space wavelengths from each PML interface, capped at half the
vacuum clearance to the nearest source or material on that side. Automatic
PML and NF2FF settings apply when the methods are omitted. To place the
surface explicitly, pass a positive distance in metres or call
``add_nf2ff(x_range=..., y_range=..., z_range=...)`` before meshing.
It must enclose all sources and non-vacuum materials with vacuum on both sides.
The mesh stores its resolved coordinates in ``mesh_data.info["nf2ff_bounds"]``.
Refine the mesh, move the Huygens
surface, and expand the radiation box separately to check convergence.

``result.far_field(directions)`` returns the complex Cartesian coefficient
of exp(−ikr)/r in volts, intensity U in W/sr, and three linear ratios.
Directivity is 4πU divided by integrated radiated power. Gain uses accepted
power, computed from net modal flux at wave ports and terminal power at
discrete ports. Realized gain uses total incident power. The three
powers are available in watts on the result. Angular integration uses
Gauss–Legendre quadrature in cos(θ) and a periodic azimuthal grid. For a
lossless radiating antenna, compare integrated radiated power with accepted
power as an independent numerical check.

Use the case selector for frequency and the component selector for S-parameters,
frequency response, total E, total H, surface current, directivity, gain or realized gain. Network plots
offer magnitude in dB and phase in degrees. Radiation surfaces have a linear
radius and dBi colour scale. Total H shows the magnetic field recovered from
the curl of E and the material permeability in each tetrahedron, in A/m.
Surface-current arrows show n × H on PEC and impedance walls in A/m, with n pointing out of the conductor. For each
vector field, set **Phase** in degrees or select **Animate phase** to cycle
through one RF period. Arrow direction, length and colour follow the field
at the displayed phase. The length and colour scales stay fixed during playback.
For example::

    sweep = solver.sweep(frequencies=[9.5e9, 10e9, 10.5e9], backend="cpu")
    print(sweep.s_parameters[:, 0, 0])
    print(sweep[1].far_field([[0,0,1]]).realized_gain)
    sweep.save("antenna.h5")
    fem.load("antenna.h5").show()

The sweep adapts at its highest frequency and reuses that final mesh for
each lower frequency. Port modes and surface impedances are recomputed at
each frequency. The drive weights apply throughout. Both ``solve()`` and
``sweep()`` allow eight updates after the initial solve. Physical convergence
requires at least three physical meshes and a maximum absolute change below
0.02 across all complex port-matrix entries. An explicit ``adaptive_tolerance``
selects the physical field-residual test unless a port convergence rule was
explicitly configured. Local residuals across all independent excitations
select cells for refinement. For an excited PML, a second check compares the
complex S-matrix while doubling shell layers on the same physical mesh. Set
``pml_tolerance`` or pass ``pml_observable`` to monitor a radiation observable.
The highest-frequency case records the refinement
history. Each case records ``mesh_adaptation_frequency`` in ``solve_info``.
The volume defaults to second-order Nédélec tetrahedra with two degrees of
freedom per edge and two per face. The port solves its modes with the same
order on the triangular cross-section. Set ``mesh(element_order=1)`` to use
first order for both. The quadratic tetrahedron has 20 local functions
in the `first-kind Nédélec space <https://defelement.org/elements/nedelec1.html>`_.
Integration defaults to degree five. The viewer
uses cell-averaged E and H. For second order, ``surface_current_vertices``
holds three one-sided current vectors per conducting triangle in A/m.
The viewer samples their linear interpolant when displaying or animating
surface current. ``surface_current`` gives the face mean at either order.
The sweep progress bar shows the frequency being solved. Adaptive diagnostics
include the pass number, physical element count, residual, tolerance and
stopping reason.
Automatic quadratic meshing selects a profile from the port apertures,
conductor geometry and bulk materials. Rectangular through-guides use
length- and filling-dependent targets. Printed conductors get finer cells
around the substrate, and thin-wire antennas get finer cells around their
arms. The general target is four elements per local wavelength at order two
and ten at order one. The highest sweep frequency sets the initial sizes.
See the `driven mesh profiles <../adaptivity.rst>`_ for the numerical settings.
Use ``mesh(wavelength_elements=8)`` to change the wavelength target while
retaining material-based sizing. An optional ``cad_point_size`` sets a length
in metres at CAD vertices. A shorter adjacent material wavelength or discrete
port gap can tighten that target. For a port-driven convergence study::

    solver.mesh(wavelength_elements=4, element_order=2)
    solver.configure_adaptivity(anisotropic=True)
    solver.configure_s_parameter_convergence(tolerance=0.02, minimum_passes=3,
                                            consecutive_checks=1)
    result = solver.solve(max_refinements=7, backend="cpu")

This compares the largest absolute change of any complex S-matrix entry at the
adaptation frequency between successive physical meshes. The initial solution
counts as pass one. A change at or above the threshold resets the consecutive
check count. The local field residual still selects cells for refinement.
For an excited absorber, shell-only checks follow physical convergence and use
the remaining refinement budget. Pass ``tolerance=None`` to restore the
residual stopping rule.

SciPy and NVIDIA cuDSS are available through ``backend="cpu"``
and ``backend="gpu"``. ``auto`` probes cuDSS and can fall back to CPU.

The `antenna examples <../../examples/antenna/README.rst>`_ cover a through-guide,
inset-fed microstrip patch, sectoral horn and gap-fed dipole.
