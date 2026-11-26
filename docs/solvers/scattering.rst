Plane-wave electromagnetic scattering
=====================================

Choose ``ScatteringSolver2D`` for targets invariant along z, or
``ScatteringSolver3D`` for a finite target. Both solve for a single plane
wave in vacuum. Lengths are in metres and frequency is in hertz.

2D cylinder
-----------

.. code-block:: python

    import fem

    solver = fem.ScatteringSolver2D(frequency=3e9, x_range=(-.163, .163),
                                    y_range=(-.163, .163), polarization="TE",
                                    direction=(1, 0))
    rod = fem.shapes.Circle(center=(0, 0), radius=.02)
    solver.add_geometry(shape=rod, material=fem.Material(epsilon=2.25))
    solver.add_pml(thickness="auto")
    solver.add_nf2ff(distance="auto")
    result = solver.solve(backend="auto")
    result.show()
    print(result.radar_cross_section([[-1, 0], [1, 0]]))  # back and forward width (m)

TE means transverse electric to z and solves for Hz (A/m). TM solves for Ez
(V/m). The incident scalar amplitude is one. The propagation direction is a
Cartesian vector in the xy plane. Phasors use exp(+iωt) and the incident
spatial phase is exp(−ik d·x), with the phase origin at the coordinate origin.
The complex nodal arrays ``incident``, ``scattered`` and ``total`` share the same
mesh. ``result.plot(quantity="real")`` returns a Matplotlib figure showing all
three fields and the polar scattering width. Use ``imag`` or ``magnitude`` for
the other field views. For each of the incident, total and scattered fields,
the viewer offers **Et** and **Hz** for TE, or **Ez** and **Ht** for TM.
Choose ``Total E vector`` or ``Total H vector`` for rotatable 2D field plots.
These show the transverse components in the mesh plane or the longitudinal
component perpendicular to it, according to polarization.
Et and Ht use arrows at triangle centres, coloured by the complex vector
magnitude. Real and imaginary views show the corresponding signed direction.
Ez and Hz use nodal heat maps. Select RCS for the polar scattering width.
For example, ``result.plot(component="Total Et", quantity="real")`` plots the
TE electric field. The complex x and y components at triangle centres are
available as ``incident_transverse``, ``total_transverse`` and ``scattered_transverse``.

3D target and surface impedance
-------------------------------

.. code-block:: python

    solver = fem.ScatteringSolver3D(frequency=3e9, x_range=(-.093, .093),
                                    y_range=(-.093, .093), z_range=(-.093, .093),
                                    direction=(1, 0, 0), electric_vector=(0, 1, 0))
    sphere = fem.shapes.Sphere(center=(0, 0, 0), radius=.015)
    solver.add_geometry(shape=sphere, material=fem.materials.copper)
    result = solver.solve(backend="gpu")
    print(result.radar_cross_section([[-1, 0, 0]]))  # monostatic RCS (m²)
    result.export_vtk("outputs/sphere.vtm")
    result.show()

The 3D solver uses Nédélec edge elements. The incident E vector
is linearly polarized, transverse to propagation, and normalized to 1 V/m.
The viewer offers E (V/m) and H (A/m) arrows for the incident, total and
scattered fields. For example, select **Total H** or call
``result.plot(component="Total H", quantity="real")``. Each complex array has
shape (tetrahedra, 3) and samples the field at physical cell centres.
``incident``, ``total`` and ``scattered`` hold E, and ``incident_h``,
``total_h`` and ``scattered_h`` hold H. The solver recovers total H from the
curl of total E using each cell's permeability. Real and imaginary views
orient the arrows by that phasor part. Magnitude uses a global phase
reference. Colour gives the full vector magnitude. Select RCS to display the radial surface. The mesh
and material toggles overlay the geometry in spatial field views. In ParaView,
open the VTM collection and select the fields or RCS block. Real and imaginary
Cartesian vectors and their magnitudes are available for each field, for
example ``total_real`` for E and ``total_h_real`` for H. Keep
the accompanying dataset directory with the VTM file.

Materials and accuracy
----------------------

Assign an isotropic ``fem.Material(epsilon=..., mu=...)``, ``fem.materials.PEC``,
a constant ``fem.SurfaceImpedance(impedance=...)`` in ohms, or a
frequency-dependent ``fem.GoodConductor`` such as ``fem.materials.copper``.
For surface impedance (SIBC) and PEC, the solver meshes the exterior and
applies the condition directly on the target surface. Bulk dielectric and
magnetic regions retain their interior fields. Later objects take precedence
in overlaps. All objects must fit inside the box with a vacuum margin.

The supplied box ranges include a Cartesian perfectly matched layer
(PML) on every side by default. The layer absorbs the scattered field with a polynomial
stretch in each exterior coordinate and a zero scalar or tangential electric
field at the outer wall. Field plots, exports and RCS sources cover the
physical domain. By default each side has five layers spanning one third of
the free-space wavelength. Use ``add_pml(thickness="auto", order=3, sigma_max=5.0,
direction="all")`` or assign profiles to individual faces such as ``"x+"``.
An explicit selection terminates uncovered sides with zero scattered scalar or
tangential electric field. Set ``pml_layers`` at construction to vary absorber
resolution and check RCS convergence. Independently check spatial
convergence by refining the physical mesh at fixed box size. Resolve the shortest
material wavelength and geometric features. To refine a target locally,
call ``add_mesh_refinement(shape=region, max_element_size=h)`` before
``mesh()``. Use a rectangle or circle in 2D, or a box or sphere in 3D.
Overlapping size controls use the smallest target.

``add_nf2ff(distance="auto")`` places a closed vacuum contour in 2D or a
surface in 3D. Each face sits 0.3 free-space wavelengths inside the PML
interface, capped at half the available clearance to the target on that side.
Pass a positive distance in metres for a uniform inset or supply explicit
``x_range``, ``y_range`` and, in 3D, ``z_range``. The surface must enclose all
targets with vacuum on both sides. Automatic PML and NF2FF settings also
apply when the methods are omitted. The mesh records the resolved monitor
coordinates in ``mesh_data.info["nf2ff_bounds"]``.

``radar_cross_section(directions)`` integrates the scattered field and its
derivative on the NF2FF contour or surface. In 2D it returns
the far-field limit of 2πr times the scattered-to-incident intensity ratio
(scattering width in m). In 3D the factor is 4πr² (area in m²). Observation
directions are Cartesian unit vectors pointing away from the target. The
negative incident direction gives backscatter. In the VTK RCS block, radius
and colour both give the linear cross section in m². The point-data array
``direction`` stores the observation vector, including at pattern nulls.

Execution and outputs
---------------------

The 2D scalar formulation defaults to quadratic triangles and the 3D vector
formulation defaults to second-order Nédélec tetrahedra. Both support
first order through ``mesh(element_order=1)`` and use degree-five
quadrature. Set the orders with ``mesh(element_order=..., quadrature_order=...)``.
The tetrahedral geometry is affine at either field order. Refine the surface
mesh to resolve the curvature of a smooth target.
By default ``solve(max_refinements=2, adaptive_tolerance=0.05)`` checks
volume equations, interface flux jumps and natural boundary conditions after
each solve. It refines the cells with the largest contributions. A zero budget uses the current mesh. The final
residual and stopping reason are available in ``result.solve_info``.
Automatic meshing targets six elements per local wavelength at order two and
ten at order one. The dielectric sizing applies inside the target so that
vacuum and its exterior PML can use the free-space wavelength.

For a 3D PEC or impedance target, select **Surface current** to inspect
the total-field conductor current. The VTK ``surface_current`` block
stores real and imaginary vectors in A/m at each triangle's corners.
Separate corner values preserve the two one-sided traces at a sharp edge.

``backend="cpu"`` uses SciPy. ``backend="gpu"`` and ``backend="hybrid"`` use
NVIDIA cuDSS for the sparse complex direct solve. ``backend="auto"`` probes
cuDSS and falls back to SciPy with a warning if unavailable. Assembly and
far-field evaluation run on the CPU. See `CUDA setup <../acceleration.md>`_.

``result.save("case.h5")`` stores the mesh, complex fields and equivalent
sources. ``fem.load("case.h5")`` supports visualization and new observation
directions. See the `examples <../../examples/scattering/README.rst>`_ and
`API reference <../api/scattering.rst>`_.
