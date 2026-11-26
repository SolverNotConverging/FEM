Plane-wave scattering examples
==============================

Run either script from the repository root after installing FEM::

    python examples/scattering/2d_dielectric_cylinder.py
    python examples/scattering/3d_copper_sphere.py

The `2D cylinder <2d_dielectric_cylinder.py>`_ example saves one Matplotlib
figure per polarization. Each figure contains the incident, total and
scattered field and the angular scattering width in metres.

The `3D sphere <3d_copper_sphere.py>`_ example applies the copper surface
impedance at 3 GHz. It opens a VTK viewer and writes a ParaView collection
containing complex electric fields on tetrahedra and bistatic RCS in m² on a
radial surface. Files are under ``outputs/<example>/`` beside each script.

Both examples use a five-layer exterior PML and automatic initial meshing.
Both examples pass ``max_refinements=2`` to ``solve()``. After the initial solve,
the solver refines up to twice, stopping early when the residual meets the tolerance.

See the `solver guide <../../docs/solvers/scattering.rst>`_ for incidence,
materials, mesh convergence and CUDA selection.
