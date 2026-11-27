"""Scattered-field Helmholtz and curl-curl finite elements in vacuum."""
import numpy as np
from dataclasses import replace
from scipy.sparse.linalg import spsolve
from skfem import Basis, BilinearForm, ElementTriP1, ElementTriP2, FacetBasis, LinearForm, asm
from skfem.helpers import curl, dot, grad

from fem.common import materials, shapes
from fem.common.contracts import bounds
from fem.common.constants import ETA_0
from fem.common.discretization import DEFAULT_QUADRATURE_ORDER, field_orders
from fem.adaptivity import run_adaptive
from fem.adaptivity.driven import AbsorberCheck, driven_residual
from fem.common.errors import ConfigurationError, GeometryError, SolverError
from fem.common.linalg import resolve_backend
from fem.common.materials import bulk_values
from fem.common.nf2ff import interface_facets
from fem.common.radiation_pml import shell_matrix
from fem.common.tet_hcurl import conductor_current, tetrahedral_element
from fem.common.validation import positive_real, real_array

from .meshing import RadiationBoxSolver, build_mesh
from .results import ScatteringResult


def _unit(vector, dimension, name):
    vector = real_array(vector, name)
    if vector.shape != (dimension,) or np.linalg.norm(vector) == 0:
        raise ConfigurationError(f"{name} must be a nonzero {dimension}-vector.")
    return vector / np.linalg.norm(vector)


def _monitor_sources(basis, box, coefficients, k, quadrature_order):
    """Recover outward flux by testing the exterior vacuum cells against the surface DOFs."""
    mesh, element = basis.mesh, basis.elem
    vector = mesh.dim() == 3
    derivative = curl if vector else grad
    product = dot if vector else np.multiply
    facets, inside_first = interface_facets(mesh, box)
    dofs = basis.get_dofs(facets=facets).all()
    physical = mesh.t[:, basis.tind]
    inside = box.contains(*mesh.p[:, physical].mean(axis=1))
    touches = np.isin(physical, mesh.facets[:, facets]).any(axis=0)
    exterior = basis.with_elements(basis.tind[touches & ~inside])
    field = exterior.interpolate(coefficients)

    @LinearForm(dtype=complex)
    def reaction(v, w):
        return dot(derivative(field), derivative(v)) - k**2 * product(field, v)

    # The exterior weak residual on the monitor gives n × curl(E) in 3D
    # or -du/dn in 2D, where n points out of the enclosed region.
    rhs = asm(reaction, exterior)[dofs]
    face = FacetBasis(mesh, element, facets=facets, intorder=quadrature_order)

    @BilinearForm
    def trace_mass(u, v, w):
        return dot(u, v) - dot(u, w.n) * dot(v, w.n) if vector else u * v

    mass = asm(trace_mass, face)[dofs][:, dofs]
    density = np.zeros(basis.N, dtype=complex)
    density[dofs] = spsolve(mass.tocsc(), rhs)
    # FacetBasis takes its normal from side zero. Orient it out of the monitor.
    normal = face.normals * np.where(inside_first, 1., -1.)[None, :, None]
    trace = face.interpolate(coefficients)
    flux = face.interpolate(density)
    if vector:
        current = flux - normal * dot(normal, flux)
        moment = -np.cross(normal, trace, axis=0)
    else:
        current = flux[None]
        moment = -normal * trace
    points = np.moveaxis(face.global_coordinates(), 0, -1).reshape(-1, mesh.dim())
    return (points, (current * face.dx).reshape(len(current), -1).T,
            (moment * face.dx).reshape(len(moment), -1).T)


class _PlaneWaveSolver(RadiationBoxSolver):
    def _configure(self, frequency, ranges, direction, pml_thickness, pml_layers):
        self.frequency = positive_real(frequency, "frequency")
        self._configure_pml(pml_thickness, pml_layers)
        self._ranges = tuple(bounds(v, f"{axis}_range") for axis, v in zip("xyz", ranges))
        self.direction = _unit(direction, len(ranges), "direction")
        self._regions = []
        self._refinements = []
        self.mesh_data = self.result = self._mesh = self._size = None

    def add_geometry(self, *, shape, material):
        """Add a bounded inclusion. Later objects take precedence in overlaps."""
        dim = len(self._ranges)
        if not isinstance(shape, shapes.Shape) or shape.dimension != dim:
            raise GeometryError(f"Expected a {dim}D shape.")
        self._check_geometry_clearance(shape.bounds, vacuum=self._vacuum(material))
        if material != materials.PEC and not isinstance(material, (materials.GoodConductor,
                                                                    materials.SurfaceImpedance)):
            eps, mu = bulk_values(material, form="scalar")
            if eps == 0 or mu == 0:
                raise ConfigurationError("Permittivity and permeability must be nonzero.")
        self._regions.append((shape, material))
        self.mesh_data = self.result = self._mesh = None

    def mesh(self, *, max_element_size=None, element_order=2, quadrature_order=DEFAULT_QUADRATURE_ORDER):
        """Build a conforming triangular or tetrahedral mesh with a target size in metres."""
        element_order, quadrature_order = field_orders(
            self._default_element_order if element_order is None else element_order,
            quadrature_order, (1, 2))
        policy, size, point, local = self._driven_mesh_settings(self.wavelength, element_order, max_element_size)
        self._surface = self._nf2ff_surface()
        material_aware = max_element_size is None or (self._adaptive_sizes is not None and self._mesh_auto_materials)
        physical, snapshot, faces = build_mesh(
            self._physical_ranges(self._pml_profiles), self._regions, size,
            interior_box=self._surface, refinements=(*self._refinements, *local), adaptive_sizes=self._adaptive_sizes,
            material_aware=material_aware, cad_point_size=point)
        snapshot = replace(snapshot, info={**snapshot.info, "mesh_profile": policy.name,
                            "automatic_local_sizes": [h for _, h in local]})
        self._attach_pml(physical, snapshot, faces, element_order, quadrature_order)
        self._apply_driven_policy(policy)
        self._mesh_auto_materials = material_aware
        self._size, self.result = size, None
        return self.mesh_data

    def solve(self, *, backend="auto", max_refinements=8, adaptive_tolerance=0.02,
              max_elements=None, max_dofs=None, pml_tolerance=None, pml_observable=None,
              max_pml_refinements=3):
        """Refine the physical mesh by residual and the PML by observable change."""
        selected = resolve_backend(backend)
        absorber = AbsorberCheck(adaptive_tolerance if pml_tolerance is None else pml_tolerance,
                                 pml_observable, max_pml_refinements)
        return run_adaptive(self, lambda: self._solve_once(selected), lambda result: self._estimate,
                            self._refine_driven, max_refinements=max_refinements, max_elements=max_elements, max_dofs=max_dofs,
                            adaptive_tolerance=adaptive_tolerance, error_type=ConfigurationError, absorber=absorber)

    def _solve_once(self, selected):
        if self._mesh is None:
            self.mesh(max_element_size=None if self._mesh_auto_materials else self._size, **self._mesh_orders)
        mesh = self._mesh
        vector = mesh.dim() == 3
        quadrature = self._mesh_orders["quadrature_order"]
        order = self._mesh_orders["element_order"]
        element = tetrahedral_element(order) if vector else (ElementTriP1() if order == 1 else ElementTriP2())
        k = 2 * np.pi / self.wavelength
        tags = self.mesh_data.metadata["element_tags"]
        basis = Basis(mesh, element, elements=np.arange(len(tags)), intorder=quadrature)
        props = [(1., 1.)] + [bulk_values(m, form="scalar") if isinstance(m, materials.Material) else (1., 1.)
                              for _, m in self._regions]
        eps, mu = np.asarray(props, dtype=complex)[tags].T
        a, b = (1 / eps, mu) if self.polarization == "TE" else (1 / mu, eps)
        a, b = a[:, None], b[:, None]
        derivative = curl if vector else grad
        product = dot if vector else np.multiply

        def phase(x):
            return np.exp(-1j * k * np.einsum("i,i...->...", self.direction, x))

        def incident(x):
            wave = phase(x)
            return self.electric_vector[:, None, None] * wave if vector else wave

        def incident_derivative(x):
            axis = np.cross(self.direction, self.electric_vector) if vector else self.direction
            return -1j * k * axis[:, None, None] * phase(x)

        @BilinearForm(dtype=complex)
        def interior(u, v, w):
            return a * dot(derivative(u), derivative(v)) - k**2 * b * product(u, v)

        @LinearForm(dtype=complex)
        def source(v, w):
            return (k**2 * (b - 1) * product(incident(w.x), v)
                    - (a - 1) * dot(incident_derivative(w.x), derivative(v)))

        @BilinearForm(dtype=complex)
        def impedance_trace(u, v, w):
            trace = dot(u, v) - dot(u, w.n) * dot(v, w.n) if vector else u * v
            return 1j * k * trace

        matrix = asm(interior, basis) + shell_matrix(
            basis, len(tags), self._physical_ranges(self._pml_profiles), self._pml_profiles, k)
        rhs = asm(source, basis)
        coeff = np.zeros(basis.N, dtype=complex)
        prescribed = list(basis.get_dofs(facets=self._pml_wall).all())
        pec_traces = []
        natural_walls = []
        for owner, facets in self._faces.items():
            material = self._regions[owner - 1][1]
            metal = FacetBasis(mesh, element, facets=facets, intorder=quadrature)
            pec = material == materials.PEC
            admittance = 0.
            if not pec:
                z = material.at_frequency(frequency=self.frequency) / ETA_0
                admittance = z if self.polarization == "TE" else 1 / z
                matrix = matrix + admittance * asm(impedance_trace, metal)
            if not pec or self.polarization == "TE":
                natural_walls.append((facets, admittance))
                @LinearForm(dtype=complex)
                def metal_load(v, w):
                    u, du = incident(w.x), incident_derivative(w.x)
                    if vector:
                        trace = dot(u, v) - dot(u, w.n) * dot(v, w.n)
                        return dot(np.cross(w.n, du, axis=0), v) - 1j*k*admittance*trace
                    return (-dot(w.n, du) - 1j*k*admittance*u) * v
                rhs = rhs + asm(metal_load, metal)
            else:
                fixed = basis.get_dofs(facets=facets).all()
                prescribed.extend(fixed)
                if vector and order == 2:
                    pec_traces.append(facets)
                elif vector:
                    # Nedelec coefficients are line integrals along oriented mesh edges.
                    edges = mesh.edges[:, fixed]
                    delta = mesh.p[:, edges[1]] - mesh.p[:, edges[0]]
                    midpoint = mesh.p[:, edges].mean(axis=1)
                    coeff[fixed] = -(self.electric_vector @ delta) * phase(midpoint) * np.sinc(
                        k * (self.direction @ delta) / (2 * np.pi))
                else:
                    coeff[fixed] = -phase(basis.doflocs[:, fixed])
        if pec_traces:
            # Project the complete incident tangential trace onto the edge
            # and face space. Assemble all touching PEC patches together.
            facets = np.unique(np.concatenate(pec_traces))
            face = FacetBasis(mesh, element, facets=facets, intorder=quadrature)
            trace_dofs = basis.get_dofs(facets=facets).all()

            @BilinearForm
            def trace_mass(u, v, w):
                return dot(u, v) - dot(u, w.n) * dot(v, w.n)

            @LinearForm(dtype=complex)
            def trace_load(v, w):
                u = incident(w.x)
                return -dot(u, v) + dot(u, w.n) * dot(v, w.n)

            mass = asm(trace_mass, face)[trace_dofs][:, trace_dofs]
            coeff[trace_dofs] = spsolve(mass.tocsc(), asm(trace_load, face)[trace_dofs])
        free = basis.complement_dofs(prescribed)
        rhs = rhs - matrix @ coeff
        if np.any(rhs[free]):
            reduced = matrix[free][:, free]
            if selected == "cpu":
                coeff[free] = spsolve(reduced.tocsc(), rhs[free])
            else:
                from fem.common.gpu import solve_sparse
                coeff[free] = solve_sparse(reduced, rhs[free])
        if not np.isfinite(coeff).all():
            raise SolverError("The plane-wave linear solve returned nonfinite coefficients.")

        self._estimate = driven_residual(basis, coeff, b.ravel(), a.ravel(), k, quadrature,
                                         incident, incident_derivative, walls=natural_walls,
                                         absorber=(self._physical_ranges(self._pml_profiles), self._pml_profiles),
                                         marking_weights=self._marking_weights)
        positions, current, moment = _monitor_sources(basis, self._surface, coeff, k, quadrature)
        incident_transverse = scattered_transverse = None
        scattered_h = None
        if vector:
            # Sample both phasors at the glyph location. For a quadratic
            # vector field its centroid value differs from its cell mean.
            sample = Basis(mesh, element, elements=basis.tind, dofs=basis.dofs,
                           quadrature=(np.full((3, 1), .25), np.array([1 / 6])), disable_doflocs=True)
            field = sample.interpolate(coeff)
            scattered = np.asarray(field)[:, :, 0].T.copy()
            centers = self.mesh_data.coordinates[self.mesh_data.elements].mean(axis=1).T
            incoming = phase(centers)[:, None] * self.electric_vector
            # Faraday's law applies to total E with the local permeability.
            # The incident plane wave is defined in vacuum throughout the mesh.
            incident_h = np.cross(self.direction, incoming) / ETA_0
            scattered_h = (1j / (k * ETA_0) * field.curl[:, :, 0].T
                           + (1 - mu[:, None]) * incident_h) / mu[:, None]
        else:
            physical_points = self.mesh_data.coordinates
            scattered, incoming = coeff[:len(physical_points)].copy(), phase(physical_points.T)
            samples = Basis(mesh, element, elements=basis.tind, dofs=basis.dofs,
                            quadrature=(np.full((2, 1), 1/3), np.array([.5])), disable_doflocs=True)
            # For exp(+iωt), the transverse field follows from the curl of
            # the longitudinal field. Evaluate each triangle on its own side
            # of material interfaces, including the quadratic basis gradients.
            factor = (-1j * ETA_0 if self.polarization == "TE" else 1j / ETA_0) / k

            def transverse(gradient):
                dx, dy = gradient[:, :, 0]
                return factor * np.column_stack((dy, -dx))

            incident_transverse = transverse(incident_derivative(samples.global_coordinates()))
            scattered_transverse = (a * transverse(samples.interpolate(coeff).grad)
                                    + (a - 1) * incident_transverse)
        current_facets = current_vertices = None
        if vector and self._faces:
            conductors = np.unique(np.concatenate(list(self._faces.values())))
            current_facets = mesh.facets[:, conductors].T
            current_vertices = conductor_current(basis, coeff, mu, k, conductors, incident_derivative)
        self.result = ScatteringResult(self.frequency, self.polarization, self.direction.copy(),
                                       self.mesh_data, incoming, scattered, positions, current,
                                       moment, {"backend": selected, "ndofs": basis.N,
                                                "electric_vector": self.electric_vector.copy() if vector else None,
                                                "method": "scipy-direct" if selected == "cpu" else "cudss-direct"},
                                       current_facets=current_facets, surface_current_vertices=current_vertices,
                                       incident_transverse=incident_transverse, scattered_transverse=scattered_transverse,
                                       scattered_h=scattered_h)
        return self.result


class ScatteringSolver2D(_PlaneWaveSolver):
    """Plane-wave scattering in the xy plane by infinite z-invariant targets."""
    def __init__(self, *, frequency, x_range, y_range, polarization="TM", direction=(1, 0),
                 pml_thickness=None, pml_layers=5):
        self._configure(frequency, (x_range, y_range), direction, pml_thickness, pml_layers)
        if not isinstance(polarization, str) or polarization.upper() not in ("TE", "TM"):
            raise ConfigurationError("Specify TE (Hz) or TM (Ez).")
        self.polarization = polarization.upper()


class ScatteringSolver3D(_PlaneWaveSolver):
    """Full-vector Maxwell scattering by bounded three-dimensional targets."""
    def __init__(self, *, frequency, x_range, y_range, z_range, direction=(1, 0, 0),
                 electric_vector=(0, 1, 0), pml_thickness=None, pml_layers=5):
        self._configure(frequency, (x_range, y_range, z_range), direction, pml_thickness, pml_layers)
        self.polarization = "vector"
        self.electric_vector = _unit(electric_vector, 3, "electric_vector")
        if abs(self.direction @ self.electric_vector) > 1e-10:
            raise ConfigurationError("electric_vector must be transverse to direction.")
