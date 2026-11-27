"""Driven three-dimensional edge-element antenna model with multimode wave ports."""

from dataclasses import dataclass, replace

import numpy as np
from scipy.sparse import bmat, eye, hstack, vstack
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu
from skfem import Basis, BilinearForm, FacetBasis, asm
from skfem.helpers import curl, dot

from fem.common import materials, shapes
from fem.common.constants import C_0, ETA_0
from fem.common.discretization import DEFAULT_QUADRATURE_ORDER, field_orders
from fem.adaptivity import run_adaptive
from fem.adaptivity.driven import AbsorberCheck, driven_residual
from fem.common.contracts import bounds
from fem.common.errors import ConfigurationError, GeometryError, SolverError
from fem.common.linalg import resolve_backend
from fem.common.materials import bulk_values
from fem.common.radiation_pml import shell_matrix
from fem.common.progress import frequency_steps
from fem.common.tet_hcurl import conductor_current, tetrahedral_element
from fem.common.validation import finite_real, positive_integer, positive_real, real_array
from fem.scattering.meshing import RadiationBoxSolver, build_mesh

from .farfield import HuygensSurface
from .ports import DiscretePort, WaveguidePort
from .results import AntennaResult, AntennaSweep


def _recover_fields(basis, coeff, mu, k, conductors):
    """Cell means and the conductor-side surface current from the FE curl."""
    electric = basis.interpolate(coeff)
    weights = basis.W / basis.W.sum()
    field = np.einsum("ieq,q->ei", electric, weights)
    mean_curl = np.einsum("ieq,q->ei", electric.curl, weights)
    magnetic = (1j / (k * ETA_0)) * mean_curl / mu[:, None]
    # Curl is affine at order two. Three one-sided corner values capture its
    # entire face trace, including discontinuities at the edges of conductors.
    vertices = conductor_current(basis, coeff, mu, k, conductors)
    return field, magnetic, vertices.mean(axis=1), vertices if basis.elem.maxdeg == 2 else None


def _drive(amplitude, phase_deg):
    amplitude = finite_real(amplitude, "amplitude")
    if amplitude < 0:
        raise ConfigurationError("Wave amplitude in sqrt(W) must be nonnegative.")
    return amplitude, finite_real(phase_deg, "phase_deg")


@dataclass(frozen=True)
class _SParameterConvergence:
    tolerance: float
    minimum_passes: int
    consecutive_checks: int

    def __call__(self, result, history, record):
        matrix = np.asarray(result.s_parameters)
        if not matrix.size or not np.isfinite(matrix).all():
            raise SolverError("S-parameter convergence requires a finite, nonempty matrix.")
        previous = history[-1] if history else None
        physical_pass, checks, delta = 1, 0, None
        if previous is not None:
            physical_pass = previous["physical_pass"]
            checks = previous["s_parameter_checks"]
            # The absorber has its own convergence check on a fixed physical mesh.
            if previous["status"] != "pml_refine":
                physical_pass += 1
                reference = previous["s_parameter_matrix"]
                if reference.shape != matrix.shape:
                    raise SolverError("Port matrix shape changed during mesh refinement.")
                delta = float(np.max(np.abs(matrix - reference)))
                checks = checks + 1 if delta < self.tolerance else 0
        record.update(physical_pass=physical_pass, s_parameter_checks=checks,
                      s_parameter_delta=delta, s_parameter_matrix=matrix.copy(),
                      s_parameter_tolerance=self.tolerance, minimum_passes=self.minimum_passes,
                      consecutive_checks=self.consecutive_checks)
        return physical_pass >= self.minimum_passes and checks >= self.consecutive_checks


class AntennaSolver3D(RadiationBoxSolver):
    """Driven antenna in a physical vacuum box surrounded by a Cartesian PML."""

    def __init__(self, *, frequency, x_range, y_range, z_range, pml_thickness=None, pml_layers=5):
        self.frequency = positive_real(frequency, "frequency")
        self._configure_pml(pml_thickness, pml_layers)
        self._ranges = tuple(
            bounds(value, f"{axis}_range")
            for axis, value in zip("xyz", (x_range, y_range, z_range))
        )
        self._regions = []
        self._refinements = []
        self._ports = {}
        self._s_parameter_check = _SParameterConvergence(.02, 3, 1)
        self._s_parameter_configured = False
        self._wavelength_elements = None
        self._cad_point_size = None
        self._mesh_frequency = None
        self._size = None
        self._mesh = self.mesh_data = self.result = None

    def add_geometry(self, *, shape, material):
        """Add an isotropic bulk, PEC, or impedance target inside the radiation box."""
        if not isinstance(shape, shapes.Shape) or shape.dimension != 3:
            raise GeometryError("Antenna geometry must be three-dimensional.")
        self._check_geometry_clearance(shape.bounds, vacuum=self._vacuum(material))
        if isinstance(material, materials.Material):
            eps, mu = bulk_values(material, form="scalar")
            if mu == 0 or not material.is_passive:
                raise ConfigurationError(
                    "Antenna power accounting requires passive media and nonzero mu."
                )
        elif material != materials.PEC and not isinstance(
            material, (materials.SurfaceImpedance, materials.GoodConductor)
        ):
            raise ConfigurationError("Expected a bulk material, PEC, or surface impedance.")
        self._regions.append((shape, material))
        self._mesh = self.mesh_data = self.result = None

    def _enclosed_bounds(self):
        return super()._enclosed_bounds() + [port.bounds for port in self._ports.values()]

    def _validate_port(self, name, axis):
        if not isinstance(name, str) or not name.strip() or name in self._ports:
            raise ConfigurationError("Each port needs a unique nonempty name.")
        if axis not in ("x", "y", "z"):
            raise ConfigurationError("Port axis must be x, y, or z.")

    def add_waveguide_port(
        self,
        *,
        name,
        axis,
        position,
        cross_section,
        mode=1,
        num_modes=None,
        amplitude=1.0,
        phase_deg=0.0,
    ):
        """Assign a numerical wave port. Mode 1 is the highest effective-index forward mode."""
        self._validate_port(name, axis)
        if not isinstance(cross_section, shapes.Shape) or cross_section.dimension != 2:
            raise GeometryError("cross_section must be a 2D aperture shape.")
        rank = positive_integer(mode, "mode")
        count = rank if num_modes is None else positive_integer(num_modes, "num_modes")
        if rank > count:
            raise ConfigurationError("Selected mode exceeds num_modes.")
        port = WaveguidePort(
            name,
            axis,
            finite_real(position, "position"),
            cross_section,
            count,
            rank,
            *_drive(amplitude, phase_deg),
        )
        self._check_geometry_clearance(port.bounds)
        self._ports[name] = port
        self._mesh = self.mesh_data = self.result = None

    def add_discrete_port(self, *, name, axis, gap, impedance=50.0, amplitude=1.0, phase_deg=0.0):
        """Assign a box or cylindrical feed gap between conductors with a real reference impedance."""
        self._validate_port(name, axis)
        port = DiscretePort(
            name, axis, gap, positive_real(impedance, "impedance"), *_drive(amplitude, phase_deg)
        )
        self._check_geometry_clearance(port.bounds)
        self._ports[name] = port
        self._mesh = self.mesh_data = self.result = None

    def solve_port_modes(self, *, port, frequency=None):
        """Compute the selected port's 2D ModeSet for inspection before a driven solve."""
        if not isinstance(self._ports.get(port), WaveguidePort):
            raise ConfigurationError(f"Expected a waveguide port, got {port!r}.")
        f = self.frequency if frequency is None else positive_real(frequency, "frequency")
        _, _, sections, _ = self._prepare(f)
        return sections[port].solve_modes(f, self._properties())

    def set_excitation(self, *, port, amplitude, phase_deg=0.0):
        """Set the incident power-wave magnitude in sqrt(W) and phase in degrees."""
        if port not in self._ports:
            raise ConfigurationError(f"Unknown port {port!r}.")
        amplitude, phase_deg = _drive(amplitude, phase_deg)
        self._ports[port] = replace(self._ports[port], amplitude=amplitude, phase_deg=phase_deg)
        self.result = None

    def configure_s_parameter_convergence(self, *, tolerance=None, minimum_passes=3, consecutive_checks=1):
        """Stop physical adaptation on the maximum change in the complex port matrix."""
        minimum = positive_integer(minimum_passes, "minimum_passes", minimum=2)
        checks = positive_integer(consecutive_checks, "consecutive_checks")
        self._s_parameter_check = (None if tolerance is None else _SParameterConvergence(
            positive_real(tolerance, "tolerance"), minimum, checks))
        self._s_parameter_configured = True

    def mesh(self, *, max_element_size=None, wavelength_elements=None, cad_point_size=None, element_order=2,
             quadrature_order=DEFAULT_QUADRATURE_ORDER):
        if wavelength_elements is not None:
            if max_element_size is not None:
                raise ConfigurationError("Choose wavelength_elements or max_element_size.")
            wavelength_elements = positive_real(wavelength_elements, "wavelength_elements")
        if cad_point_size is not None:
            cad_point_size = positive_real(cad_point_size, "cad_point_size")
        if self._adaptive_sizes is not None:
            return self._mesh_at(self._mesh_frequency, max_element_size, element_order, quadrature_order,
                                 self._wavelength_elements, self._cad_point_size)
        return self._mesh_at(self.frequency, max_element_size, element_order, quadrature_order,
                             wavelength_elements, cad_point_size)

    def _mesh_at(self, frequency, max_element_size, element_order, quadrature_order,
                 wavelength_elements, cad_point_size):
        element_order, quadrature_order = field_orders(element_order, quadrature_order, (1, 2))
        if not self._regions:
            raise GeometryError("Add antenna and feed geometry before meshing.")
        self._surface = self._nf2ff_surface()
        policy, size, point, local = self._driven_mesh_settings(
            C_0 / frequency, element_order, max_element_size, wavelength_elements, cad_point_size)
        material_aware = max_element_size is None or (self._adaptive_sizes is not None and self._mesh_auto_materials)
        physical, snapshot, faces = build_mesh(
            self._physical_ranges(self._pml_profiles),
            self._regions,
            size,
            interior_box=self._surface,
            refinements=(*self._refinements, *local),
            adaptive_sizes=self._adaptive_sizes,
            material_aware=material_aware,
            cad_point_size=point,
            port_gaps={p.name: p.gap for p in self._ports.values() if isinstance(p, DiscretePort)},
        )
        snapshot = replace(snapshot, info={**snapshot.info, "mesh_profile": policy.name,
                           "cad_point_size": point, "automatic_local_sizes": [h for _, h in local]})
        self._attach_pml(physical, snapshot, faces, element_order, quadrature_order)
        self._apply_driven_policy(policy)
        self._mesh_auto_materials = material_aware
        self._wavelength_elements, self._mesh_frequency = wavelength_elements, frequency
        self._cad_point_size = cad_point_size
        self._size, self.result = size, None
        return self.mesh_data

    def _prepare(self, max_frequency):
        if not self._ports:
            raise ConfigurationError("Add at least one waveguide or discrete port.")
        if self._mesh is None or (self._mesh_auto_materials and self._mesh_frequency < max_frequency):
            self._mesh_at(max_frequency, None if self._mesh_auto_materials else self._size,
                          wavelength_elements=self._wavelength_elements, cad_point_size=self._cad_point_size,
                          **self._mesh_orders)
        quadrature = self._mesh_orders["quadrature_order"]
        basis = Basis(
            self._mesh, tetrahedral_element(self._mesh_orders["element_order"]),
            elements=np.arange(len(self.mesh_data.elements)), intorder=quadrature
        )
        wall_groups = [(self._regions[i - 1][1], f) for i, f in self._faces.items() if i]
        metal = np.concatenate([f for _, f in wall_groups] or [np.empty(0, dtype=int)])
        conductors = metal.copy()
        sections = {}
        boundaries = dict(self.mesh_data.metadata["boundary_facets"])
        tags = self.mesh_data.metadata["element_tags"]
        length = max(b - a for a, b in self._ranges)
        for port in self._ports.values():
            if isinstance(port, DiscretePort):
                cells = self.mesh_data.metadata["port_gap_cells"][port.name]
                sections[port.name] = port.voltage_trace(basis, metal, cells)
                continue
            # A reference face replaces its cap condition. The edges inherit the
            # side-wall material from the faces that continue along the guide.
            axis = "xyz".index(port.axis)
            impedance_edges = {}
            for wall, faces in wall_groups:
                if wall == materials.PEC:
                    continue
                vertices = self._mesh.p[axis, self._mesh.facets[:, faces]]
                side_faces = faces[np.ptp(vertices, axis=0) > length * 1e-9]
                edges = np.unique(self._mesh.f2e[:, side_faces])
                impedance_edges.update(
                    (tuple(sorted(pair)), wall) for pair in self._mesh.edges[:, edges].T
                )
            section = port.section(
                basis,
                conductors,
                length,
                tags,
                self.mesh_data.metadata["physical_names"],
                impedance_edges,
                quadrature,
            )
            facets = section.face.find
            conductors = np.setdiff1d(conductors, facets)
            boundaries[f"port:{port.name}"] = self._mesh.facets[:, facets].T
            sections[port.name] = section
        pec = np.concatenate(
            [f for mat, f in wall_groups if mat == materials.PEC] or [np.empty(0, dtype=int)]
        )
        remaining = np.intersect1d(pec, conductors)
        boundaries["pec"] = self._mesh.facets[:, remaining].T
        boundaries["sibc"] = self._mesh.facets[:, np.setdiff1d(conductors, remaining)].T
        self.mesh_data = replace(
            self.mesh_data, metadata={**self.mesh_data.metadata, "boundary_facets": boundaries}
        )
        fixed = basis.get_dofs(facets=np.union1d(remaining, self._pml_wall)).all()
        return basis, basis.complement_dofs(fixed), sections, conductors

    def _adapt_frequency(self, frequency, selected, max_refinements, adaptive_tolerance,
                         max_elements, max_dofs, pml_tolerance, pml_observable, max_pml_refinements):
        prepared = None
        check = self._s_parameter_check
        if adaptive_tolerance is None:
            adaptive_tolerance = check.tolerance if check is not None else .02
        elif not self._s_parameter_configured:
            # An explicit discretization tolerance requests the field-residual criterion.
            check = None

        def solve_once():
            nonlocal prepared
            prepared = self._prepare(frequency)
            self.result = self._solve(frequency, selected, *prepared, estimate=True)
            return self.result

        absorber = AbsorberCheck(adaptive_tolerance if pml_tolerance is None else pml_tolerance,
                                 pml_observable, max_pml_refinements)
        result = run_adaptive(self, solve_once, lambda result: self._estimate, self._refine_driven,
                              max_refinements=max_refinements, adaptive_tolerance=adaptive_tolerance,
                              max_elements=max_elements, max_dofs=max_dofs,
                              error_type=ConfigurationError, absorber=absorber,
                              physical_check=check)
        return result, prepared

    def solve(self, *, backend="auto", max_refinements=8, adaptive_tolerance=None,
              max_elements=None, max_dofs=None, pml_tolerance=None, pml_observable=None,
              max_pml_refinements=3):
        """Adapt against all independent port excitations, then report their coherent superposition."""
        self.result, _ = self._adapt_frequency(
            self.frequency, resolve_backend(backend), max_refinements, adaptive_tolerance,
            max_elements, max_dofs, pml_tolerance, pml_observable, max_pml_refinements
        )
        return self.result

    def sweep(self, *, frequencies, backend="auto", max_refinements=8, adaptive_tolerance=None,
              max_elements=None, max_dofs=None, pml_tolerance=None, pml_observable=None,
              max_pml_refinements=3):
        """Adapt at the highest frequency in hertz and reuse that mesh at every lower frequency."""
        fs = real_array(frequencies, "frequencies")
        if fs.ndim != 1 or not len(fs) or (fs <= 0).any() or (np.diff(fs) <= 0).any():
            raise ConfigurationError("Sweep frequencies must be positive and strictly increasing.")
        selected = resolve_backend(backend)
        lower = []
        for frequency in frequency_steps(fs):
            if frequency == fs[-1]:
                highest, prepared = self._adapt_frequency(
                    frequency, selected, max_refinements, adaptive_tolerance, max_elements, max_dofs,
                    pml_tolerance, pml_observable, max_pml_refinements)
                highest.metadata["mesh_adaptation_frequency"] = frequency
            else:
                case = self._solve(frequency, selected, *prepared)
                case.metadata["mesh_adaptation_frequency"] = highest.frequency
                lower.append(case)
        self.result = AntennaSweep((*lower, highest))
        return self.result

    def _properties(self):
        return [(1.0, 1.0)] + [
            bulk_values(m, form="scalar") if isinstance(m, materials.Material) else (1.0, 1.0)
            for _, m in self._regions
        ]

    def _solve(self, frequency, backend, basis, free, sections, conductors, *, estimate=False):
        k = 2 * np.pi * frequency / C_0
        quadrature = self._mesh_orders["quadrature_order"]
        tags = self.mesh_data.metadata["element_tags"]
        props = self._properties()
        eps, mu = np.asarray(props, dtype=complex)[tags].T

        @BilinearForm(dtype=complex)
        def interior(u, v, w):
            return dot(curl(u), curl(v)) / mu[:, None] - k * k * eps[:, None] * dot(u, v)

        @BilinearForm(dtype=complex)
        def impedance(u, v, w):
            return 1j * k * (dot(u, v) - dot(u, w.n) * dot(v, w.n))

        matrix = asm(interior, basis) + shell_matrix(
            basis, len(tags), self._physical_ranges(self._pml_profiles), self._pml_profiles, k
        )
        natural_walls = []
        for owner, facets in self._faces.items():
            material = self._regions[owner - 1][1]
            if material == materials.PEC:
                continue
            facets = np.intersect1d(facets, conductors)
            if not len(facets):
                continue
            ratio = ETA_0 / material.at_frequency(frequency=frequency)
            face = FacetBasis(self._mesh, basis.elem, facets=facets, intorder=quadrature)
            matrix += ratio * asm(impedance, face)
            natural_walls.append((facets, ratio))

        ports = tuple(self._ports.values())
        rhs = np.zeros((basis.N, len(ports)), dtype=complex)
        outputs = {}
        port_modes, flux_operators = {}, {}
        loads, traces = [], []
        factor = 1j * k * ETA_0
        port_profiles = {}
        for col, port in enumerate(ports):
            section = sections[port.name]
            if isinstance(port, DiscretePort):
                voltage = section
                scale = np.sqrt(2 * port.impedance)
                loads.append((factor / port.impedance) * voltage)
                traces.append(voltage.T)
                rhs[:, col] = (2 * factor * scale / port.impedance) * voltage.toarray().ravel()
                outputs[port.name, port.mode] = voltage.T / scale
                continue
            modes = section.solve_modes(frequency, props)
            if modes[port.mode - 1].metadata["classification"] != "propagating":
                raise ConfigurationError(
                    f"{port.name} mode {port.mode} is below cutoff at {frequency:g} Hz."
                )
            port_modes[port.name] = modes
            load, projector, flux, profiles = section.operators(modes)
            port_profiles[col] = profiles
            loads.append(factor * load)
            traces.append(projector)
            rhs[:, col] = 2 * factor * load[:, port.mode - 1].toarray().ravel()
            flux_operators[col] = projector, flux
            for rank, mode in enumerate(modes, 1):
                if mode.metadata["classification"] == "propagating":
                    outputs[port.name, rank] = projector[rank - 1, :]

        # Each port imposes a low-rank update L P E to the volume equation.
        # Introduce its modal amplitudes y = P E as unknowns. Eliminating y
        # recovers exactly the DtN/gap operator, without a dense clique across
        # every pair of port DOFs (especially expensive on a refined gap).
        load = hstack(loads, format="csr")[free]
        trace = vstack(traces, format="csc")[:, free]
        reduced = bmat([[matrix[free][:, free], load],
                        [-trace, eye(trace.shape[0])]], format="csr")
        reduced_rhs = np.vstack((rhs[free], np.zeros((trace.shape[0], len(ports)), dtype=complex)))
        # A closed PEC guide can isolate its ports completely from the exterior.
        # Source-free disconnected blocks have zero driven field, including in the PML.
        count, component = connected_components(reduced, directed=False)
        if count > 1:
            driven = np.any(reduced_rhs != 0, axis=1)
            active = np.isin(component, np.unique(component[driven]))
            reduced = reduced[active][:, active]
            reduced_rhs = reduced_rhs[active]
            free = free[active[:len(free)]]
        if backend == "cpu":
            solutions = splu(reduced.tocsc()).solve(reduced_rhs)
        else:
            from fem.common.gpu import solve_sparse

            solutions = solve_sparse(reduced, reduced_rhs)
        if not np.isfinite(solutions).all():
            raise SolverError(f"Nonfinite antenna solution at {frequency:g} Hz.")
        solutions = solutions[:len(free)]
        incident = np.array([p.drive for p in ports])
        coeff = np.zeros(basis.N, dtype=complex)
        if estimate:
            combined = np.zeros(len(tags))
            worst = 0.
            worst_components = {}
            absorber = self._physical_ranges(self._pml_profiles), self._pml_profiles
            domain_report = {}
            cache = {} if len(ports) > 1 else None
            for excitation, column in enumerate(solutions.T):
                coeff[free] = column
                port_traces, gap_sources = [], []
                for col, port in enumerate(ports):
                    if isinstance(port, DiscretePort):
                        cells = self.mesh_data.metadata["port_gap_cells"][port.name]
                        axis = "xyz".index(port.axis)
                        span = np.diff(np.asarray(port.bounds).reshape(3, 2)[axis]).item()
                        area = basis.dx[cells].sum() / span
                        voltage = (sections[port.name].T @ coeff).item()
                        drive = 2 * np.sqrt(2 * port.impedance) if col == excitation else 0.
                        source = 1j * ETA_0 * (drive - voltage) / (port.impedance * area)
                        gap_sources.append((cells, axis, source))
                    else:
                        amplitudes = flux_operators[col][0] @ coeff
                        if col == excitation:
                            amplitudes[port.mode - 1] -= 2
                        prescribed = 1j * ETA_0 * np.einsum("m,mifq->ifq", amplitudes, port_profiles[col])
                        port_traces.append((sections[port.name].face, prescribed))
                indicators, residual, components, report = driven_residual(
                    basis, coeff, eps, 1 / mu, k, quadrature, walls=natural_walls,
                    gap_sources=gap_sources, port_traces=port_traces, absorber=absorber, cache=cache,
                    marking_weights=self._marking_weights)
                for key, value in report.items():
                    domain_report[key] = max(value, domain_report.get(key, 0.))
                if not np.isfinite(residual):
                    raise SolverError("Nonfinite adaptive residual for an antenna port excitation.")
                combined = np.maximum(combined, indicators)
                if residual >= worst:
                    worst, worst_components = residual, components
            self._estimate = combined, worst, worst_components, domain_report
        coeff[free] = solutions @ incident
        channels = tuple((p.name, p.mode) for p in ports)
        modal_s = {}
        for key, projector in outputs.items():
            row = (projector[:, free] @ solutions).ravel()
            if key in channels:
                row[channels.index(key)] -= 1.0
            modal_s[key] = row
        accepted = 0.0
        for col, port in enumerate(ports):
            if isinstance(port, DiscretePort):
                b = modal_s[channels[col]] @ incident
                accepted += abs(incident[col]) ** 2 - abs(b) ** 2
            else:
                projector, flux = flux_operators[col]
                total = projector @ coeff
                inward = -total.copy()
                inward[port.mode - 1] += 2 * incident[col]
                accepted += np.vdot(inward, flux @ total).real
        field, magnetic, current, current_vertices = _recover_fields(basis, coeff, mu, k, conductors)
        surface = HuygensSurface.sample(self._mesh, basis.elem, self._surface, coeff, k, quadrature)
        return AntennaResult(
            frequency,
            channels,
            incident,
            modal_s,
            surface,
            surface.radiated_power(),
            self.mesh_data,
            field,
            port_modes,
            self._mesh.facets[:, conductors].T,
            current,
            float(accepted),
            {
                "backend": backend,
                "ndofs": basis.N,
                "reference_impedances": {
                    p.name: p.impedance for p in ports if isinstance(p, DiscretePort)
                },
                "method": "scipy-direct" if backend == "cpu" else "cudss-direct",
            },
            magnetic_field=magnetic,
            surface_current_vertices=current_vertices,
        )
