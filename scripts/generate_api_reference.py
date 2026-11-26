"""Render the curated public API inventory from installed package signatures."""

from importlib import import_module
import argparse
import inspect
import json
import re
from pathlib import Path

from fem import __version__

ROOT = Path(__file__).resolve().parents[1]
DESCRIPTIONS = {
    "neff": "Effective indices with shape (frequencies, modes), or a vector for a single mode.",
    "show": "Open the plot window when true.",
    "connect": "Connect adjacent frequency samples when true; preserve mode ordering before connecting.",
    "frequency": "Operating frequency in hertz; finite and positive.",
    "x_range": "Physical x extent or increasing bounds, in metres.",
    "y_range": "Physical y extent or increasing bounds, in metres.",
    "z_range": "Physical z extent or increasing bounds, in metres.",
    "epsilon": "Relative permittivity; supported scalar/tensor forms are described below.",
    "mu": "Relative permeability; supported scalar/diagonal forms are described below.",
    "background_material": "Predefined bulk Material assigned to unfilled space.",
    "boundary": "Predefined PEC or PMC exterior-boundary material.",
    "material": "A predefined bulk, ideal-boundary, or supported SIBC material.",
    "geometry": "A handle returned by add_geometry or a geometry convenience method.",
    "max_element_size": "Maximum initial element edge length in metres.",
    "resolution": "Initial node counts; use instead of a maximum element size.",
    "wavelength_elements": "Minimum number of initial elements per local wavelength.",
    "element_order": "Finite-element polynomial order supported by this backend.",
    "quadrature_order": "Polynomial integration degree for finite-element assembly. Defaults to five.",
    "solve_options": "Additional keyword arguments accepted by this solver's solve method.",
    "num_modes": "Number of modes requested; positive integer.",
    "neff_guess": "Dimensionless complex effective-index search target.",
    "eigensolver_tolerance": "Algebraic eigensolver convergence tolerance.",
    "linear_solver_tolerance": "Algebraic linear-system residual tolerance.",
    "residual_tolerance": "Maximum accepted eigenproblem residual, separate from adaptation.",
    "divergence_tolerance": "Maximum accepted discrete Gauss-law residual.",
    "max_refinements": "Maximum mesh refinements after the initial solve; zero means one solve.",
    "adaptive_tolerance": "Relative discretization-residual stopping threshold.",
    "thickness": "Thickness in metres inside each selected outer face. Auto uses one third of the free-space wavelength.",
    "order": "Polynomial order of the PML profile.",
    "direction": "Propagation direction for solve; selected coordinate direction for PML.",
    "sigma_max": "Backend PML-strength magnitude; the outgoing stretch has negative imaginary sign.",
    "target_reflection": "Desired PML amplitude reflection ratio in (0, 1).",
    "component": "Field component to display, such as Ey; electrostatics also accepts potential or mesh.",
    "quantity": "Displayed field quantity: real, imag, magnitude/abs, or phase; static fields support real or magnitude.",
    "mode": "Zero-based mode index.",
    "case": "Zero-based sweep case index.",
    "block": "Wait for the interactive viewer to close when true.",
    "path": "Destination/source HDF5 path. Saving is atomic; loading does not run a solver.",
    "frequencies": "Strictly increasing positive frequencies in hertz.",
    "potential": "Prescribed electric potential in volts.",
    "outer_potential": "Exterior potential in volts; None permits natural boundaries.",
    "density": "Volume charge density in coulombs per cubic metre.",
    "region": "Geometry primitive or supported boundary name.",
    "shape": "A predefined fem.common.shapes object in metres.",
    "clip": "Intersect the shape with the solver domain; otherwise out-of-bounds objects raise GeometryError.",
    "name": "Optional name used for later identification and diagnostics.",
    "center": "Physical centre coordinates in metres.",
    "radius": "Positive radius in metres.",
    "points": "Ordered polygon vertex coordinates in metres.",
    "material_aware": "Use material-dependent initial mesh sizing.",
    "dim": "Electrostatic mesh dimension: 1 or 2.",
    "max_elements": "Adaptive mesh element budget.",
    "marking_fraction": "Fraction of squared error indicators marked for refinement.",
    "Zs": "Surface impedance in ohms; alternatively select a metal preset.",
    "preset": "Metal name for the good-conductor impedance model.",
    "results": "Nonempty sequence of completed periodic mode sets, in sweep order.",
    "number": "Zero-based mode index.",
    "num_points": "Number of field sampling points; positive integer.",
}
DESCRIPTIONS.update(
    {
        "backend": "Execution backend: auto (default) selects hybrid GPU Arnoldi/CPU Ritz when available, otherwise warns and selects SciPy CPU; hybrid, gpu (resident CUDA Ritz), and cpu explicitly select a backend.",
        "boundary_refinement": "Boundary mesh-size multiplier in (0, 1]; None disables boundary refinement.",
        "boundary_refinement_width": "Width in metres over which boundary mesh refinement transitions to the bulk size; None selects the mesher default.",
        "interface_refinement": "Material-interface mesh-size multiplier in (0, 1]; None disables interface refinement.",
        "interface_refinement_width": "Width in metres over which interface refinement transitions to the bulk size; None selects the mesher default.",
        "dense_limit": "Largest system dimension allowed for dense eigensolution.",
        "dense_linearization_limit": "Largest linearized eigenproblem dimension allowed for dense eigensolution.",
        "eigensolver": "Eigenproblem algorithm: auto, dense, or sparse where supported.",
        "factor": "Mesh refinement factor greater than one.",
        "polarization": "TE, TM, or both for the 2D periodic solver; result filtering selects the named polarization.",
        "transition_width": "Distance in metres over which local mesh sizing transitions to the background size.",
        "propagation_ratio_tolerance": "Minimum real-power fraction used to classify a mode as propagating.",
        "symmetry": "2D waveguide mode identification: None disables symmetry, auto detects a sampled Cnv subgroup (n=1..8), or declare a verified group such as C4v about the domain centre.",
        "degeneracy_tolerance": "Relative effective-index diameter of a candidate modal cluster, scaled by max(1, abs(neff)); basis rotations must also pass the original residual and divergence checks.",
        "track_modes": "Retain overlap-matched mode IDs and phases across 2D waveguide solves and remeshing; physical scene changes reset the reference.",
        "track_options": "Optional partial tracking configuration: symmetry (auto), degeneracy_tolerance (1e-5), overlap_threshold (0.5, squared overlap), and symmetry_tolerance (0.25). Omitted keys use defaults; track_modes=False disables tracking.",
        "max_pml_fraction": "Maximum permitted modal PML fraction; None disables PML-mode filtering.",
        "ncv": "Arnoldi subspace dimension; None lets the eigensolver choose.",
        "max_iterations": "CPU: maximum ARPACK iterations. GPU: maximum Arnoldi operator applications; None allows thirteen subspace cycles.",
        "random_seed": "Seed for reproducible eigensolver starting vectors.",
        "refine_interfaces": "Enable local refinement at dielectric interfaces.",
        "dielectric_refinement_factor": "Mesh-size multiplier near dielectric interfaces.",
        "pec_refinement_factor": "Mesh-size multiplier near PEC boundaries.",
        "pec_refinement_distance": "Distance in metres over which PEC mesh refinement applies.",
    }
)

TYPES = {
    "block": "bool",
    "case": "int",
    "mode": "int",
    "number": "int",
    "center": "tuple[float, ...] / m",
    "points": "sequence[tuple[float, ...]] / m",
    "component": "str | None",
    "quantity": "str",
    "path": "str | PathLike",
    "epsilon": "float | complex | array-like / relative",
    "mu": "float | complex | array-like / relative",
    "x_range": "float | tuple[float, float] / m",
    "y_range": "float | tuple[float, float] / m",
    "results": "Sequence[PeriodicModeSet]",
}


def section(name, character="-"):
    return name + "\n" + character * len(name) + "\n\n"


def entry(name, obj, returned, descriptions=None):
    signature = inspect.signature(obj)
    parameters = [p for p in signature.parameters.values() if p.name not in ("self", "cls")]
    signature = signature.replace(parameters=parameters)
    result = section(f"``{name}``", "~")
    result += ".. code-block:: python\n\n    " + name + str(signature) + "\n\n"
    doc = inspect.getdoc(obj)
    if doc:
        first = doc.split("\n\n", 1)[0].replace("\n", " ")
        first = re.sub(r"(?<!`)\|([^|`\n]+)\|(?!`)", r"``|\1|``", first)
        if not any(s in first for s in (":meth:", "legacy", "compatibility")):
            result += first + "\n\n"
    if parameters:
        result += ".. list-table:: Arguments\n   :header-rows: 1\n   :widths: 16 20 12 16 36\n\n   * - Argument\n     - Type / units\n     - Required / optional\n     - Default\n     - Meaning\n"
        for p in parameters:
            annotation = (
                str(p.annotation).replace("typing.", "")
                if p.annotation is not inspect.Parameter.empty
                else TYPES.get(
                    p.name,
                    type(p.default).__name__
                    if p.default is not inspect.Parameter.empty and p.default is not None
                    else "array-like or scalar",
                )
            )
            if annotation.startswith("<class '"):
                annotation = annotation[8:-2]
            default = "—" if p.default is inspect.Parameter.empty else repr(p.default)
            meaning = (descriptions or {}).get(p.name) or DESCRIPTIONS[p.name]
            required = p.default is inspect.Parameter.empty and p.kind != inspect.Parameter.VAR_KEYWORD
            result += f"   * - ``{p.name}``\n     - ``{annotation}``\n     - {'Required' if required else 'Optional'}\n     - ``{default}``\n     - {meaning}\n"
        result += "\n"
    return result + "Returns: " + returned + ".\n\n"


def render_references():
    inventory = json.loads((ROOT / "tests/public_api.json").read_text())
    for package, spec in inventory.items():
        module = import_module(package)
        out = section(package + " user API", "=")
        out += f"Version {__version__}. This reference covers the deliberately supported user API.\nImplementation helpers are documented in their source modules.\n\n"
        if package in ("fem.scattering", "fem.antenna"):
            family = package.split(".")[1]
            out += section("Plane-wave scattering" if family == "scattering" else "Waveguide-fed antennas")
            out += f"See the `solver guide <../solvers/{family}.rst>`_ for field conventions and examples.\n\n"
            out += ", ".join(f"``{name}``" for name in spec["exports"]) + ".\n\n"
            descriptions = {
                "direction": "Cartesian propagation direction, normalized internally.",
                "directions": "Cartesian observation directions, one per row.",
                "electric_vector": "Real incident electric polarization perpendicular to propagation.",
                "polarization": "TE (Hz) or TM (Ez) relative to the invariant z axis.",
                "pml_thickness": "Exterior PML thickness per side in metres. Defaults to one third of the free-space wavelength at the constructor frequency.",
                "pml_layers": "Number of mesh layers normal to each face of the physical box.",
                "distance": "Inset in metres from each PML interface, or 'auto' for a wavelength-based vacuum clearance.",
                "backend": "SciPy CPU or NVIDIA cuDSS (gpu or hybrid). Auto probes cuDSS.",
                "quantity": "Real part, imaginary part, or magnitude.",
                "path": "HDF5 result path for save, or .vtm path for export_vtk.",
                "component": "Incident, total or scattered field. Choose Et/Hz for 2D TE, Ez/Ht for 2D TM, or E/H in 3D (e.g. 'Total H'). RCS in either dimension, or Surface current on a 3D conductor.",
                "material": "Isotropic bulk material, PEC, surface impedance or good conductor.",
            }
            returns = {"mesh": "the physical mesh", "solve": "a ScatteringResult",
                       "add_geometry": "None", "add_mesh_refinement": "None", "add_pml": "the PML profile", "add_nf2ff": "None", "show": "the displayed figure or plotter",
                       "plot": "a Matplotlib figure in 2D or a PyVista plotter in 3D",
                       "vtk_datasets": "a VTK multiblock dataset",
                       "radar_cross_section": "an array of widths (m) or areas (m²)",
                       "save": "the HDF5 path", "export_vtk": "the VTM path"}
            if family == "antenna":
                descriptions.update(
                    name="Unique port name.", axis="Waveguide axis (x, y or z).",
                    position="Reference plane coordinate in metres along the guide axis.",
                    u_range="First transverse Cartesian coordinate interval in metres.",
                    v_range="Second transverse Cartesian coordinate interval in metres.",
                    cross_section="2D shape selecting the conforming aperture.",
                    num_modes="Retained mode count. Defaults to the selected mode rank.",
                    mode="One-based effective-index rank of the launched mode.",
                    gap="Box or cylinder spanning the two feed terminals along the port axis.",
                    impedance="Positive real reference impedance in ohms.",
                    amplitude="Incident power-wave magnitude in square roots of watts.",
                    phase_deg="Phase in degrees for the drive or displayed vector phasor.",
                    port="Name of the configured waveguide port.",
                    path="HDF5 result path.", quantity="Network magnitude (dB) or phase; real or imaginary vector field.",
                    component="S-parameters, Total E, Total H, Surface current, Directivity, Gain or Realized gain.",
                )
                returns.update(solve="an AntennaResult", sweep="an AntennaSweep",
                               add_waveguide_port="None", add_discrete_port="None", solve_port_modes="a 2D waveguide ModeSet",
                               show_port_modes="the 2D waveguide mode viewer", set_excitation="None",
                               far_field="a RadiationPattern", plot="a Matplotlib figure or VTK plotter")
            for cls_name, methods in (spec["solvers"] | spec["results"]).items():
                cls = getattr(module, cls_name)
                for method in methods:
                    if method == "__init__":
                        out += entry(cls_name, cls, "the configured solver", descriptions)
                    else:
                        out += entry(f"{cls_name}.{method}", getattr(cls, method),
                                     returns[method], ({**descriptions, "direction": "One face, an axis pair, or all exterior faces."}
                                                       if method == "add_pml" else descriptions))
            yield ROOT / spec["documentation"] / f"{family}.rst", out.rstrip() + "\n"
            continue
        out += section("Workflow and units")
        out += (
            "All lengths are in metres, frequencies in hertz, and constitutive values are relative.\n"
            "Construction and configuration use keyword arguments. Call ``mesh()`` to build\n"
            "the initial mesh, ``solve()`` to obtain and store a typed result, and ``show()``\n"
            "to inspect it interactively. ``mesh_data`` and ``result`` are initially None.\n"
            "Geometry edits invalidate both; remeshing invalidates the result. Automatic\n"
            "meshing reuses the last explicit settings. Calling ``show()`` without a result\n"
            "raises ``NoResultError``.\n\n"
            "Adaptive defaults are ``max_refinements=2`` and ``adaptive_tolerance=0.05``.\n"
            "Zero refinements performs one solve. ``solve_info`` distinguishes algebraic\n"
            "residuals from adaptive discretization residuals and records stopping reasons.\n"
            "Python mode and case indices are zero-based.\n\n"
            "``solve()`` and frequency sweeps neither save nor open windows. Call\n"
            "``result.save(path)`` explicitly. ``load_result(path)`` returns inspection-ready\n"
            "results without solving. Archives use ``fem-results`` schema ``1.0``;\n"
            "old or convention-incompatible files raise ``PersistenceError``. Loaded\n"
            "results can plot, show, and save; they do not restart a solver or restore callbacks.\n"
            "``plot()`` returns a Matplotlib Figure without opening a window.\n"
            "``show(block=True)`` waits for the viewer; ``block=False`` returns immediately.\n\n"
        )
        if package == "fem.electrostatics":
            out += "Electrostatic fields are static. ``epsilon`` accepts a positive scalar, a\npositive diagonal, or a real symmetric positive-definite tensor. Potential is\nnodal; ``element_electric_field`` and ``element_displacement_field`` are cell\nfields. ``electric_field`` and ``displacement_field`` are nodal averages.\n\n"
        else:
            out += "The time convention is ``exp(+i*omega*t)`` with guided propagation\n``exp(-i*beta*z)``. Passive materials have nonpositive imaginary constitutive\nvalues; passive forward attenuation is ``-Im(beta)``.\n\n"
            out += "Materials support scalar or diagonal epsilon and mu. Unsupported off-diagonal\nforms raise an explicit configuration error. Periodic archives store periodic\nfield envelopes; the Bloch phase is recorded separately from these fields.\n\n"
        out += (
            section("Supported exports")
            + ", ".join(f"``{name}``" for name in spec["exports"])
            + ".\n\n"
        )
        out += section("Solver construction and operations")
        for solver, methods in spec["solvers"].items():
            cls = getattr(module, solver)
            for method in methods:
                name = solver if method == "__init__" else solver + "." + method
                obj = cls if method == "__init__" else getattr(cls, method)
                returned = {
                    "__init__": f"a configured ``{solver}``",
                    "mesh": "the initial mesh stored in ``mesh_data``",
                    "solve": "the physics-specific result stored in ``result``",
                    "show": "the interactive viewer window",
                    "sweep": "an iterator of mode sets in increasing frequency order",
                }.get(
                    method,
                    "the configured geometry/excitation handle, or None for in-place configuration",
                )
                out += entry(name, obj, returned)
        out += section("Returned results")
        out += "Result objects are returned by solving or loading. Their constructors are\nimplementation details. Inspect ``mesh_data``, ``metadata``, ``solve_info``, and\n``frequency`` where applicable. Modal sets support iteration, indexing, ``neff``\nand ``beta``; each mode provides sampled fields and residual diagnostics.\n\n"
        for kind, methods in spec["results"].items():
            cls = getattr(module, kind)
            for method in methods:
                out += entry(
                    kind + "." + method,
                    getattr(cls, method),
                    {
                        "save": "the written Path",
                        "plot": "a Matplotlib Figure",
                        "show": "the interactive viewer window",
                        "mode": "the selected mode",
                        "deembed": "a result with updated reference planes",
                        "conductor_charge": "conductor charge in coulombs (per metre for a 2D cross-section)",
                    }.get(method, "the selected data or diagnostic report"),
                )
        if package == "fem.waveguide":
            out += section("ModeSolver2D tracking options")
            out += """``ModeSolver2D.solve(track_modes=True)`` enables symmetry identification,
degenerate-subspace handling, and overlap matching using the defaults below.
``track_options`` is optional: omit it, pass ``None`` or ``{}``, or supply a
mapping containing only the settings to override. The mapping is not modified.

.. list-table:: Supported track_options keys
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning and allowed values
   * - ``symmetry``
     - ``"auto"``
     - Sample the continuous scene for a centred Cnv subgroup with coordinate-aligned mirrors (n=1..8). Use ``None`` for overlap tracking without symmetry, or declare ``"C1v"`` through ``"C8v"``; declared groups are also verified.
   * - ``degeneracy_tolerance``
     - ``1e-5``
     - Finite positive cluster diameter in effective index, scaled by max(1, abs(neff)). This groups candidate modes; it does not replace eigenproblem accuracy tolerances.
   * - ``overlap_threshold``
     - ``0.5``
     - Minimum squared individual or principal-angle overlap in (0, 1] for retaining mode IDs and aligning modal subspaces. Unmatched modes receive new IDs.
   * - ``symmetry_tolerance``
     - ``0.25``
     - Finite positive upper bound on the sampled symmetry-subspace residual for assigning a resolved irrep and mirror partner. Unresolved labels remain ``None``.

.. code-block:: python

    # Defaults; no track_options argument is required.
    modes = solver.solve(num_modes=24, track_modes=True)

    # Partial configuration; all other tracking settings keep their defaults.
    modes = solver.solve(num_modes=24, track_modes=True, track_options={
        "symmetry": "C4v",
        "degeneracy_tolerance": 1e-4,
        "overlap_threshold": 0.8,
    })

    # Disable identification and clear the tracking reference.
    modes = solver.solve(num_modes=24, track_modes=False)

``track_modes=False`` disables tracking even if ``track_options`` is supplied.
Options are still validated; unknown keys, invalid values, or non-mapping
configurations raise ``ConfigurationError``. Effective settings are stored in
``result.metadata["track_options"]`` when tracking is enabled.

Reuse the same solver to retain mode IDs and phases across remeshing and adaptive
passes. Geometry changes reset the reference. Loading an archive restores its
diagnostics, not a live tracking session. Each tracked mode records ``mode_id``,
``symmetry_irrep``, ``symmetry_partner``, ``tracking_overlap``, ``subspace_overlap``,
and longitudinal electric/magnetic fractions in ``mode.metadata``.
Multidimensional irreps retain their common label and separate mirror partners.

Rotated eigenvectors must pass the ordinary ``residual_tolerance`` and
``divergence_tolerance`` arguments to ``solve()``. Otherwise valid eigenvectors
are retained and ``symmetry_projection="residual_rejected"`` is recorded.
See the `waveguide guide <../solvers/waveguide.rst>`_ for detector and sampling
limits and the meaning of unresolved or truncated multiplets.

"""
        if package in ('fem.waveguide', 'fem.periodic'):
            out += section("Transverse vector views")
            out += """In the interactive viewer opened by ``show()``, select ``Et`` or ``Ht``
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

"""
        out += section("Result data and diagnostics")
        out += "``mesh_data.coordinates`` stores physical nodes in metres; ``elements`` stores\nzero-based connectivity. ``axes`` identifies physical coordinate order.\n``mesh_data.metadata['context']`` records material and boundary configuration.\nThe result is an inspection snapshot; editing it cannot restart a solver.\n\n"
        if package == "fem.electrostatics":
            out += "``potential`` is nodal potential in V. ``electric_field`` and\n``displacement_field`` are recovered nodal arrays (V/m and C/m²);\n``element_electric_field`` and ``element_displacement_field`` retain cell fields.\n``conductor_charges`` maps configured conductor names to charge, and\n``energy`` records electrostatic energy (per unit transverse area in 1D,\nper unit length in 2D). ``residual_norm`` is the algebraic residual;\n``adaptive_history`` records mesh sizes, error indicators, and stopping status.\n\n"
        else:
            out += "Each selected mode exposes ``neff`` (dimensionless), ``beta`` (rad/m),\n``fields`` (sampled E/H in V/m and A/m), ``coefficients`` (FEM expansion),\nand ``residual`` (algebraic eigenproblem residual). Field component names\ninclude Ex, Ey, Ez, Hx, Hy, and Hz. Coordinate and cell ownership arrays\npreserve field locations; they must not be interpreted as interchangeable\nnodal or cell values. ``solve_info`` records adaptive history separately.\n\n"
            if package == "fem.periodic":
                out += "Periodic fields are Bloch envelopes. ``period`` is in m; each mode also\nprovides ``gamma``, ``bloch_multiplier``, folded propagation quantities,\nand Gauss-law/PML filtering diagnostics. Combine solved cases with\n``PeriodicSweepResult.from_results(results)`` and call ``save(path)``.\nLoaded multi-case archives index cases lazily.\n\n"
        out += section("Geometry and material values")
        out += (
            "Define reusable materials and shapes with ``fem.common`` before assigning them.\n"
            "Use ``Material(name=..., epsilon=..., mu=...)`` for bulk media,\n"
            "``materials.PEC`` or ``materials.PMC`` for ideal boundaries, and the\n"
            "documented ``materials.copper``-style presets where SIBC is supported.\n"
            "Continuous primitives and Boolean/transformed shapes live in\n"
            "``fem.common.shapes``. Solver packages also re-export Material, materials, and shapes.\n\n"
        )
        out += entry("fem.load", import_module("fem").load,
                     "a typed result; multi-case archives provide lazy case access")
        out += (
            section("Errors")
            + "Invalid inputs raise ``ConfigurationError`` or ``GeometryError`` where available.\nMesh and numerical failures raise the corresponding ``MeshError`` or\n``SolverError``. ``NoResultError`` requires a successful solve first.\n``PersistenceError`` identifies an incompatible or unreadable archive.\nViewer errors include the executable path or installation setting needed to\ncorrect a launch failure. Saving and loading do not require an active GUI.\n\n"
        )
        out += section("Backend selection")
        out += (
            '``backend="auto"`` tries hybrid GPU Arnoldi/CPU Ritz, then warns before SciPy CPU fallback.\n'
            '``hybrid`` explicitly selects that path; ``gpu`` explicitly requires resident CUDA Ritz.\n'
            "Both accelerated choices use cuDSS for direct scattering solves.\n\n"
        )
        yield ROOT / spec["documentation"] / (package.split(".")[1] + ".rst"), out.rstrip() + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="Fail if generated references differ; write nothing."
    )
    args = parser.parse_args()
    stale = []
    for path, content in render_references():
        if args.check:
            if not path.exists() or path.read_text(encoding="utf-8") != content:
                stale.append(str(path.relative_to(ROOT)))
        else:
            path.write_text(content, encoding="utf-8")
    if stale:
        raise SystemExit("Run python scripts/generate_api_reference.py: " + ", ".join(stale))


if __name__ == "__main__":
    main()
