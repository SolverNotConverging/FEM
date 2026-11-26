"""Frequency sweeps of guided and Bloch modes on one physical mesh."""
from copy import copy
from dataclasses import replace

import numpy as np

from .constants import C_0
from .contracts import ElectromagneticSolverMixin
from .errors import ConfigurationError
from .progress import frequency_steps
from .validation import real_array


class ModalSolverMixin(ElectromagneticSolverMixin):
    def _frequency_case(self, frequency):
        case = copy(self)
        # Frequency is immutable for the public solver. Each working copy owns
        # its frequency-dependent operators but shares the physical mesh.
        object.__setattr__(case, "frequency", float(frequency))
        case.omega = 2 * np.pi * frequency
        case.k0 = case.omega / C_0
        case.geometry = copy(self.geometry)
        impedances = {native.name: record.material.at_frequency(frequency=frequency)
                      for record, native in self._objects.values()
                      if getattr(native, "kind", None) == "impedance"}
        case.geometry.boundaries = [
            replace(boundary, impedance=impedances[boundary.name])
            if boundary.name in impedances else boundary for boundary in self.geometry.boundaries
        ]
        if case.mesh_data is not None:
            case._reassemble_frequency()
        return case

    def sweep(self, *, frequencies, max_refinements=2, adaptive_tolerance=0.05, **solve_options):
        """Yield mode sets in increasing frequency order, adapting only at the highest frequency."""
        fs = real_array(frequencies, "frequencies")
        if fs.ndim != 1 or not len(fs) or np.any(fs <= 0) or np.any(np.diff(fs) <= 0):
            raise ConfigurationError("Sweep frequencies must be positive and strictly increasing.")
        high = self if self.frequency == fs[-1] else self._frequency_case(fs[-1])
        for frequency in frequency_steps(fs):
            if frequency == high.frequency:
                last = high.solve(max_refinements=max_refinements, adaptive_tolerance=adaptive_tolerance,
                                  **solve_options)
            else:
                case = high._frequency_case(frequency)
                modes = case.solve(max_refinements=0, adaptive_tolerance=adaptive_tolerance, **solve_options)
                yield case._finish_result(replace(modes, metadata={
                    **modes.metadata, "mesh_adaptation_frequency": high.frequency,
                }))
        yield high._finish_result(replace(last, metadata={
            **last.metadata, "mesh_adaptation_frequency": high.frequency,
        }))
