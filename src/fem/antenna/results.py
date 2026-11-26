"""Modal power waves and far-field metrics for a coherent antenna feed."""

from dataclasses import dataclass

import numpy as np

from fem.common.constants import ETA_0
from fem.common.contracts import MeshSnapshot, ResultMixin
from .farfield import HuygensSurface


@dataclass
class RadiationPattern:
    electric: np.ndarray
    intensity: np.ndarray
    directivity: np.ndarray
    gain: np.ndarray
    realized_gain: np.ndarray


@dataclass
class AntennaResult(ResultMixin):
    frequency: float
    channels: tuple
    incident: np.ndarray
    modal_s_parameters: dict
    huygens: HuygensSurface
    radiated_power: float
    _mesh_snapshot: MeshSnapshot
    electric_field: np.ndarray
    port_modes: dict
    current_facets: np.ndarray
    surface_current: np.ndarray
    accepted_power: float
    metadata: dict
    magnetic_field: np.ndarray | None = None
    surface_current_vertices: np.ndarray | None = None
    dimension = 3

    @property
    def s_parameters(self):
        """Selected-mode S matrix. Rows receive, columns excite independently."""
        return np.stack([self.modal_s_parameters[key] for key in self.channels])

    @property
    def outgoing(self):
        """Outgoing coherent power waves across all retained propagating modes."""
        return {key: row @ self.incident for key, row in self.modal_s_parameters.items()}

    @property
    def active_s_parameters(self):
        """Active reflection b_i/a_i in each selected mode, NaN at unpowered ports."""
        b = self.s_parameters @ self.incident
        return np.divide(
            b, self.incident, out=np.full_like(b, complex(np.nan, np.nan)), where=self.incident != 0
        )

    @property
    def incident_power(self):
        return float(np.vdot(self.incident, self.incident).real)

    @property
    def port_voltage(self):
        """Peak complex terminal voltage for each discrete feed, in volts."""
        z = self.metadata["reference_impedances"]
        outgoing = self.outgoing
        return {
            name: np.sqrt(2 * z[name]) * (a + outgoing[name, mode])
            for (name, mode), a in zip(self.channels, self.incident)
            if name in z
        }

    @property
    def port_current(self):
        """Peak complex current delivered into each discrete feed, in amperes."""
        z = self.metadata["reference_impedances"]
        outgoing = self.outgoing
        return {
            name: np.sqrt(2 / z[name]) * (a - outgoing[name, mode])
            for (name, mode), a in zip(self.channels, self.incident)
            if name in z
        }

    def show_port_modes(self, *, port, block=True):
        """Inspect the numerical port modes in the 2D waveguide viewer."""
        return self.port_modes[port].show(block=block)

    def far_field(self, directions):
        """Coherent far electric amplitude (V), intensity (W/sr), and linear power ratios."""
        electric = self.huygens.amplitude(directions)
        intensity = np.sum(abs(electric) ** 2, axis=1) / (2 * ETA_0)

        def normalized(power):
            return 4 * np.pi * intensity / power if power > 0 else np.full_like(intensity, np.nan)

        return RadiationPattern(
            electric,
            intensity,
            normalized(self.radiated_power),
            normalized(self.accepted_power),
            normalized(self.incident_power),
        )

    def plot(self, *, component="S-parameters", quantity=None, phase_deg=0.0):
        from .visualization import plot_result

        return plot_result(self, component, quantity, phase_deg)


@dataclass
class AntennaSweep(ResultMixin):
    cases: tuple
    dimension = 3

    def __len__(self):
        return len(self.cases)

    def __getitem__(self, index):
        return self.cases[index]

    @property
    def mesh_data(self):
        return self.cases[0].mesh_data

    @property
    def metadata(self):
        return dict(self.cases[-1].metadata)

    @property
    def frequencies(self):
        return np.array([r.frequency for r in self.cases])

    @property
    def s_parameters(self):
        return np.stack([r.s_parameters for r in self.cases])

    @property
    def active_s_parameters(self):
        return np.stack([r.active_s_parameters for r in self.cases])

    def plot(self, *, component="S-parameters", case=0, quantity=None, phase_deg=0.0):
        if component != "S-parameters":
            return self.cases[case].plot(
                component=component, quantity=quantity, phase_deg=phase_deg
            )
        from .visualization import plot_sweep

        return plot_sweep(self, quantity)
