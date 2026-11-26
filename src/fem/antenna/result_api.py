"""Antenna result archives and shared viewer integration."""

from pathlib import Path
import h5py

from fem.common.contracts import MeshSnapshot
from fem.common.errors import PersistenceError
from fem.common.persistence import (
    atomic_h5,
    read_value,
    validate_envelope,
    write_envelope,
    write_value,
)
from fem.waveguide import ModeSet, Mode
from fem.waveguide.result_api import _sampled_fields
from .farfield import HuygensSurface
from .results import AntennaResult, AntennaSweep


def save_result(result, path):
    with atomic_h5(path) as archive:
        write_envelope(
            archive,
            family="antenna",
            kind="sweep" if isinstance(result, AntennaSweep) else "single",
            dimension=3,
            representation="power-waves-and-huygens-currents",
        )
        write_value(archive, "result", result)
    return Path(path)


def load_result(path):
    registry = {c.__name__: c for c in (MeshSnapshot, AntennaResult, AntennaSweep, HuygensSurface)}

    def restore_modes(**values):
        modes = ModeSet(**values)
        object.__setattr__(modes, "_mesh_snapshot", modes.metadata["port_mesh"])
        return modes

    registry.update(ModeSet=restore_modes, Mode=Mode, SampledFields=_sampled_fields)
    try:
        with h5py.File(path) as archive:
            validate_envelope(archive, family="antenna")
            result = read_value(archive["result"], registry)
            if not isinstance(result, (AntennaResult, AntennaSweep)):
                raise PersistenceError("Invalid antenna result.")
            return result
    except (KeyError, OSError, ValueError, TypeError) as exc:
        raise PersistenceError(f"Cannot load antenna archive {path}: {exc}") from exc


def show_result(result, *, block=True):
    from fem.common.viewer import show

    return show(result, block=block)
