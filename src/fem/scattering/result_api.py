"""Data-only HDF5 persistence for plane-wave solutions."""
from pathlib import Path

import h5py

from fem.common.contracts import MeshSnapshot
from fem.common.errors import PersistenceError
from fem.common.persistence import atomic_h5, read_value, validate_envelope, write_envelope, write_value
from .results import ScatteringResult

FAMILY = "plane_wave_scattering"


def save_result(result, path):
    with atomic_h5(path) as handle:
        write_envelope(handle, family=FAMILY, kind="single", dimension=result.dimension,
                       representation="plane-wave-phasors")
        write_value(handle, "result", result)
    return Path(path)


def load_result(path):
    try:
        with h5py.File(path) as handle:
            validate_envelope(handle, family=FAMILY)
            result = read_value(handle["result"], {"ScatteringResult": ScatteringResult,
                                                  "MeshSnapshot": MeshSnapshot})
            if not isinstance(result, ScatteringResult) or result.dimension != handle.attrs["dimension"]:
                raise PersistenceError("Invalid plane-wave result.")
            return result
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise PersistenceError(f"Cannot load scattering result from {path}: {exc}") from exc


def show_result(result, *, block=True):
    from fem.common.viewer import show
    return show(result, block=block)
