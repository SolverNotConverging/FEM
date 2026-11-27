"""Console progress for frequency sweeps and mesh refinement."""
import sys
from itertools import chain

from tqdm.auto import tqdm


def frequency_steps(frequencies):
    """Visit the adaptation frequency first; advance after each completed solve."""
    with tqdm(total=len(frequencies), desc="Frequency sweep", unit="freq") as bar:
        for frequency in chain(frequencies[-1:], frequencies[:-1]):
            bar.set_postfix_str(f"{frequency:g} Hz")
            yield float(frequency)
            bar.update()


def report_refinement(owner, step, limit, cells, residual, tolerance, status, field_dofs=None,
                      *, criterion=None):
    if limit:
        dofs = f"  field_dofs={field_dofs}" if field_dofs is not None else ""
        target = criterion if criterion is not None else f"tolerance={tolerance:g}"
        tqdm.write(f"{type(owner).__name__}  refinement {step}/{limit}  "
                   f"elements={cells}{dofs}  residual={residual:.4g}  "
                   f"{target}  {status}", file=sys.stderr)
