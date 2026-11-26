"""Isolate the automatic quadratic target change for scattering and antenna meshes."""

import argparse
import csv
from contextlib import nullcontext
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

import fem
import numpy as np
from matplotlib.figure import Figure

from amr_observables import inclusion, multiport


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    # Only replace the legacy background size. Both runs use the public
    # automatic path, including its material-aware CAD size field.
    for name, build, cls in [("cylinder", inclusion, fem.ScatteringSolver2D),
                             ("two-port", multiport, fem.AntennaSolver3D)]:
        options = {"element_order": 2} if name == "two-port" else {}
        for denominator in (6, 4):
            control = (patch.object(cls, "_wavelength_size", staticmethod(lambda wavelength, order: wavelength / 6))
                       if denominator == 6 else nullcontext())
            start = perf_counter()
            with control:
                result, observable, measured, exact = build(None, dict(backend="cpu", max_refinements=0), **options)
            elapsed = perf_counter() - start
            if not np.isclose(result.mesh_data.info["size"], 1/denominator):
                raise RuntimeError("The solver used a different background target than the requested comparison.")
            record = result.solve_info["adaptive_history"][0]
            rows.append(dict(case=name, cells_per_wavelength=denominator, observable=observable,
                             relative_error=np.linalg.norm(measured-exact)/np.linalg.norm(exact),
                             physical_cells=record["elements"], field_dofs=record["field_dofs"], elapsed_s=elapsed))
            print(rows[-1], flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    fig = Figure(figsize=(7.6, 3.6), layout="constrained")
    for ax, key, label in zip(fig.subplots(1, 2), ("relative_error", "field_dofs"),
                              ("Relative observable error (%)", "Field DOFs")):
        for name in ("cylinder", "two-port"):
            case = [r for r in rows if r["case"] == name]
            ax.semilogy(["λ/6", "λ/4"], [r[key] * (100 if key == "relative_error" else 1) for r in case],
                        marker="o", label=name)
        ax.set(ylabel=label, xlabel="Automatic size target")
        ax.grid(alpha=.25)
        ax.legend()
    fig.savefig(args.output / "comparison.png", dpi=150)


if __name__ == "__main__":
    main()
