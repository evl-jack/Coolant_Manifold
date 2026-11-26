"""Follow a fixed binary fixture with an interface-resolved diagnostic mesh."""

import argparse
import json
import os
from pathlib import Path
from time import perf_counter

import numpy as np

from .density import DesignGrid
from .domain import build_domain
from .evaluator import DesignTransfer
from .gate2 import _baseline_record, _continuation, write_report
from .gate3_screen import obstruction, PENALTIES
from .leakage import binary_leakage
from .local_refinement import build_local_binary_domain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("sharp", "obstructed"))
    parser.add_argument("alpha", type=float, choices=PENALTIES)
    parser.add_argument("fitted_reference", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--base-spacing", type=float, default=0.125)
    parser.add_argument("--max-triangles", type=int, default=170_000)
    parser.add_argument("--mesh-only", action="store_true")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical screen as a Sol Slurm job")
    args.output.mkdir(parents=True)
    grid = DesignGrid(0.25)
    raw = grid.base_raw if args.case == "sharp" else obstruction(grid)
    base = build_domain(grid, args.base_spacing)
    domain = None

    def mesh_report():
        nonlocal domain
        start = perf_counter()
        domain, details = build_local_binary_domain(base, raw, args.alpha,
                                                    max_triangles=args.max_triangles)
        return dict(case=args.case, alpha_max=args.alpha, base_spacing_mm=args.base_spacing,
                    sol_job_id=os.environ["SLURM_JOB_ID"], elapsed_seconds=perf_counter() - start,
                    study=details)

    if not write_report(args.output, "mesh", mesh_report):
        raise SystemExit(1)
    if args.mesh_only:
        return

    def solve_report():
        fitted = json.loads(args.fitted_reference.read_text(encoding="utf-8"))["fixtures"][args.case]
        start = perf_counter()
        drag = DesignTransfer(domain.cell_index, grid.size).forward(args.alpha * (1 - raw))
        state = _continuation(domain.fem, domain.fixed, domain.prescribed, drag)
        immersed = _baseline_record(domain, state)
        leakage = binary_leakage(domain, state.values, raw)
        split_error = 100 * abs(immersed["side_fraction"] - fitted["side_fraction"])
        power_error = abs(immersed["power_w_per_m"] / fitted["power_w_per_m"] - 1)
        np.testing.assert_allclose(41 + grid.cell_area_mm2 * raw.sum(), 83., rtol=1e-12)
        return dict(case=args.case, alpha_max=args.alpha, immersed=immersed, fitted=fitted,
                    leakage=leakage, flow_error_pp=split_error, power_error_relative=power_error,
                    per_fixture_thresholds_passed=bool(split_error <= 0.02 and power_error <= 0.005
                                                       and leakage["passed"]),
                    elapsed_seconds=perf_counter() - start, scope="fixed binary fixture; no production mesh approval")

    if not write_report(args.output, "solution", solve_report):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
