"""Screen unfiltered binary Brinkman fixtures against independent fitted solves."""

import argparse
import json
import os
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace

import numpy as np

from ..fem import TaylorHood
from ..junction import boundary_data
from ..project43_baseline import Geometry, fitted_mesh as sharp_mesh
from .density import DesignGrid
from .domain import build_domain
from .evaluator import DesignTransfer, VISCOSITY
from .gate2 import _baseline_record, _continuation
from .leakage import binary_leakage
from .reconstruction import fitted_mesh


PENALTIES = (1e5, 1e6, 1e7)
ANALYSIS_SPACINGS_MM = (0.25, 0.125, 0.0625)


def obstruction(grid):
    """Move 1 mm² from a bottom-wall throttle into a shoulder next to the branch."""
    x, y = grid.centers_mm.T
    raw = grid.base_raw.copy()
    raw[(x > 13) & (x < 14) & (y < 1)] = 0
    raw[(x > 10.5) & (x < 11.5) & (y > 4) & (y < 5)] = 1
    return raw


def fitted_reference(case, spacing_mm):
    if case == "sharp":
        mesh = sharp_mesh(Geometry(), spacing_mm)
    else:
        from shapely.geometry import box

        polygon = box(0, 0, 20, 4).union(box(9.5, 4, 10.5, 7))
        polygon = polygon.difference(box(13, 0, 14, 1)).union(box(10.5, 4, 11.5, 5))
        np.testing.assert_allclose(polygon.area, 83.0, rtol=1e-12)
        mesh = fitted_mesh(polygon, spacing_mm)
    fem = TaylorHood(*mesh)
    fixed, prescribed, ports, flux = boundary_data(fem)
    state = _continuation(fem, fixed, prescribed, 0.)
    return _baseline_record(SimpleNamespace(fem=fem, ports=ports, flux=flux), state)


def screen(alpha_max, spacing_mm, fitted_reports):
    grid = DesignGrid(0.25)
    domain = build_domain(grid, spacing_mm)
    transfer = DesignTransfer(domain.cell_index, grid.size)
    fixture_reports = {}
    for case, raw in (("sharp", grid.base_raw), ("obstructed", obstruction(grid))):
        start = perf_counter()
        np.testing.assert_allclose(41 + grid.cell_area_mm2 * raw.sum(), 83., rtol=1e-12)
        fitted = fitted_reports[case]
        resistance = transfer.forward(alpha_max * (1 - raw))
        state = _continuation(domain.fem, domain.fixed, domain.prescribed, resistance)
        immersed = _baseline_record(domain, state)
        leakage = binary_leakage(domain, state.values, raw)
        flow_error_pp = 100 * abs(immersed["side_fraction"] - fitted["side_fraction"])
        power_error_relative = abs(immersed["power_w_per_m"] / fitted["power_w_per_m"] - 1)
        fixture_reports[case] = dict(
            fitted=fitted, immersed=immersed, leakage=leakage,
            flow_error_pp=flow_error_pp, power_error_relative=power_error_relative,
            per_mesh_thresholds_passed=bool(flow_error_pp <= 0.02 and power_error_relative <= 0.005
                                           and leakage["passed"]),
            elapsed_seconds=perf_counter() - start,
        )
    layer_mm = 1000 * np.sqrt(VISCOSITY / alpha_max)
    return dict(alpha_max=alpha_max, brinkman_length_mm=layer_mm,
                analysis_spacing_mm=spacing_mm, spacing_over_brinkman_length=spacing_mm/layer_mm,
                design_spacing_mm=grid.spacing_mm, fixtures=fixture_reports)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("alpha", type=float, choices=PENALTIES)
    parser.add_argument("spacing", type=float, choices=ANALYSIS_SPACINGS_MM)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference", type=Path, help="Archived fitted result on this mesh level")
    parser.add_argument("--fitted-only", action="store_true", help="Solve the independent fitted references only")
    arguments = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical screen as a Sol Slurm job")
    if arguments.fitted_only:
        result = dict(analysis_spacing_mm=arguments.spacing,
                      fixtures={case: fitted_reference(case, arguments.spacing) for case in ("sharp", "obstructed")})
    else:
        if arguments.reference is None:
            parser.error("Provide the independent fitted result with --reference")
        reference = json.loads(arguments.reference.read_text(encoding="utf-8"))
        if reference["analysis_spacing_mm"] != arguments.spacing:
            parser.error("The fitted reference mesh level must match the screen level")
        result = screen(arguments.alpha, arguments.spacing, reference["fixtures"])
    result["sol_job_id"] = os.environ["SLURM_JOB_ID"]
    with arguments.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
