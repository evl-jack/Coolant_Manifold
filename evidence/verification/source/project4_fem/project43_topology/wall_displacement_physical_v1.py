"""Independently evaluate the converged wall-displacement v1 endpoint."""

import argparse
import json
import os
from pathlib import Path

import numpy as np

from .density import DesignGrid
from .full_density_continuation import FITTED_SPACINGS_MM, fitted_screen, summarize_fitted_replays
from .gate2 import write_report
from .gate2_restoration import measured
from .reconstruction import solve_fitted, write_reconstruction
from .wall_displacement_density_v1 import WallDisplacementDensityV1
from .wall_displacement_geometry_v1 import audit_wall_displacement, reconstruct_clipped


FLOW_BAND = (.10, .90)


def run(output, optimization, mesh_level=None, level_records=None):
    summary = json.loads((optimization / "summary.json").read_text(encoding="utf-8"))
    candidate = summary["results"].get("final_candidate") if summary["completed"] else None
    if candidate is None:
        return {"attempted": False, "reason": "Optimization has no converged final-stage candidate"}
    source = optimization / candidate / "checkpoint.json"
    checkpoint = json.loads(source.read_text(encoding="utf-8"))["results"]
    grid = DesignGrid(.25)
    pipeline = WallDisplacementDensityV1(grid)
    projected = np.asarray(checkpoint["projected_density"], dtype=float)
    geometry = reconstruct_clipped(grid, projected, pipeline.free)
    geo = write_reconstruction(geometry, output / "geometry.json")
    geometry_hash = geo["geometry_sha256"]
    if mesh_level is not None:
        spacing = FITTED_SPACINGS_MM[mesh_level]
        if geometry.port_component is None:
            return {"geometry_sha256": geometry_hash, "mesh_level": mesh_level,
                    "response": None, "reason": "inlet and both outlets are not connected"}
        (response, _, _), cost = measured(lambda: solve_fitted(geometry, spacing))
        return {"geometry_sha256": geometry_hash, "mesh_level": mesh_level,
                "response": dict(response, cost=cost, in_physical_flow_band=
                                 FLOW_BAND[0] <= response["side_share"] <= FLOW_BAND[1])}
    audit = audit_wall_displacement(grid, projected, pipeline.free, geometry)
    if level_records is None:
        mesh = fitted_screen(geometry, FLOW_BAND, area_cap_mm2=None)
    elif geometry.port_component is None:
        mesh = {"attempted": False, "reason": "inlet and both outlets are not connected"}
    else:
        records = [json.loads(path.read_text(encoding="utf-8"))["results"] for path in level_records]
        replays = []
        for index, record in enumerate(records):
            if (record["geometry_sha256"] != geometry_hash or record["mesh_level"] != index
                    or record["response"]["spacing_mm"] != FITTED_SPACINGS_MM[index]):
                raise ValueError("Fitted level records do not match the endpoint and prescribed mesh order")
            replays.append(record["response"])
        mesh = summarize_fitted_replays(replays, FLOW_BAND)
    power = checkpoint["power_w_per_m"]
    side = checkpoint["side_fraction"]
    result = {"attempted": True, "checkpoint": str(source), "geometry": geo, "wall_geometry_audit": audit,
              "immersed_power_W_per_m": power, "immersed_side_share": side, "fitted": mesh}
    if mesh.get("completed"):
        final = mesh["levels"][-1]
        result.update(physical_fluid_area_mm2=geometry.area_mm2, side_share=final["side_share"],
                      mechanical_power_W_per_m=final["mechanical_power_W_per_m"],
                      fitted_minus_immersed_power_W_per_m=final["mechanical_power_W_per_m"] - power,
                      fitted_minus_immersed_side_pp=100 * (final["side_share"] - side))
    result["physical_screen_passed"] = bool(mesh.get("fitted_diagnostic_passed")
                                            and audit["exact_clipping_passed"] and audit["connected_ports"]
                                            and audit["hole_count"] == 0 and audit["closed_pocket_count"] == 0)
    result["smoothness_passed"] = bool(audit["turning_angle_passed"])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--optimization", type=Path, required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--mesh-level", type=int, choices=range(len(FITTED_SPACINGS_MM)))
    modes.add_argument("--level-records", type=Path, nargs=len(FITTED_SPACINGS_MM))
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical verification on Sol")
    args.output.mkdir(parents=True)
    report = "fitted_mesh_level" if args.mesh_level is not None else "physical_audit"
    if not write_report(args.output, report,
                        lambda: run(args.output, args.optimization, args.mesh_level, args.level_records)):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
