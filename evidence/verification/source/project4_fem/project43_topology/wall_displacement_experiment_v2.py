"""Execute the single Sol experiment with a reconstructed wall-angle constraint."""

import argparse
from dataclasses import asdict, replace
from hashlib import sha256
import json
import os
from pathlib import Path
import traceback

import numpy as np
from shapely.geometry import LineString, Point, box
from shapely.ops import unary_union

from .density import DesignGrid
from .domain import build_domain
from .evaluator import REYNOLDS, Stage
from .full_density_continuation import fitted_screen
from .gate2 import write_report
from .gate2_restoration import FEASIBILITY_TOLERANCE
from .optimize import optimize
from .relaxed_flow import save
from .restoration import restore
from .wall_displacement_fem_v1 import WallDisplacementEvaluatorV1
from .wall_displacement_optimize_v1 import WorkLog, record_design, sharp_reference
from .wall_displacement_smooth_v2 import (
    GeometricConstraint, LIMIT, ProtectedDensity, SAMPLES, SCALE, SECT, STEP, TAU, arc_table,
    sample_angles, wall_chains,
)
from .wall_displacement_verify_v1 import verify


STAGE = Stage(.10, 8., 1e5)
GEOMETRIC_FD_LIMIT = 1e-4


def setup(checkpoint):
    payload = checkpoint.read_bytes()
    record = json.loads(payload)["results"]
    if (record["stage"]["q"], record["stage"]["beta"]) != (.10, 4.) or record.get("pde", "navier-stokes") != "navier-stokes":
        raise ValueError("The starting archive must be the nonlinear q=0.10, beta=4 checkpoint")
    grid = DesignGrid(.25)
    pipeline = ProtectedDensity(grid)
    prior = np.asarray(record["raw_density"], dtype=float)
    design = prior[pipeline.free].copy()
    if grid.size != 960 or not np.all((design >= 0) & (design <= 1)):
        raise ValueError("Invalid archived raw design")
    return pipeline, design, sha256(payload).hexdigest()


def independent_value(chains):
    residuals = []
    for points in chains:
        line = LineString(points)
        for s in np.linspace(SECT, line.length - SECT, SAMPLES):
            a, b, c = (np.asarray(line.interpolate(float(s + offset)).coords[0]) for offset in (-SECT, 0., SECT))
            u, v = b - a, c - b
            residuals.append(LIMIT - np.clip(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)), -1., 1.))
    peak = max(residuals)
    return float((peak + TAU * np.log(np.exp((np.array(residuals) - peak) / TAU).sum())) / SCALE)


def verify_geometry(output, checkpoint, common):
    pipeline, design, digest = setup(checkpoint)
    constraint = GeometricConstraint(pipeline, STAGE)
    initial = {"editable_cells": pipeline.size, "checkpoint_sha256": digest, "step": STEP,
               "samples_per_chain": SAMPLES, "tau": TAU, "tolerance": GEOMETRIC_FD_LIMIT}
    try:
        base = constraint.sample(design)
        independent = independent_value(base[4])
        initial.update(geometric_constraint=base[0], independent_geometric_constraint=independent,
                       chain_spacing_mm=base[2], geometry=physical_geometry(constraint, base))
        # The archived constraints at the seed are useful even if a geometry stencil fails.
        reference, _ = sharp_reference(common)
        fluid = WallDisplacementEvaluatorV1(domain_for(pipeline), pipeline, STAGE, reference).evaluate(design)
        initial.update(seed_inequalities=np.r_[fluid.constraints, base[0]], seed_power_W_per_m=fluid.power_w_per_m,
                       seed_side_fraction=fluid.side_fraction)
        save(output, "seed", initial)
        coarse = constraint.gradient(design, base, STEP)
        fine = constraint.gradient(design, base, STEP / 2.)
        errors = {"max_coordinate_difference": float(np.max(np.abs(coarse - fine))),
                  "independent_value_difference": abs(base[0] - independent)}
        directions = []
        rng = np.random.default_rng(43579)
        for _ in range(3):
            direction = rng.uniform(-1., 1., size=design.size)
            # Include genuine inward moves at numerical bounds. Use a feasible
            # one-sided direction at the mixed active/interior archive.
            direction[design < STEP] = np.abs(direction[design < STEP])
            direction[design > 1. - STEP] = -np.abs(direction[design > 1. - STEP])
            direction /= np.linalg.norm(direction)
            h = STEP / 2.
            perturbed = constraint.sample(design + h * direction)
            if perturbed[1] != base[1]:
                raise ValueError("Independent directional stencil changed contour or clipping signature")
            measured = (independent_value(perturbed[4]) - independent) / h
            predicted = fine @ direction
            directions.append(dict(independent_derivative=measured, coordinate_gradient_action=predicted,
                                   absolute_difference=abs(measured - predicted)))
        passed = (errors["max_coordinate_difference"] <= GEOMETRIC_FD_LIMIT
                  and errors["independent_value_difference"] <= 1e-10
                  and all(item["absolute_difference"] <= GEOMETRIC_FD_LIMIT for item in directions))
        return dict(initial, passed=bool(passed), errors=errors, directional_checks=directions,
                    coordinate_gradient=coarse, halved_step_gradient=fine)
    except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
        return dict(initial, passed=False, error=f"{type(error).__name__}: {error}", traceback=traceback.format_exc())


def domain_for(pipeline):
    domain = build_domain(pipeline.grid, .25)
    domain.fem.inertia = REYNOLDS
    return domain


class ConstrainedEvaluator(WallDisplacementEvaluatorV1):
    def __init__(self, domain, pipeline, reference, constraint):
        super().__init__(domain, pipeline, STAGE, reference)
        self.geometric = constraint

    def evaluate(self, design):
        if self.last_design is not None and np.array_equal(design, self.last_design):
            return self.last_evaluation
        geometric = self.geometric.sample(design)
        derivative = self.geometric.gradient(design, geometric, STEP)
        fluid = super().evaluate(design)
        self.last_evaluation = replace(fluid, constraints=np.r_[fluid.constraints, geometric[0]],
                                       jacobian=np.vstack([fluid.jacobian, derivative]))
        return self.last_evaluation


def physical_geometry(constraint, sample):
    value, _, spacing, geometry, chains = sample
    vertices, _ = sample_angles(chains, verification=True)
    turn = float(np.degrees(np.arccos(vertices.min())))
    original_wall = unary_union([LineString(points) for points in wall_chains(
        replace(geometry, geometry=constraint.original, polygons=(constraint.original,), port_components=(0, 0, 0)))])
    displacement = 0.
    for chain in chains:
        arc, chain = arc_table(chain)
        positions = np.unique(np.r_[arc, np.linspace(0., arc[-1], int(np.ceil(arc[-1] / .025)) + 1)])
        for x, y in zip(*(np.interp(positions, arc, axis) for axis in chain.T)):
            displacement = max(displacement, original_wall.distance(Point(x, y)))
    polygon = geometry.port_component
    protected = unary_union([box(*bounds) for bounds in ProtectedDensity.core_boxes])
    # The complete outlet port is 1 mm wide and the two protected strips
    # certify all normal through-channel sections are at least this wide.
    protected_intact = bool(polygon.covers(protected))
    return {"geometric_constraint": value, "max_turn_degrees": turn, "sample_spacing_mm": spacing,
            "connected_ports": True, "protected_core_covered": protected_intact,
            "minimum_normal_channel_width_mm": 1. if protected_intact else None,
            "width_basis": "Exact protected strips and 1 mm outlet, normal to each channel axis",
            "fluid_area_mm2": geometry.area_mm2,
            "added_fluid_area_mm2": polygon.difference(constraint.original).area,
            "removed_fluid_area_mm2": constraint.original.difference(polygon).area,
            "max_wall_displacement_sampled_mm": displacement,
            "max_wall_displacement_upper_bound_mm": displacement + .0125,
            "geometric_screen_passed": bool(value <= FEASIBILITY_TOLERANCE and turn <= 45. + 1e-10 and protected_intact)}


def run(output, checkpoint, common):
    pipeline, design, digest = setup(checkpoint)
    reference, _ = sharp_reference(common)
    constraint = GeometricConstraint(pipeline, STAGE)
    domain = domain_for(pipeline)
    evaluator = ConstrainedEvaluator(domain, pipeline, reference, constraint)
    work = WorkLog()
    save(output, "setup", {"seed_archive": str(checkpoint), "seed_sha256": digest, "free": pipeline.free,
                          "sharp_power_W_per_m": reference, "stage": asdict(STAGE), "tau": TAU,
                          "normalization_W_per_m": reference, "flow_band": (.10, .90), "flow_constraint_scale": .001})
    attempts = []
    chosen = None
    for index, margin in enumerate((.05, 0.)):
        result = work.measure("restoration", evaluator, lambda: restore(evaluator, design, margin=margin, max_evaluations=250))
        attempts.append(asdict(result))
        save(output, f"restoration_{index}", attempts[-1])
        if result.accepted_design is not None:
            chosen = np.asarray(result.accepted_design, dtype=float)
            break
    if chosen is None:
        records = [row for attempt in attempts for row in attempt["trial_history"]
                   if row.get("inequalities") is not None]
        best = min(records, key=lambda row: max(row["inequalities"]), default=None)
        return {"passed": False, "stop_reason": "Both Phase I attempts failed", "best_evaluated": best,
                "work": work.totals()}
    result = work.measure("mma", evaluator, lambda: optimize(evaluator, chosen, algorithm="mma", max_evaluations=500))
    save(output, "mma", asdict(result))
    terminal = np.asarray(result.terminal_design, dtype=float)
    folder = output / "terminal"
    folder.mkdir()
    replay = ConstrainedEvaluator(domain, pipeline, reference, constraint)
    record = record_design(replay, terminal, folder, stationarity=True, geometry_constraint=constraint)
    sample = constraint.sample(terminal)
    geometric = physical_geometry(constraint, sample)
    save(folder, "checkpoint", record)
    native = result.optimizer_success and result.status in (1, 2, 3, 4)
    fitted = None
    if geometric["geometric_screen_passed"]:
        fitted = fitted_screen(sample[3], (.10, .90), area_cap_mm2=None)
    return {"passed": bool(native and record["kkt"]["verified"] and fitted and fitted["fitted_diagnostic_passed"]),
            "mma_native_converged": native, "mma_status": result.status, "kkt_verified": record["kkt"]["verified"],
            "physical_geometry": geometric, "fitted": fitted, "work": work.totals()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("fluid", "geometry", "optimize"))
    parser.add_argument("output", type=Path)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--common", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run numerical work in a Sol Slurm job")
    args.output.mkdir(parents=True, exist_ok=False)
    if args.mode == "fluid":
        function = lambda: verify(args.common / "sharp_fitted_baseline.json", protected_v2=True)
    elif args.mode == "geometry":
        function = lambda: verify_geometry(args.output, args.checkpoint, args.common)
    else:
        function = lambda: run(args.output, args.checkpoint, args.common)
    completed = write_report(args.output, "summary", function)
    passed = completed and json.loads((args.output / "summary.json").read_text())["results"]["passed"]
    raise SystemExit(0 if passed else 2)


if __name__ == "__main__":
    main()
