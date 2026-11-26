"""Run the topology verification fixtures on ASU Sol and archive their evidence."""

import argparse
from dataclasses import asdict
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import traceback
from types import SimpleNamespace

import numpy as np

from ..fem import TaylorHood
from ..junction import boundary_data
from ..project43_baseline import Geometry, fitted_mesh
from .density import DesignGrid, DensityPipeline
from .domain import build_domain
from .evaluator import DesignTransfer, Evaluator, REYNOLDS, Stage, baseline_reference, response_values
from .kkt_audit import audit_kkt
from .optimize import optimize
from .verify_density import run_density_verification
from .verify_flow import run_flow_verification


def environment():
    versions = {}
    for package in ("numpy", "scipy", "nlopt", "shapely", "triangle"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return {"python": platform.python_version(), "host": platform.node(), "packages": versions}


def _continuation(fem, fixed, prescribed, resistance):
    state = None
    load = fem.load()
    try:
        for inertia in (0.0, REYNOLDS / 2, REYNOLDS):
            fem.inertia = inertia
            state = fem.solve(fixed, prescribed, load, resistance=resistance,
                              initial=None if state is None else state.values)
    finally:
        fem.inertia = REYNOLDS
    return state


def _baseline_record(domain, state):
    flux, imbalance, power, fraction, _ = response_values(domain, state.values)
    return dict(triangles=len(domain.fem.triangles), dofs=domain.fem.size,
                flux=flux.tolist(), mass_imbalance=imbalance, power_w_per_m=power,
                side_fraction=fraction, nonlinear_residual=state.residual_history[-1])


def baseline_checks(spacing_mm=0.5, analysis_spacing_mm=0.25, alpha_max=1e5):
    """Compare three baseline representations at the same nominal analysis spacing."""
    grid = DesignGrid(spacing_mm)
    pipeline = DensityPipeline(grid)
    domain = build_domain(grid, analysis_spacing_mm)
    stage = Stage(q=0.01, beta=0, alpha_max=alpha_max)
    density = pipeline.evaluate(grid.base_raw, **asdict(stage))
    _, filtered = baseline_reference(domain, pipeline, alpha_max)
    transfer = DesignTransfer(domain.cell_index, grid.size)
    binary = _continuation(domain.fem, domain.fixed, domain.prescribed,
                           transfer.forward(alpha_max * (1 - grid.base_raw)))
    fitted_fem = TaylorHood(*fitted_mesh(Geometry(), analysis_spacing_mm))
    fixed, prescribed, ports, flux = boundary_data(fitted_fem)
    fitted = _continuation(fitted_fem, fixed, prescribed, 0.0)
    fitted_domain = SimpleNamespace(fem=fitted_fem, ports=ports, flux=flux)

    x, y = grid.centers_mm.T
    half = spacing_mm / 2
    seam_masks = {
        "inlet_side": np.isclose(x, 6 + half) & (y < 4),
        "main_outlet": np.isclose(x, 16 - half) & (y < 4),
        "side_outlet": np.isclose(y, 6 - half) & (x > 9.5) & (x < 10.5),
    }
    seams = {}
    for name, mask in seam_masks.items():
        if not np.any(mask):
            raise RuntimeError(f"No design cells sampled at seam {name}")
        seams[name] = dict(cell_centers_mm=grid.centers_mm[mask].tolist(),
                           raw=grid.base_raw[mask].tolist(), filtered=density.filtered[mask].tolist(),
                           projected=density.projected[mask].tolist(),
                           resistance=density.resistance[mask].tolist(),
                           fixed_fluid_projected=1.0, fixed_fluid_resistance=0.0)
    return dict(stage=asdict(stage), design_spacing_mm=spacing_mm,
                analysis_spacing_mm=analysis_spacing_mm, filtered_area_mm2=density.area_mm2,
                fitted=_baseline_record(fitted_domain, fitted),
                binary=_baseline_record(domain, binary), filtered=_baseline_record(domain, filtered),
                seams=seams)


def _optimizer_result(evaluator, outcome):
    final = evaluator.evaluate(outcome.terminal_design) if outcome.failed_design is None else None
    witness = None if final is None else audit_kkt(
        outcome.terminal_design, final.gradient, final.constraints, final.jacobian)
    return {
        "status": outcome.status, "message": outcome.message,
        "terminal_design": outcome.terminal_design.tolist(),
        "terminal_source": outcome.terminal_source,
        "best_feasible_design": None if outcome.best_feasible_design is None
        else outcome.best_feasible_design.tolist(),
        "failed_design": None if outcome.failed_design is None else outcome.failed_design.tolist(),
        "objective_calls": outcome.objective_calls,
        "evaluator_requests": outcome.evaluator_requests,
        "unique_evaluations": outcome.unique_evaluations,
        "iterations": outcome.iterations,
        "objective": None if final is None else final.objective,
        "constraints": None if final is None else final.constraints.tolist(),
        "kkt": None if witness is None else dict(
            verified=witness.verified, multipliers=witness.multipliers.tolist(),
            primal=witness.primal_violation, dual=witness.dual_violation,
            projected_stationarity=witness.stationarity_inf,
            complementarity=witness.complementarity_inf),
    }


def optimizer_smoke(algorithm):
    class KnownProblem:
        def evaluate(self, design):
            delta = design - 0.8
            total = design.sum()
            return SimpleNamespace(objective=float(delta @ delta), gradient=2 * delta,
                                   constraints=np.array((0.6 - total, total - 1, -1.0)),
                                   jacobian=np.array(((-1., -1.), (1., 1.), (0., 0.))))

    evaluator = KnownProblem()
    outcome = optimize(evaluator, np.zeros(2), algorithm=algorithm, max_evaluations=300)
    report = _optimizer_result(evaluator, outcome)
    error = float(np.max(np.abs(outcome.terminal_design - 0.5)))
    report.update(known_design=[0.5, 0.5], known_objective=0.18, design_error=error,
                  passed=error <= 1e-4 and report["kkt"] is not None and report["kkt"]["verified"])
    return report


def topology_pilot(algorithm, max_evaluations):
    """Matched 240-cell pilot on a verification fixture, with no production choice."""
    grid = DesignGrid(0.5)
    pipeline = DensityPipeline(grid)
    domain = build_domain(grid, 0.5)
    stage = Stage(0.01, 0.0, 1e5)
    power, initial_state = baseline_reference(domain, pipeline, stage.alpha_max)
    evaluator = Evaluator(domain, pipeline, stage, power, initial_state=initial_state)
    initial = evaluator.evaluate(grid.base_raw)
    outcome = optimize(evaluator, grid.base_raw, algorithm=algorithm, max_evaluations=max_evaluations)
    report = _optimizer_result(evaluator, outcome)
    report.update(initial_constraints=initial.constraints.tolist(), initial_side_fraction=initial.side_fraction,
                  forward_attempts=evaluator.forward_attempts, forward_solves=evaluator.forward_solves,
                  adjoint_rhs=evaluator.adjoint_rhs, cells=grid.size, triangles=len(domain.fem.triangles),
                  stage=asdict(stage), power_reference_w_per_m=power, max_evaluations=max_evaluations,
                  evaluation_history=[dict(design=item.design.tolist(), objective=item.objective,
                                           inequalities=item.inequalities.tolist(), feasible=item.feasible)
                                      for item in outcome.evaluation_history])
    return report


def reconstruction_checks(output):
    test_file = Path(__file__).with_name("test_reconstruction.py")
    run = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", str(test_file)],
                         capture_output=True, text=True, check=False)
    (output / "reconstruction-fixtures.log").write_text(run.stdout + run.stderr, encoding="utf-8")
    return {"passed": run.returncode == 0, "exit_code": run.returncode, "summary": run.stdout[-3000:]}


def write_report(output, name, function):
    def numpy_value(value):
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
        raise TypeError(f"Unsupported report value: {type(value).__name__}")

    try:
        result = {"completed": True, "results": function()}
        formatted = json.dumps(result, allow_nan=False, indent=2, default=numpy_value)
    except Exception:
        result = {"completed": False, "error": traceback.format_exc()}
        formatted = json.dumps(result, indent=2)
    target = output / f"{name}.json"
    target.write_text(formatted + "\n", encoding="utf-8")
    print(f"{name}: {'completed' if result['completed'] else 'failed'}", flush=True)
    if not result["completed"]:
        print(result["error"], flush=True)
    return result["completed"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--flow-only", action="store_true", help="Replay the discrete flow verification only")
    modes.add_argument("--optimizers-only", action="store_true", help="Replay the optimizer checks only")
    parser.add_argument("--pilot-evaluations", type=int, default=80)
    arguments = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical verification as a Sol Slurm job")
    arguments.output.mkdir(parents=True)
    write_report(arguments.output, "environment", environment)
    if not arguments.optimizers_only:
        if not arguments.flow_only:
            write_report(arguments.output, "density", run_density_verification)
            write_report(arguments.output, "baselines", baseline_checks)
        write_report(arguments.output, "flow", run_flow_verification)
        if not arguments.flow_only:
            write_report(arguments.output, "reconstruction", lambda: reconstruction_checks(arguments.output))
    if not arguments.flow_only:
        for algorithm in ("slsqp", "mma"):
            write_report(arguments.output, f"optimizer_smoke_{algorithm}", lambda: optimizer_smoke(algorithm))
            write_report(arguments.output, f"topology_pilot_{algorithm}",
                         lambda: topology_pilot(algorithm, arguments.pilot_evaluations))


if __name__ == "__main__":
    main()
