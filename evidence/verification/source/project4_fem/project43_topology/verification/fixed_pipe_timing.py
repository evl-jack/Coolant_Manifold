"""Measure fixed-pipe adjoint and centered difference costs on Sol."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import scipy
from scipy.sparse.linalg import splu

from ... import fem as fem_module
from ..evaluator import Stage, response_values
from ..fixed_pipe_pilot_model import SIDE_CONSTRAINT_SCALE, FixedPipeEvaluator, build_fixed_pipe


STAGE = Stage(q=.10, beta=8., alpha_max=1e5)
CASES = {
    "pilot-128": (.25, .25),
    "scale-32": (.125, .50),
    "scale-128": (.125, .25),
    "scale-512": (.125, .125),
}
THREADS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
PROBE_STEPS = (.01, .001, .0001, .00001, .000001)
RESPONSES = ("normalized_power", "side_fraction")


@dataclass(frozen=True)
class Elapsed:
    wall_seconds: float
    cpu_seconds: float

    def __add__(self, other):
        return Elapsed(self.wall_seconds + other.wall_seconds, self.cpu_seconds + other.cpu_seconds)


@dataclass
class ForwardDiagnostics:
    lu: Elapsed = Elapsed(0., 0.)
    substitution: Elapsed = Elapsed(0., 0.)
    factor_count: int = 0


def stamp():
    return Elapsed(perf_counter(), process_time())


def since(start):
    end = stamp()
    return Elapsed(end.wall_seconds - start.wall_seconds, end.cpu_seconds - start.cpu_seconds)


def box_reference(path):
    """Read the archived unchanged immersed beta=8 responses."""
    document = json.loads(path.read_text(encoding="utf-8"))
    reference = document.get("results", document)
    return float(reference["power_w_per_m"]), float(reference["side_fraction"])


def cpu_model():
    info = Path("/proc/cpuinfo")
    if info.exists():
        for line in info.read_text(encoding="utf-8").splitlines():
            key, _, value = line.partition(":")
            if key.strip() == "model name":
                return value.strip()
    return platform.processor()


def peak_rss_bytes():
    try:
        import resource
    except ImportError:
        return None
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def write_json(path, value):
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def gray_field(pipeline):
    """Sample one interior physical field at each design grid's cell centers."""
    x, y = pipeline.grid.centers_mm[pipeline.free].T
    return .28 + .46 * np.exp(-((x - 10.) / .75) ** 2) * (1. - .12 * (y - 4.)) + .05 * np.cos(np.pi * y)


def forward(evaluator, design, initial, *, differentiate=False):
    timing = ForwardDiagnostics()

    def factor(*args, **kwargs):
        start = stamp()
        lu = splu(*args, **kwargs)
        timing.lu += since(start)
        timing.factor_count += 1
        original_solve = lu.solve

        def solve(*solve_args, **solve_kwargs):
            start_solve = stamp()
            result = original_solve(*solve_args, **solve_kwargs)
            timing.substitution += since(start_solve)
            return result

        # TaylorHood.solve only uses the return value's solve method.
        return SimpleNamespace(solve=solve)

    start = stamp()
    density = evaluator.density(design)
    resistance = evaluator.transfer.forward(density.resistance)
    filter_transfer = since(start)
    solver_start = stamp()
    with patch.object(fem_module, "splu", factor):
        state = evaluator.solve(resistance, initial)
    solver_time = since(solver_start)
    response_start = stamp()
    _, _, power, fraction, state_gradients = response_values(evaluator.domain, state.values,
                                                           differentiate=differentiate)
    responses = np.array([power / evaluator.power_reference, fraction])
    if differentiate:
        state_gradients[0] /= evaluator.power_reference
    response_time = since(response_start)
    return state, responses, state_gradients, density, {
        "total": asdict(since(start)),
        "nonlinear_state": asdict(solver_time),
        "filter_and_transfer": asdict(filter_transfer),
        "port_response_evaluation": asdict(response_time),
        "newton_factorization": asdict(timing.lu),
        "newton_substitution": asdict(timing.substitution),
        "newton_iterations": len(state.residual_history) - 1,
        "newton_factorizations": timing.factor_count,
        "residual": state.residual_history[-1],
        "forward_solves": 1,
        "peak_rss_bytes": peak_rss_bytes(),
    }


def adjoint(evaluator, state, density, state_gradients):
    start = stamp()
    free = state.free
    tangent = state.tangent[free][:, free].tocsc()
    slicing = since(start)
    factor_start = stamp()
    factor = splu(tangent)
    factor_time = since(factor_start)
    fields_start = stamp()
    fem = evaluator.domain.fem
    velocity = fem.velocity_fields(state.values)[0]
    fields_time = since(fields_start)
    gradients = []
    records = []
    for label, derivative in zip(RESPONSES, state_gradients):
        response_start = stamp()
        dual = np.zeros(fem.size)
        solve_start = stamp()
        dual[free] = factor.solve(derivative[free], trans="T")
        solve_time = since(solve_start)
        reverse_start = stamp()
        dual_velocity = fem.velocity_fields(dual)[0]
        cotangent = -fem.weights * np.sum(velocity * dual_velocity, axis=2)
        gradients.append(evaluator.pipeline.pullback(density, evaluator.transfer.transpose(cotangent)))
        records.append({"response": label, "total_excluding_shared_setup": asdict(since(response_start)),
                        "transpose_substitution": asdict(solve_time), "gradient_assembly": asdict(since(reverse_start))})
    total = since(start)
    return np.asarray(gradients), {
        "total": asdict(total), "total_tangent_slicing": asdict(slicing), "sparse_factorization": asdict(factor_time),
        "shared_state_fields": asdict(fields_time), "adjoint_rhs": 2, "adjoint_factorizations": 1,
        "responses": records, "peak_rss_bytes": peak_rss_bytes(),
    }


def make_evaluator(mesh_spacing, density_spacing, reference):
    domain, pipeline = build_fixed_pipe(max_spacing_mm=mesh_spacing, density_spacing_mm=density_spacing)
    return FixedPipeEvaluator(domain, pipeline, STAGE, *reference)


_worker_model = None
_worker_initial = None
_worker_design = None


def init_worker(mesh_spacing, density_spacing, reference, initial, design):
    global _worker_model, _worker_initial, _worker_design
    _worker_model = make_evaluator(mesh_spacing, density_spacing, reference)
    _worker_initial = initial
    _worker_design = design


def difference(evaluator, initial, design, coordinate, step):
    start = stamp()
    results = {}
    for label, signed_step in (("minus", -step), ("plus", step)):
        perturbed = design.copy()
        perturbed[coordinate] += signed_step
        if not 0. < perturbed[coordinate] < 1.:
            raise ValueError("A centered benchmark step left the interior density box")
        _, response, _, _, timing = forward(evaluator, perturbed, initial)
        results[label] = {"responses": response.tolist(), "timing": timing}
    delta = (np.asarray(results["plus"]["responses"]) - np.asarray(results["minus"]["responses"])) / (2. * step)
    return {"coordinate": int(coordinate), "step": step, "gradient": delta.tolist(),
            "points": results, "task_time": asdict(since(start)), "worker_pid": os.getpid(),
            "peak_rss_bytes": peak_rss_bytes()}


def worker_difference(coordinate, step):
    return difference(_worker_model, _worker_initial, _worker_design, coordinate, step)


def summarize_fd(rows, elapsed, workers):
    keys = ("total", "nonlinear_state", "filter_and_transfer", "port_response_evaluation",
            "newton_factorization", "newton_substitution")
    points = [point["timing"] for row in rows for point in row["points"].values()]
    aggregate = {
        key: {field: sum(timing[key][field] for timing in points)
              for field in ("cpu_seconds", "wall_seconds")}
        for key in keys
    }
    task = {field: sum(row["task_time"][field] for row in rows) for field in ("cpu_seconds", "wall_seconds")}
    rss = [row["peak_rss_bytes"] for row in rows if row["peak_rss_bytes"] is not None]
    return {
        "parallel_elapsed": asdict(elapsed), "aggregate_perturbed_forward_work": aggregate,
        "aggregate_coordinate_task_work": task, "parent_peak_rss_bytes": peak_rss_bytes(),
        "max_worker_peak_rss_bytes": max(rss, default=None),
        "configured_workers": workers, "observed_worker_count": len({row["worker_pid"] for row in rows}),
        "newton_iterations": sum(point["newton_iterations"] for point in points),
        "newton_factorizations": sum(point["newton_factorizations"] for point in points),
        "forward_solves": len(points), "adjoint_rhs": 0,
    }


def errors(adjoint_gradient, finite_difference, coordinates):
    selected = adjoint_gradient[:, coordinates]
    delta = finite_difference - selected
    rows = {}
    for i, label in enumerate(RESPONSES):
        absolute = np.abs(delta[i])
        gradnorm = np.linalg.norm(selected[i])
        relnorm = np.linalg.norm(delta[i]) / gradnorm if gradnorm > 0 else None
        rows[label] = {"absolute_error_linf": float(np.max(absolute)),
                       "absolute_error_l2": float(np.linalg.norm(delta[i])),
                       "relative_error_l2": None if relnorm is None else float(relnorm),
                       "worst_coordinate": int(coordinates[np.argmax(absolute)])}
    return rows


def collect_fd(args, evaluator, density_spacing, mesh_spacing, reference, initial, design, coordinates, output):
    """Save each successful coordinate immediately, including both shared responses."""
    complete = []
    start = stamp()
    with (output / "full_fd_progress.jsonl").open("x", encoding="utf-8", buffering=1) as log:
        def record(row):
            complete.append(row)
            log.write(json.dumps(row, allow_nan=False) + "\n")
            print(f"Centered FD coordinate {row['coordinate']} completed ({len(complete)}/{len(coordinates)})",
                  flush=True)

        if args.workers == 1:
            for coordinate in coordinates:
                record(difference(evaluator, initial, design, coordinate, args.step))
        else:
            with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn"),
                                     initializer=init_worker,
                                     initargs=(mesh_spacing, density_spacing, reference, initial, design)) as pool:
                futures = [pool.submit(worker_difference, coordinate, args.step) for coordinate in coordinates]
                for future in as_completed(futures):
                    record(future.result())
    complete.sort(key=lambda row: row["coordinate"])
    return complete, summarize_fd(complete, since(start), args.workers)


def study(args):
    if args.workers < 1 or args.coordinate_limit is not None and args.coordinate_limit < 1:
        raise ValueError("Worker count and any smoke-test coordinate limit must be positive")
    if not np.isfinite(args.step) or args.step <= 0.:
        raise ValueError("The finite-difference step must be finite and positive")
    start = stamp()
    args.output.mkdir(parents=True, exist_ok=False)
    reference = box_reference(args.reference_metrics)
    mesh_spacing, density_spacing = CASES[args.case]
    evaluator = make_evaluator(mesh_spacing, density_spacing, reference)
    design = gray_field(evaluator.pipeline)
    setup = since(start)
    start = stamp()
    state, responses, state_gradients, density, forward_timing = forward(evaluator, design, None, differentiate=True)
    gradient, adjoint_timing = adjoint(evaluator, state, density, state_gradients)
    total_adjoint = since(start)

    mesh = evaluator.domain.fem
    mesh_digest = hashlib.sha256(mesh.points.tobytes() + mesh.triangles.tobytes()).hexdigest()
    header = {
        "case": args.case, "stage": asdict(STAGE), "mesh_spacing_mm": mesh_spacing,
        "density_spacing_mm": density_spacing, "design_variables": design.size, "triangles": len(mesh.triangles),
        "fem_dofs": mesh.size, "mesh_sha256": mesh_digest,
        "reference_power_w_per_m": reference[0], "reference_side_fraction": reference[1],
        "reference_metrics_sha256": hashlib.sha256(args.reference_metrics.read_bytes()).hexdigest(),
        "response_order": list(RESPONSES),
        "side_constraint_gradient_factor": -1. / SIDE_CONSTRAINT_SCALE,
        "gray_design": design.tolist(), "gray_responses": responses.tolist(), "adjoint_gradients": gradient.tolist(),
        "forward": forward_timing, "adjoint": adjoint_timing, "complete_adjoint": asdict(total_adjoint),
        "setup_and_model_assembly_excluded": asdict(setup),
        "initialization": "Cold common baseline and identical baseline-state warm start for every centered perturbation",
        "timing": "Wall and process CPU seconds; workers limited to one math-library thread",
        "memory": "Per-process peak RSS, cumulative from process launch, reported as max worker and parent separately",
        "newton_residual_tolerance": 1e-11, "newton_iteration_limit": 35,
        "hardware": {"hostname": platform.node(), "cpu": cpu_model(), "python": platform.python_version(),
                     "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform(),
                     "threads": {key: os.getenv(key) for key in THREADS},
                     "sol_job_id": os.getenv("SLURM_JOB_ID"), "slurm_cpus": os.getenv("SLURM_CPUS_PER_TASK")},
    }
    write_json(args.output / "adjoint.json", header)
    # Centered FD uses the same cold state to warm-start each perturbation.
    # Its baseline needs response values only, so charge it the separately timed
    # value-only postprocessing on that same state and count just one PDE solve.
    fd_response_start = stamp()
    response_values(evaluator.domain, state.values, differentiate=False)
    fd_response = since(fd_response_start)
    response_with_gradients = forward_timing["port_response_evaluation"]
    fd_baseline = Elapsed(
        forward_timing["total"]["wall_seconds"] - response_with_gradients["wall_seconds"] + fd_response.wall_seconds,
        forward_timing["total"]["cpu_seconds"] - response_with_gradients["cpu_seconds"] + fd_response.cpu_seconds,
    )
    coordinates = np.arange(design.size if args.coordinate_limit is None else min(args.coordinate_limit, design.size))
    rows, fd_timing = collect_fd(args, evaluator, density_spacing, mesh_spacing, reference, state.values, design,
                                 coordinates, args.output)
    fd_gradient = np.asarray([row["gradient"] for row in rows]).T
    fd_wall = fd_timing["parallel_elapsed"]["wall_seconds"] + fd_baseline.wall_seconds
    fd_work = fd_timing["aggregate_coordinate_task_work"]["wall_seconds"] + fd_baseline.wall_seconds
    fd_cpu = fd_timing["aggregate_coordinate_task_work"]["cpu_seconds"] + fd_baseline.cpu_seconds
    write_json(args.output / "benchmark.json", {
        **header, "complete_coordinate_gradient": coordinates.size == design.size,
        "fd_step": args.step, "fd_coordinates": coordinates.tolist(), "fd_gradients": fd_gradient.tolist(),
        "gradient_error": errors(gradient, fd_gradient, coordinates), "finite_difference": fd_timing,
        "including_common_baseline": {
            "fd_baseline": asdict(fd_baseline),
            "fd_forward_solve_count": 1 + fd_timing["forward_solves"], "adjoint_forward_solve_count": 1,
            "fd_elapsed_wall_seconds": fd_wall, "fd_aggregate_compute_wall_seconds": fd_work,
            "fd_aggregate_compute_cpu_seconds": fd_cpu,
            "measured_serial_fd_over_adjoint_speedup": fd_wall / total_adjoint.wall_seconds
            if args.workers == 1 and coordinates.size == design.size else None,
            "measured_elapsed_fd_over_adjoint_ratio": fd_wall / total_adjoint.wall_seconds,
            "aggregate_fd_over_adjoint_cpu_ratio": fd_cpu / total_adjoint.cpu_seconds,
        },
    })
    if args.case == "pilot-128" and not args.skip_sweep:
        # Include the largest sensitivity coordinate for each response and a reproducible
        # independent coordinate. Reuse the fixed-step FD solves when already present.
        rng = np.random.default_rng(432026)
        probes = sorted(set([*np.argmax(np.abs(gradient), axis=1).tolist(), int(rng.integers(design.size))]))
        by_coordinate = {row["coordinate"]: row for row in rows}
        sweep = []
        for step in PROBE_STEPS:
            probe_rows = [by_coordinate[index] if index in by_coordinate and step == args.step
                          else difference(evaluator, state.values, design, index, step)
                          for index in probes]
            probe_gradient = np.asarray([row["gradient"] for row in probe_rows]).T
            sweep.append({"step": step, "gradient_error": errors(gradient, probe_gradient, probes),
                          "results": probe_rows})
        write_json(args.output / "step_sweep.json", {"probe_coordinates": probes, "steps": sweep})


def plot(root, output):
    import matplotlib.pyplot as plt

    candidates = [json.loads(path.read_text(encoding="utf-8")) for path in root.rglob("benchmark.json")]
    complete = {row["case"]: row for row in candidates if row["complete_coordinate_gradient"]}
    scaling_names = ("scale-32", "scale-128", "scale-512")
    scaling = [complete[name] for name in scaling_names]
    if len({row["mesh_sha256"] for row in scaling}) != 1:
        raise RuntimeError("Scaling cases must use the same FEM mesh")
    if len({row["hardware"]["cpu"] for row in scaling}) != 1:
        raise RuntimeError("Scaling cases use different CPU models")
    if len({json.dumps(row["hardware"]["threads"], sort_keys=True) for row in scaling}) != 1:
        raise RuntimeError("Scaling cases use different math-library thread settings")
    if len({row["finite_difference"]["configured_workers"] for row in scaling}) != 1:
        raise RuntimeError("Scaling cases use different FD worker counts")
    sweep_path = next(root.rglob("step_sweep.json"))
    sweep = json.loads(sweep_path.read_text(encoding="utf-8"))["steps"]
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.8), layout="constrained")
    counts = [row["design_variables"] for row in scaling]
    totals = [row["including_common_baseline"] for row in scaling]
    ax = axes[0, 0]
    ax.loglog(counts, [row["complete_adjoint"]["wall_seconds"] for row in scaling], "o-", label="Adjoint elapsed")
    ax.loglog(counts, [row["fd_elapsed_wall_seconds"] for row in totals], "s-", label="Centered FD elapsed")
    ax.loglog(counts, [row["fd_aggregate_compute_wall_seconds"] for row in totals], "^-", label="Centered FD aggregate work")
    ax.set_ylabel("Runtime (s)")
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    ax.loglog(counts, [row["fd_forward_solve_count"] for row in totals], "s-", label="Centered FD")
    ax.loglog(counts, [row["adjoint_forward_solve_count"] for row in totals], "o-", label="Adjoint")
    ax.set_ylabel("Forward solves, including baseline")
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    ax.loglog(counts, [row["measured_elapsed_fd_over_adjoint_ratio"] for row in totals], "s-", label="Elapsed wall time")
    ax.loglog(counts, [row["aggregate_fd_over_adjoint_cpu_ratio"] for row in totals], "o-", label="Aggregate CPU time")
    ax.set_ylabel("Centered FD / adjoint time")
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    for label, display, marker in (("normalized_power", "Normalized power", "o-"), ("side_fraction", "Side flow", "s-")):
        ax.loglog([row["step"] for row in sweep],
                  [row["gradient_error"][label]["absolute_error_linf"] for row in sweep], marker, label=display)
    ax.set(xlabel="Centered step size", ylabel="Largest sampled derivative error")
    ax.legend(fontsize=8)
    for ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
        ax.set_xticks(counts, labels=[str(count) for count in counts])
        ax.set_xlabel("Editable density variables")
    for ax in axes.flat:
        ax.grid(alpha=.22, which="major")
    for suffix in (".pdf", ".png"):
        fig.savefig(output / f"fixed_pipe_sensitivity_cost{suffix}", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("case", choices=CASES)
    run.add_argument("output", type=Path)
    run.add_argument("--reference-metrics", type=Path, required=True)
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--step", type=float, default=.0001)
    run.add_argument("--coordinate-limit", type=int)
    run.add_argument("--skip-sweep", action="store_true")
    figures = sub.add_parser("plot")
    figures.add_argument("root", type=Path)
    figures.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "run":
        study(args)
    else:
        plot(args.root, args.output)


if __name__ == "__main__":
    main()
