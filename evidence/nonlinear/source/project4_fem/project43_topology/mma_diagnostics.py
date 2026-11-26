"""Replay the known MMA problem with native scalar, native vector and shared adapters."""

import argparse
import json
import os
from pathlib import Path
from types import SimpleNamespace

import nlopt
import numpy as np

from .kkt_audit import audit_kkt
from .optimize import optimize


FEASIBILITY_TOLERANCE = 1e-6
DESIGN_TOLERANCE = 1e-4
FD_TOLERANCE = 1e-8
EXPECTED_DESIGN = np.array([0.5, 0.5])
CONSTRAINTS = {
    "gate2": ("minimum_sum", "maximum_sum", "inactive_constant"),
    "two_inequalities": ("minimum_sum", "maximum_sum"),
    "maximum_sum_only": ("maximum_sum",),
}


def independent_values(design, names):
    x, y = map(float, design)
    values = {"minimum_sum": 0.6 - x - y, "maximum_sum": x + y - 1.0,
              "inactive_constant": -1.0}
    return np.array([(x - 0.8)**2 + (y - 0.8)**2, *(values[name] for name in names)])


class KnownProblem:
    def __init__(self, names):
        self.names = names
        self.history = []

    def evaluate(self, design):
        delta = design - 0.8
        total = design.sum()
        constraints = {
            "minimum_sum": (0.6 - total, (-1.0, -1.0)),
            "maximum_sum": (total - 1.0, (1.0, 1.0)),
            "inactive_constant": (-1.0, (0.0, 0.0)),
        }
        values, jacobian = zip(*(constraints[name] for name in self.names))
        return SimpleNamespace(objective=float(delta @ delta), gradient=2 * delta,
                               constraints=np.array(values), jacobian=np.array(jacobian))

    def objective(self, design, gradient):
        value = self.evaluate(design)
        if gradient.size:
            gradient[:] = value.gradient
        self.history.append(design.copy())
        return value.objective

    def scalar_constraint(self, index):
        def constraint(design, gradient):
            value = self.evaluate(design)
            if gradient.size:
                gradient[:] = value.jacobian[index]
            return value.constraints[index]

        return constraint

    def vector_constraint(self, result, design, gradient):
        value = self.evaluate(design)
        result[:] = value.constraints
        if gradient.size:
            gradient[:] = value.jacobian


def callback_parity(names):
    problem = KnownProblem(names)
    records = []
    for design in [np.zeros(2), np.array([0.23, 0.37]), EXPECTED_DESIGN, np.array([0.9, 0.1])]:
        references = independent_values(design, names)
        step = 1e-6
        differences = []
        for index in range(design.size):
            direction = np.eye(design.size)[index] * step
            # The origin is on the box boundary. The three-point forward
            # stencil stays in the box and is exact for this quadratic.
            if design[index] == 0:
                once = independent_values(design + direction, names)
                twice = independent_values(design + 2 * direction, names)
                differences.append((-3 * references + 4 * once - twice) / (2 * step))
            else:
                above = independent_values(design + direction, names)
                below = independent_values(design - direction, names)
                differences.append((above - below) / (2 * step))
        fd = np.column_stack(differences)

        direct = problem.evaluate(design)
        objective_gradient = np.full(design.size, np.nan)
        objective = problem.objective(design, objective_gradient)
        scalar_gradients = np.full_like(direct.jacobian, np.nan)
        scalars = [problem.scalar_constraint(row)(design, scalar_gradients[row])
                   for row in range(len(names))]
        vector_gradients = np.full_like(direct.jacobian, np.nan)
        vector = np.full(len(names), np.nan)
        problem.vector_constraint(vector, design, vector_gradients)
        values = {"project_evaluator": np.r_[direct.objective, direct.constraints],
                  "native_scalar": np.r_[objective, scalars],
                  "native_vector": np.r_[objective, vector]}
        jacobians = {"project_evaluator": np.vstack((direct.gradient, direct.jacobian)),
                     "native_scalar": np.vstack((objective_gradient, scalar_gradients)),
                     "native_vector": np.vstack((objective_gradient, vector_gradients))}
        records.append({
            "design": design.tolist(), "references": references.tolist(), "reference_fd": fd.tolist(),
            "max_value_error": {key: max_error(value, references)
                                for key, value in values.items()},
            "max_derivative_error": {key: max_error(value, fd)
                                     for key, value in jacobians.items()},
        })

    errors = [error for row in records for field in ("max_value_error", "max_derivative_error")
              for error in row[field].values()]
    maximum = None if any(error is None for error in errors) else max(errors)
    return {"passed": maximum is not None and maximum <= FD_TOLERANCE, "threshold": FD_TOLERANCE,
            "max_error": maximum, "points": records}


def max_error(actual, reference):
    return float(np.max(np.abs(actual - reference))) if np.all(np.isfinite(actual)) else None


def describe_design(problem, design):
    if design is None:
        return None
    values = independent_values(design, problem.names)
    x, y = map(float, design)
    independent_gradient = np.array([2 * x - 1.6, 2 * y - 1.6])
    rows = {"minimum_sum": (-1.0, -1.0), "maximum_sum": (1.0, 1.0),
            "inactive_constant": (0.0, 0.0)}
    independent_jacobian = np.array([rows[name] for name in problem.names])
    report = audit_kkt(design, independent_gradient, values[1:], independent_jacobian)
    design_error = float(np.max(np.abs(design - EXPECTED_DESIGN)))
    return {"design": design.tolist(), "objective": float(values[0]),
            "inequalities": values[1:].tolist(), "max_constraint": float(np.max(values[1:])),
            "design_error": design_error, "feasible": report.primal_violation <= FEASIBILITY_TOLERANCE,
            "kkt_verified": report.verified, "kkt_multipliers": report.multipliers.tolist(),
            "kkt_stationarity": report.stationarity_inf,
            "passed": design_error <= DESIGN_TOLERANCE and report.verified}


def summarize_trials(problem, designs):
    trials = []
    best = None
    for design in designs:
        values = independent_values(design, problem.names)
        feasible = bool(np.max(values[1:]) <= FEASIBILITY_TOLERANCE and
                        np.all((design >= 0) & (design <= 1)))
        trials.append({"design": design.tolist(), "objective": float(values[0]),
                       "inequalities": values[1:].tolist(), "feasible": feasible})
        if feasible and (best is None or values[0] < best[0]):
            best = values[0], design
    return {"trials": trials, "feasible_trials": sum(row["feasible"] for row in trials),
            "best_feasible": describe_design(problem, None if best is None else best[1])}


def run_native(problem, initial, mode, *, max_evaluations, verbosity):
    config = nlopt.opt(nlopt.LD_MMA, initial.size)
    config.set_lower_bounds(0.0)
    config.set_upper_bounds(1.0)
    config.set_min_objective(problem.objective)
    if mode == "native_scalar":
        for row in range(len(problem.names)):
            config.add_inequality_constraint(problem.scalar_constraint(row), FEASIBILITY_TOLERANCE)
    else:
        config.add_inequality_mconstraint(
            problem.vector_constraint, np.full(len(problem.names), FEASIBILITY_TOLERANCE)
        )
    config.set_maxeval(max_evaluations)
    config.set_ftol_abs(1e-10)
    config.set_xtol_abs(1e-8)
    config.set_param("verbosity", verbosity)
    configure_exceptions = getattr(config, "set_exceptions_enabled", None)
    if configure_exceptions is not None:
        configure_exceptions(False)
    design = None
    error = None
    try:
        design = config.optimize(initial.copy())
    except (RuntimeError, ValueError, nlopt.RoundoffLimited, nlopt.ForcedStop) as exception:
        error = f"{type(exception).__name__}: {exception}"
    objective = float(config.last_optimum_value())
    return {"status": int(config.last_optimize_result()), "exception": error,
            "objective_calls": int(config.get_numevals()),
            "reported_objective": objective if np.isfinite(objective) else None,
            "terminal": describe_design(problem, design),
            **summarize_trials(problem, problem.history)}


def run_project(problem, initial, algorithm, *, max_evaluations, verbosity):
    outcome = optimize(problem, initial, algorithm=algorithm, max_evaluations=max_evaluations,
                       mma_verbosity=verbosity)
    return {"status": outcome.status, "message": outcome.message,
            "optimizer_success": outcome.optimizer_success, "terminal_source": outcome.terminal_source,
            "objective_calls": outcome.objective_calls, "evaluator_requests": outcome.evaluator_requests,
            "iterations": outcome.iterations, "unique_evaluations": outcome.unique_evaluations,
            "terminal": describe_design(problem, outcome.terminal_design),
            **summarize_trials(problem, [row.design for row in outcome.evaluation_history])}


def run(max_evaluations, verbosity):
    cases = []
    for fixture, names in CONSTRAINTS.items():
        for start_name, start in {
            "zero": np.zeros(2), "feasible": np.array([0.4, 0.4]),
            "upper_infeasible": np.array([0.8, 0.8]),
        }.items():
            initial = describe_design(KnownProblem(names), start)
            outcomes = {}
            for mode in ("native_scalar", "native_vector", "project_mma", "project_slsqp"):
                print(f"{fixture} {start_name} {mode}", flush=True)
                problem = KnownProblem(names)
                options = {"max_evaluations": max_evaluations, "verbosity": verbosity}
                try:
                    if mode.startswith("native"):
                        outcomes[mode] = run_native(problem, start, mode, **options)
                    else:
                        outcomes[mode] = run_project(problem, start, mode.removeprefix("project_"), **options)
                except (RuntimeError, ValueError) as error:
                    outcomes[mode] = {"exception": f"{type(error).__name__}: {error}"}
            cases.append({"fixture": fixture, "start": start_name, "initial": initial,
                          "constraint_names": names, "outcomes": outcomes})
    return {"nlopt_version": f"{nlopt.version_major()}.{nlopt.version_minor()}.{nlopt.version_bugfix()}",
            "job_id": os.environ["SLURM_JOB_ID"], "max_evaluations": max_evaluations,
            "expected_design": EXPECTED_DESIGN.tolist(), "expected_objective": 0.18,
            "design_tolerance": DESIGN_TOLERANCE, "feasibility_tolerance": FEASIBILITY_TOLERANCE,
            "callback_parity": {fixture: callback_parity(names) for fixture, names in CONSTRAINTS.items()},
            "cases": cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New JSON report file")
    parser.add_argument("--max-evaluations", type=int, default=300)
    parser.add_argument("--verbosity", type=int, default=0, help="Native MMA verbosity")
    arguments = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical verification as a Sol Slurm job")
    if arguments.max_evaluations < 1:
        parser.error("--max-evaluations must be positive")
    with arguments.output.open("x", encoding="utf-8") as output:
        json.dump(run(arguments.max_evaluations, arguments.verbosity), output, indent=2, allow_nan=False)
        output.write("\n")


if __name__ == "__main__":
    main()
