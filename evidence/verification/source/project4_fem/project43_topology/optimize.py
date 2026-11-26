"""Adapters for bounded optimization with shared, nonpositive inequalities."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize


@dataclass(frozen=True)
class EvaluationRecord:
    design: np.ndarray
    objective: float
    inequalities: np.ndarray
    feasible: bool


@dataclass(frozen=True)
class OptimizationResult:
    algorithm: str
    terminal_design: np.ndarray
    terminal_source: str
    best_feasible_design: np.ndarray | None
    status: int | None
    message: str
    optimizer_success: bool
    objective_calls: int
    evaluator_requests: int
    unique_evaluations: int
    iterations: int | None
    iteration_history: tuple[np.ndarray, ...]
    evaluation_history: tuple[EvaluationRecord, ...]
    feasible_history: tuple[EvaluationRecord, ...]
    failed_design: np.ndarray | None


@dataclass(frozen=True)
class _Sample:
    objective: float
    gradient: np.ndarray
    constraints: np.ndarray
    jacobian: np.ndarray


class _EvaluationLimit(RuntimeError):
    pass


class _EvaluationFailure(RuntimeError):
    pass


_MMA_STATUS = {
    -5: "forced stop",
    -4: "roundoff limited progress",
    -3: "out of memory",
    -2: "invalid arguments",
    -1: "optimization failure",
    1: "success",
    2: "stop value reached",
    3: "objective tolerance reached",
    4: "design tolerance reached",
    5: "objective evaluation limit reached",
    6: "time limit reached",
}


def optimize(
    evaluator,
    initial,
    *,
    algorithm,
    max_evaluations,
    max_iterations=None,
    ftol=1e-10,
    xtol=1e-8,
    constraint_tolerance=1e-6,
    mma_verbosity=1,
    progress=None,
):
    """Run one stage without substituting a trial point for the optimizer result.

    The evaluator returns objective, gradient, constraints, and jacobian. The
    constraints have the common convention g <= 0, and the Jacobian has one
    row per constraint. The evaluator owns PDE work counters. This adapter
    tracks distinct successfully sampled designs and native objective calls.

    SciPy exposes major-iteration callbacks. NLopt's Python API exposes
    objective evaluations, and its native verbosity log reports internal
    progress. This adapter leaves MMA iteration counts unset.
    """
    algorithm = algorithm.upper()
    if algorithm not in {"SLSQP", "MMA"}:
        raise ValueError("algorithm must be SLSQP or MMA")
    if max_evaluations < 1 or (max_iterations is not None and max_iterations < 1):
        raise ValueError("Evaluation and iteration limits must be positive")
    if algorithm == "MMA" and max_iterations is not None:
        raise ValueError("NLopt MMA supports an evaluation limit, not an iteration limit")
    if not np.all(np.isfinite([ftol, xtol, constraint_tolerance])) or min(
        ftol, xtol, constraint_tolerance
    ) < 0:
        raise ValueError("Optimizer tolerances must be finite and nonnegative")

    initial = np.array(initial, dtype=float, copy=True)
    if (
        initial.ndim != 1
        or not initial.size
        or not np.all(np.isfinite(initial))
        or np.any((initial < 0) | (initial > 1))
    ):
        raise ValueError("Expected a finite initial design in the unit box")

    evaluate = getattr(evaluator, "evaluate", evaluator)
    history = []
    iteration_history = []
    cache = {}
    constraint_count = None
    objective_calls = 0
    evaluator_requests = 0
    failed_design = None

    def sample(design):
        nonlocal constraint_count, evaluator_requests, failed_design
        evaluator_requests += 1
        # SLSQP can exceed the bounds by floating-point roundoff at a callback.
        tolerance = 8 * np.finfo(float).eps
        if not np.all(np.isfinite(design)) or np.any(
            (design < -tolerance) | (design > 1 + tolerance)
        ):
            failed_design = design.copy()
            raise _EvaluationFailure("Optimizer supplied a nonfinite or out-of-box design")
        design = np.clip(design, 0.0, 1.0)
        key = design.tobytes()
        if key in cache:
            value, record = cache[key]
            if progress is not None:
                progress(record, len(cache), evaluator_requests)
            return value
        if len(cache) >= max_evaluations:
            raise _EvaluationLimit("Unique design evaluation limit reached")

        try:
            value = evaluate(design)
            objective = float(value.objective)
            gradient = np.array(value.gradient, dtype=float, copy=True)
            constraints = np.array(value.constraints, dtype=float, copy=True)
            jacobian = np.array(value.jacobian, dtype=float, copy=True)
            if constraint_count is None:
                constraint_count = constraints.size
            if (
                gradient.shape != design.shape
                or constraints.shape != (constraint_count,)
                or jacobian.shape != (constraint_count, design.size)
                or not np.isfinite(objective)
                or not np.all(np.isfinite(gradient))
                or not np.all(np.isfinite(constraints))
                or not np.all(np.isfinite(jacobian))
            ):
                raise ValueError("Nonfinite or inconsistent evaluator response")
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            failed_design = design.copy()
            raise _EvaluationFailure(f"{type(error).__name__}: {error}") from error

        feasible = bool(np.all(constraints <= constraint_tolerance))
        record = EvaluationRecord(design.copy(), objective, constraints, feasible)
        history.append(record)
        value = _Sample(objective, gradient, constraints, jacobian)
        cache[key] = value, record
        if progress is not None:
            progress(record, len(cache), evaluator_requests)
        return value

    def slsqp_objective(design):
        nonlocal objective_calls
        objective_calls += 1
        value = sample(design)
        return value.objective, value.gradient

    def mma_objective(design, gradient):
        nonlocal objective_calls
        objective_calls += 1
        value = sample(design)
        if gradient.size:
            gradient[:] = value.gradient
        return value.objective

    def mma_constraints(inequalities, design, jacobian):
        value = sample(design)
        inequalities[:] = value.constraints
        if jacobian.size:
            jacobian[:] = value.jacobian

    terminal = initial
    terminal_source = "initial design"
    status = None
    iterations = None
    success = False
    try:
        sample(initial)
        if algorithm == "SLSQP":
            constraints = []
            if constraint_count:
                constraints.append(
                    {
                        "type": "ineq",
                        "fun": lambda design: -sample(design).constraints,
                        "jac": lambda design: -sample(design).jacobian,
                    }
                )
            result = minimize(
                slsqp_objective,
                initial,
                jac=True,
                method="SLSQP",
                bounds=[(0.0, 1.0)] * initial.size,
                constraints=constraints,
                callback=lambda design: iteration_history.append(design.copy()),
                options={"maxiter": max_iterations or max_evaluations, "ftol": ftol},
            )
            terminal = result.x.copy()
            status = int(result.status)
            message = str(result.message)
            iterations = int(result.nit)
            success = bool(result.success)
        else:
            try:
                import nlopt
            except ImportError as error:
                raise ImportError("NLopt with LD_MMA is required for the MMA stage") from error

            optimizer = nlopt.opt(nlopt.LD_MMA, initial.size)
            optimizer.set_lower_bounds(0.0)
            optimizer.set_upper_bounds(1.0)
            optimizer.set_min_objective(mma_objective)
            if constraint_count:
                optimizer.add_inequality_mconstraint(
                    mma_constraints, np.full(constraint_count, constraint_tolerance)
                )
            optimizer.set_maxeval(max_evaluations)
            optimizer.set_ftol_abs(ftol)
            optimizer.set_xtol_abs(xtol)
            optimizer.set_param("verbosity", mma_verbosity)
            # Some supported Python bindings always raise on native early stops.
            configure_exceptions = getattr(optimizer, "set_exceptions_enabled", None)
            if configure_exceptions is not None:
                configure_exceptions(False)
            try:
                terminal = optimizer.optimize(initial).copy()
            except (nlopt.RoundoffLimited, nlopt.ForcedStop) as error:
                status = int(optimizer.last_optimize_result())
                message = _MMA_STATUS.get(status, str(error))
                if history:
                    terminal = history[-1].design.copy()
                    terminal_source = "last evaluated MMA trial after native exception"
            else:
                status = int(optimizer.last_optimize_result())
                message = _MMA_STATUS.get(status, f"NLopt result {status}")
                success = status in {1, 2, 3, 4}
                terminal_source = "optimizer result"
        if algorithm == "SLSQP":
            terminal_source = "optimizer result"
    except (_EvaluationLimit, _EvaluationFailure) as error:
        message = str(error)
        if iteration_history:
            terminal = iteration_history[-1].copy()
            terminal_source = "last SLSQP major iteration"
            iterations = len(iteration_history)

    feasible_history = tuple(record for record in history if record.feasible)
    best = min(feasible_history, key=lambda record: record.objective, default=None)
    return OptimizationResult(
        algorithm=algorithm,
        terminal_design=terminal,
        terminal_source=terminal_source,
        best_feasible_design=None if best is None else best.design.copy(),
        status=status,
        message=message,
        optimizer_success=success,
        objective_calls=objective_calls,
        evaluator_requests=evaluator_requests,
        unique_evaluations=len(history),
        iterations=iterations,
        iteration_history=tuple(iteration_history),
        evaluation_history=tuple(history),
        feasible_history=feasible_history,
        failed_design=failed_design,
    )
