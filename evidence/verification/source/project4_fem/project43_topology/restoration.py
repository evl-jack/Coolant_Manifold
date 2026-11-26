"""SLSQP Phase I epigraph restoration for the engineering inequalities."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize


FEASIBILITY_TOLERANCE = 1e-6


@dataclass(frozen=True)
class RestorationPoint:
    design: list[float]
    slack: float
    physical_objective: float
    power_w_per_m: float | None
    inequalities: list[float]
    engineering_violation: float
    margin_violation: float
    epigraph_violation: float
    engineering_feasible: bool
    margin_feasible: bool
    epigraph_feasible: bool


@dataclass(frozen=True)
class RestorationResult:
    margin: float
    feasibility_tolerance: float
    max_evaluations: int
    initial: RestorationPoint | None
    terminal_design: list[float] | None
    terminal_slack: float | None
    terminal_source: str
    engineering_violation: float | None
    margin_violation: float | None
    epigraph_violation: float | None
    accepted_design: list[float] | None
    accepted_source: str | None
    accepted: RestorationPoint | None
    status: int | None
    message: str
    optimizer_success: bool
    objective_calls: int
    physical_attempts: int
    unique_evaluations: int
    iterations: int | None
    trial_history: tuple[RestorationPoint, ...]
    iteration_history: tuple[RestorationPoint, ...]
    failed_design: list[float] | None


class _RestorationStopped(RuntimeError):
    pass


def _point(augmented, inequalities, margin, physical_objective, power, slack_weights):
    design, slack = augmented[:-1], float(augmented[-1])
    box = float(max(0.0, -design.min(), design.max() - 1.0))
    maximum = float(inequalities.max())
    engineering = max(box, maximum)
    tightened = inequalities + margin * slack_weights
    margin_violation = max(box, float(tightened.max()))
    epigraph = max(box, -slack, float((tightened - slack * slack_weights).max()))
    return RestorationPoint(
        design=design.tolist(), slack=slack, physical_objective=physical_objective,
        power_w_per_m=power, inequalities=inequalities.tolist(),
        engineering_violation=engineering, margin_violation=margin_violation,
        epigraph_violation=epigraph,
        engineering_feasible=engineering <= FEASIBILITY_TOLERANCE,
        margin_feasible=margin_violation <= FEASIBILITY_TOLERANCE,
        epigraph_feasible=epigraph <= FEASIBILITY_TOLERANCE,
    )


def restore(evaluator, initial, *, margin=0.05, max_evaluations=250, hard_tail=0, progress=None):
    """Minimize t with g(m) + margin - t <= 0, 0 <= m <= 1, and t >= 0.

    A feasible engineering design and a feasible epigraph point are separate
    outcomes. A design is eligible for handoff only when its engineering and
    tightened inequalities independently pass at t=0. The returned native
    endpoint and any accepted evaluated trial are identified separately.
    The last hard_tail constraints stay outside the epigraph slack and margin.
    """
    initial = np.array(initial, dtype=float, copy=True)
    if (
        initial.ndim != 1 or not initial.size or not np.all(np.isfinite(initial))
        or np.any((initial < 0) | (initial > 1))
    ):
        raise ValueError("Expected a finite initial design in the unit box")
    if not np.isfinite(margin) or margin < 0 or max_evaluations < 1:
        raise ValueError("Margin must be finite and nonnegative and budget must be positive")
    if not isinstance(hard_tail, int) or hard_tail < 0:
        raise ValueError("The number of hard constraints must be a nonnegative integer")

    evaluate = getattr(evaluator, "evaluate", evaluator)
    physical = {}
    trials = {}
    iteration_history = []
    objective_gradient = np.r_[np.zeros(initial.size), 1.0]
    physical_attempts = 0
    physical_requests = 0
    objective_calls = 0
    failed_design = None
    constraint_count = None
    slack_weights = None

    def sample(design):
        nonlocal physical_attempts, physical_requests, failed_design, constraint_count, slack_weights
        physical_requests += 1
        roundoff = 8 * np.finfo(float).eps
        if not np.all(np.isfinite(design)):
            raise _RestorationStopped("SLSQP supplied a nonfinite physical design")
        if np.any((design < -roundoff) | (design > 1.0 + roundoff)):
            failed_design = design.tolist()
            raise _RestorationStopped("SLSQP supplied a physical design outside its bounds")
        bounded = np.clip(design, 0.0, 1.0)
        key = bounded.tobytes()
        if key in physical:
            return physical[key]
        if physical_attempts >= max_evaluations:
            raise _RestorationStopped("Unique physical design evaluation limit reached")

        physical_attempts += 1
        try:
            value = evaluate(bounded)
            physical_objective = float(value.objective)
            power = getattr(value, "power_w_per_m", None)
            power = None if power is None else float(power)
            inequalities = np.array(value.constraints, dtype=float, copy=True)
            jacobian = np.array(value.jacobian, dtype=float, copy=True)
            if constraint_count is None:
                constraint_count = inequalities.size
                if hard_tail >= constraint_count:
                    raise ValueError("At least one constraint must be eligible for restoration slack")
                slack_weights = np.ones(constraint_count)
                slack_weights[constraint_count - hard_tail:] = 0.
            if (
                not constraint_count or inequalities.shape != (constraint_count,)
                or jacobian.shape != (constraint_count, design.size)
                or not np.all(np.isfinite(inequalities)) or not np.all(np.isfinite(jacobian))
                or not np.isfinite(physical_objective) or (power is not None and not np.isfinite(power))
            ):
                raise ValueError("Nonfinite or inconsistent restoration response")
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            failed_design = bounded.tolist()
            raise _RestorationStopped(f"{type(error).__name__}: {error}") from error

        physical[key] = inequalities, jacobian, physical_objective, power
        return physical[key]

    def observe(augmented):
        if not np.isfinite(augmented[-1]):
            raise _RestorationStopped("SLSQP supplied a nonfinite slack")
        inequalities, jacobian, physical_objective, power = sample(augmented[:-1])
        key = augmented.tobytes()
        if key not in trials:
            trials[key] = _point(augmented, inequalities, margin, physical_objective, power, slack_weights)
        if progress is not None:
            progress(trials[key], physical_attempts, physical_requests)
        return trials[key], inequalities, jacobian

    def objective(augmented):
        nonlocal objective_calls
        objective_calls += 1
        return float(augmented[-1]), objective_gradient

    def slsqp_inequalities(augmented):
        _, inequalities, _ = observe(augmented)
        return (augmented[-1] - margin) * slack_weights - inequalities

    def slsqp_jacobian(augmented):
        _, _, jacobian = observe(augmented)
        return np.column_stack((-jacobian, slack_weights))

    initial_record = terminal_record = None
    terminal_design = initial.tolist()
    terminal_slack = None
    terminal_source = "initial design"
    status = iterations = None
    success = False
    try:
        initial_constraints = sample(initial)[0]
        t0 = max(0.0, float(initial_constraints[slack_weights > 0].max()) + margin) + 1.0
        start = np.r_[initial, t0]
        initial_record, _, _ = observe(start)
        terminal_record = initial_record
        terminal_slack = t0
        outcome = minimize(
            objective, start, method="SLSQP", jac=True,
            bounds=[(0.0, 1.0)] * initial.size + [(0.0, None)],
            constraints=[{"type": "ineq", "fun": slsqp_inequalities, "jac": slsqp_jacobian}],
            callback=lambda augmented: iteration_history.append(observe(augmented)[0]),
            # Slack-only steps may advance SLSQP without consuming a physical evaluation.
            options={"ftol": 1e-10, "maxiter": max(100, 4 * max_evaluations)},
        )
        status, iterations = int(outcome.status), int(outcome.nit)
        success, message = bool(outcome.success), str(outcome.message)
        terminal_source = "optimizer result"
        terminal_record = None
        terminal_design = outcome.x[:-1].tolist() if np.all(np.isfinite(outcome.x[:-1])) else None
        terminal_slack = float(outcome.x[-1]) if np.isfinite(outcome.x[-1]) else None
        terminal_record, _, _ = observe(outcome.x)
    except _RestorationStopped as error:
        message = str(error) if status is None else f"{message} | endpoint assessment: {error}"
        if terminal_source != "optimizer result" and iteration_history:
            terminal_source = "last SLSQP major iteration"
            terminal_record = iteration_history[-1]
            terminal_design, terminal_slack = terminal_record.design, terminal_record.slack
            iterations = len(iteration_history)

    accepted = None
    accepted_source = None
    if terminal_record and terminal_record.engineering_feasible and terminal_record.margin_feasible:
        accepted = terminal_record
        accepted_source = terminal_source
    else:
        feasible = (point for point in trials.values() if point.engineering_feasible and point.margin_feasible)
        accepted = min(feasible, key=lambda point: max(point.inequalities), default=None)
        if accepted is not None:
            accepted_source = "margin-feasible evaluated trial"

    return RestorationResult(
        margin=float(margin), feasibility_tolerance=FEASIBILITY_TOLERANCE,
        max_evaluations=max_evaluations, initial=initial_record, terminal_design=terminal_design,
        terminal_slack=terminal_slack, terminal_source=terminal_source,
        engineering_violation=None if terminal_record is None else terminal_record.engineering_violation,
        margin_violation=None if terminal_record is None else terminal_record.margin_violation,
        epigraph_violation=None if terminal_record is None else terminal_record.epigraph_violation,
        accepted_design=None if accepted is None else accepted.design, accepted_source=accepted_source,
        accepted=accepted, status=status, message=message, optimizer_success=success,
        objective_calls=objective_calls, physical_attempts=physical_attempts, unique_evaluations=len(physical),
        iterations=iterations, trial_history=tuple(trials.values()), iteration_history=tuple(iteration_history),
        failed_design=failed_design,
    )
