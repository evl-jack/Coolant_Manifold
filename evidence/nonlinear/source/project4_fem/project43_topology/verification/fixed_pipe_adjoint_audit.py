"""Replay the fixed-pipe pilot and verify the full reduced derivatives."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ...fem import State
from ..evaluator import Stage
from ..fixed_pipe_pilot_model import SIDE_CONSTRAINT_SCALE, FixedPipeEvaluator, build_fixed_pipe
from ..verify_flow import (
    ADJOINT_LIMIT, EQUIVALENCE_TOLERANCE, MASS_LIMIT, RESIDUAL_DERIVATIVE_LIMIT,
    RESIDUAL_LIMIT, TAYLOR_ORDER_LIMIT, TRANSFER_TOLERANCE, _power_and_flux,
)


STEP_SIZES = (1e-2, 5e-3, 2e-3, 1e-3, 1e-4, 1e-5, 1e-6)
RESPONSES = ("normalized_power", "side_fraction")
CASES = ("original", "gray", "beta_0", "beta_2", "beta_8")
SEED = 432618


def relative_error(actual, expected, floor=1e-10):
    difference = float(np.max(np.abs(actual - expected)))
    return difference / max(float(np.max(np.abs(actual))), float(np.max(np.abs(expected))), floor)


def independent_filter(pipeline):
    """Build the 5-point cell-neighbor operator and its eliminated anchors."""
    grid = pipeline.grid
    centers = grid.centers_mm
    rows = np.flatnonzero(pipeline.free)
    row_at_cell = {cell: row for row, cell in enumerate(rows)}
    step = grid.spacing_mm
    radius_squared = pipeline.filter_radius_mm ** 2
    mass = grid.cell_area_mm2
    matrix = mass * np.eye(len(rows))
    fixed_rhs = np.zeros(len(rows))
    ny, nx = grid.shape
    for row, cell in enumerate(rows):
        iy, ix = np.unravel_index(cell, grid.shape)
        for y, x in ((iy, ix - 1), (iy, ix + 1), (iy - 1, ix), (iy + 1, ix)):
            matrix[row, row] += radius_squared
            if y == ny:
                near_collar = 9.5 < centers[cell, 0] < 10.5 and 6 < centers[cell, 1] + step < 7
                fixed_rhs[row] += radius_squared * near_collar
                continue
            if not (0 <= y < ny and 0 <= x < nx):
                raise ValueError("Unexpected editable cell beside the design-grid perimeter")
            neighbor = y * nx + x
            if neighbor in row_at_cell:
                matrix[row, row_at_cell[neighbor]] = -radius_squared
            else:
                fixed_rhs[row] += radius_squared * grid.base_raw[neighbor]
    return matrix, fixed_rhs


def density_action(pipeline, matrix, density, stage, direction):
    filtered_action = np.linalg.solve(matrix, pipeline.grid.cell_area_mm2 * direction)
    filtered = density.filtered[pipeline.free]
    if stage.beta:
        slope = stage.beta / (2 * np.tanh(stage.beta / 2) * np.cosh(stage.beta * (filtered - .5)) ** 2)
    else:
        slope = np.ones_like(filtered)
    projected = density.projected[pipeline.free]
    resistance_slope = -stage.alpha_max * stage.q * (stage.q + 1) / (stage.q + projected) ** 2
    resistance_action = np.zeros(pipeline.grid.size)
    resistance_action[pipeline.free] = resistance_slope * slope * filtered_action
    return filtered_action, resistance_action


def filter_transfer_checks(evaluator, matrix, fixed_rhs, design, rng):
    pipeline, domain = evaluator.pipeline, evaluator.domain
    stage = evaluator.stage
    density = evaluator.density(design)
    editable = pipeline.free
    anchors = pipeline.grid.base_raw[~editable]
    mass = pipeline.grid.cell_area_mm2
    direction = rng.uniform(-1, 1, pipeline.size)
    cotangent = rng.standard_normal(pipeline.grid.size)
    exact_state = np.linalg.solve(matrix, mass * design + fixed_rhs)
    action, resistance_action = density_action(pipeline, matrix, density, stage, direction)
    step = 1e-5
    plus, minus = (evaluator.density(design + sign * step * direction) for sign in (1, -1))
    numerical_filter_action = (plus.filtered[editable] - minus.filtered[editable]) / (2 * step)
    numerical_resistance_action = (plus.resistance - minus.resistance) / (2 * step)
    production_reverse = pipeline.filter_pullback(cotangent)
    exact_reverse = mass * np.linalg.solve(matrix.T, cotangent[editable])

    cell_index = domain.cell_index
    active = cell_index >= 0
    scale = .004 ** 2 / .001002
    element_resistance = np.zeros(len(cell_index))
    element_resistance[active] = scale * density.resistance[cell_index[active]]
    products = rng.standard_normal(domain.fem.weights.shape)
    reverse_reference = np.zeros(pipeline.grid.size)
    np.add.at(reverse_reference, np.repeat(cell_index[active], products.shape[1]),
              scale * products[active].ravel())
    forward = evaluator.transfer.forward(density.resistance)
    reverse = evaluator.transfer.transpose(products)
    comparisons = {
        "filter_forward": relative_error(density.filtered[editable], exact_state),
        "filter_directional": relative_error(numerical_filter_action, action),
        "filter_transpose": relative_error(production_reverse, exact_reverse),
        "filter_dot_product": relative_error(direction @ production_reverse, action @ cotangent[editable], 1),
        "resistance_directional": relative_error(numerical_resistance_action, resistance_action),
        "transfer_forward": relative_error(forward[:, 0], element_resistance),
        "transfer_transpose": relative_error(reverse, reverse_reference),
        "transfer_dot_product": relative_error(np.sum(forward * products), density.resistance @ reverse, 1),
    }
    anchor_preservation = all(np.array_equal(getattr(density, key)[~editable], anchors)
                              for key in ("raw", "filtered", "projected"))
    fixed_fluid_drag_zero = bool(np.all(density.resistance[~editable][anchors == 1] == 0)
                                 and np.all(forward[~active] == 0))
    bounds = {key: [float(getattr(density, key).min()), float(getattr(density, key).max())]
              for key in ("raw", "filtered", "projected")}
    bounds_pass = all(-1e-14 <= lo <= hi <= 1 + 1e-14 for lo, hi in bounds.values())
    return {
        "experiment_beta": stage.beta, "condition_2": float(np.linalg.cond(matrix)),
        "errors": comparisons, "fixed_anchor_preservation": anchor_preservation,
        "fixed_fluid_zero_resistance": fixed_fluid_drag_zero, "density_bounds": bounds,
        "finite_difference_step": step,
        "passed": anchor_preservation and fixed_fluid_drag_zero and bounds_pass and all(
            value <= (RESIDUAL_DERIVATIVE_LIMIT if name.endswith("directional") else TRANSFER_TOLERANCE)
            for name, value in comparisons.items()),
    }


def independent_state_actions(fem, state, direction, resistance):
    """Apply each weak-form tangent term by direct quadrature and scatter."""
    velocity, gradient = fem.velocity_fields(state)
    delta_velocity, delta_gradient = fem.velocity_fields(direction)
    delta_pressure = direction[2 * fem.n_u:][fem.triangles] @ fem.p_basis.T
    symmetric_gradient = delta_gradient + delta_gradient.swapaxes(-1, -2)
    weighted_gradient = fem.weights[..., None, None] * fem.gradient
    viscous = fem.viscosity * np.einsum("eqaj,eqij->eai", weighted_gradient, symmetric_gradient)
    pressure = -np.einsum("eqai,eq->eai", weighted_gradient, delta_pressure)
    divergence = np.trace(delta_gradient, axis1=-2, axis2=-1)
    continuity = np.zeros(fem.size)
    np.add.at(continuity, 2 * fem.n_u + fem.triangles,
              -(fem.weights * divergence) @ fem.p_basis)
    transport = np.einsum("eqij,eqj->eqi", gradient, delta_velocity)
    convection = np.einsum("eqij,eqj->eqi", delta_gradient, velocity)

    def velocity_action(force):
        return fem.scatter_velocity(np.einsum("eq,qa,eqi->eai", fem.weights, fem.basis, force))

    return {
        "viscous": fem.scatter_velocity(viscous), "pressure": fem.scatter_velocity(pressure),
        "continuity": continuity, "brinkman": velocity_action(resistance[..., None] * delta_velocity),
        "convection_transported": velocity_action(fem.inertia * transport),
        "convection_transporting": velocity_action(fem.inertia * convection),
    }


def state_jacobian_checks(evaluator, state, resistance, rng):
    domain, fem = evaluator.domain, evaluator.domain.fem
    free = state.free
    operator = fem.viscous_pressure + fem.drag_matrix(resistance)
    kinds = {"velocity": free[free < 2 * fem.n_u],
             "pressure": free[free >= 2 * fem.n_u], "combined": free}
    out = {}
    for kind, indices in kinds.items():
        direction = np.zeros(fem.size)
        direction[indices] = rng.uniform(-1, 1, len(indices))
        pieces = independent_state_actions(fem, state.values, direction, resistance)
        exact = sum(pieces.values())[free]
        tangent = state.tangent[free][:, free] @ direction[free]
        records = []
        for step in (1e-2, 1e-3, 1e-4):
            plus, minus = (fem.residual(state.values + sign * step * direction, operator, evaluator.load,
                                        tangent=False)[0][free] for sign in (1, -1))
            fd = (plus - minus) / (2 * step)
            records.append({"step": step, "fd_relative_error": relative_error(tangent, fd)})
        velocity = free < 2 * fem.n_u
        drag = fem.drag_matrix(resistance) @ direction
        block_errors = {
            "momentum": relative_error(tangent[velocity], exact[velocity]),
            "continuity": relative_error(tangent[~velocity], exact[~velocity]),
            "brinkman": relative_error(drag[free], pieces["brinkman"][free]),
            "complete": relative_error(tangent, exact),
        }
        out[kind] = {"errors": block_errors, "records": records,
                     "term_inf_norms": {key: float(np.linalg.norm(value[free], np.inf))
                                        for key, value in pieces.items()},
                     "passed": max(block_errors.values()) <= RESIDUAL_DERIVATIVE_LIMIT and all(
                         row["fd_relative_error"] <= RESIDUAL_DERIVATIVE_LIMIT for row in records)}
    return {"free_dofs": len(free), "eliminated_dofs": len(domain.fixed),
            "prescribed_values_preserved": bool(np.array_equal(state.values[domain.fixed], domain.prescribed)),
            "directions": out, "passed": all(row["passed"] for row in out.values())
            and bool(np.array_equal(state.values[domain.fixed], domain.prescribed))}


def feasible_direction(design, direction):
    """Use centered steps where possible, inward steps at (near-)active bounds."""
    direction = direction / np.linalg.norm(direction, np.inf)
    use = direction != 0
    centered_limit = np.min(np.minimum(design[use], 1 - design[use]) / np.abs(direction[use]))
    if centered_limit >= STEP_SIZES[0]:
        return direction, False
    direction = np.where(design < .5, np.abs(direction), -np.abs(direction))
    return direction, True


def directional_check(evaluator, matrix, design, value, state, direction, steps):
    direction, inward = feasible_direction(design, direction)
    fem = evaluator.domain.fem
    free = state.free
    density = evaluator.density(design)
    _, resistance_direction = density_action(evaluator.pipeline, matrix, density, evaluator.stage, direction)
    action = fem.drag_matrix(evaluator.transfer.forward(resistance_direction)) @ state.values
    adjoint = np.array([value.gradient @ direction, -SIDE_CONSTRAINT_SCALE * value.jacobian[0] @ direction])
    worst_residual = max(state.residual_history)
    final_residuals, balances = [], []

    def primal(candidate, solve):
        if np.min(candidate) < 0 or np.max(candidate) > 1:
            raise ValueError("Finite-difference perturbation leaves the raw-density bounds")
        resistance = evaluator.transfer.forward(evaluator.density(candidate).resistance)
        residual = fem.drag_matrix(resistance) @ state.values
        current = evaluator.solve(resistance, state.values) if solve else state
        power, flux = _power_and_flux(evaluator.domain, current.values)
        final_residuals.append(current.residual_history[-1])
        balances.append(abs(flux[1:].sum() - flux[0]) / flux[0])
        return np.array([power / evaluator.power_reference, flux[2] / flux[0]]), residual[free]

    base, base_residual = primal(design, False)
    independent_constraint = (evaluator.side_reference - base[1]) / SIDE_CONSTRAINT_SCALE
    response_parity = np.abs(base - [value.objective, value.side_fraction])
    constraint_parity = float(abs(independent_constraint - value.constraints[0]))
    records = []
    for step in steps:
        plus, residual_plus = primal(design + step * direction, True)
        if inward:
            farther, residual_farther = primal(design + 2 * step * direction, True)
            fd = (4 * (plus - base) - (farther - base)) / (2 * step)
            residual_fd = (4 * (residual_plus - base_residual) - (residual_farther - base_residual)) / (2 * step)
        else:
            minus, residual_minus = primal(design - step * direction, True)
            fd = (plus - minus) / (2 * step)
            residual_fd = (residual_plus - residual_minus) / (2 * step)
        absolute = np.abs(fd - adjoint)
        scale = np.maximum(np.maximum(np.abs(fd), np.abs(adjoint)), 1e-4)
        records.append({"step": step, "fd": fd.tolist(), "absolute_error": absolute.tolist(),
                        "relative_error": [float(error / max(abs(a), abs(b))) if max(abs(a), abs(b)) >= 1e-8 else None
                                           for error, a, b in zip(absolute, fd, adjoint)],
                        "screen_scaled_error": (absolute / scale).tolist(),
                        "fixed_state_resistance_tangent_error": relative_error(action[free], residual_fd),
                        "e0": np.abs(plus - base).tolist(),
                        "e1": np.abs(plus - base - step * adjoint).tolist()})

    # Require two consecutive slopes before the nominal response-noise floor.
    noise_floor = 5e-11
    orders = []
    for coarse, fine in zip(records, records[1:]):
        orders.append([float(np.log(old / new) / np.log(coarse["step"] / fine["step"]))
                       if min(old, new) > noise_floor else None for old, new in zip(coarse["e1"], fine["e1"])])
    windows = {}
    for i, response in enumerate(RESPONSES):
        windows[response] = next(({
            "steps": [records[j]["step"], records[j + 2]["step"]], "orders": [left[i], right[i]],
            "three_point_order": float(np.log(records[j]["e1"][i] / records[j + 2]["e1"][i])
                                       / np.log(records[j]["step"] / records[j + 2]["step"])),
        } for j, (left, right) in enumerate(zip(orders, orders[1:]))
            if all(v is not None and TAYLOR_ORDER_LIMIT <= v <= 2.2 for v in (left[i], right[i]))), None)

    derivative_passed = sum(max(row["screen_scaled_error"]) <= ADJOINT_LIMIT
                            and row["fixed_state_resistance_tangent_error"] <= RESIDUAL_DERIVATIVE_LIMIT
                            for row in records) >= 2
    return {"difference": "second-order inward" if inward else "centered", "direction": direction.tolist(),
            "direction_normalization": "infinity norm one", "adjoint": adjoint.tolist(),
            "base_response_absolute_parity": response_parity.tolist(), "base_constraint_abs_parity": constraint_parity,
            "baseline_peak_residual": worst_residual, "maximum_final_residual": max(final_residuals),
            "maximum_mass_imbalance": float(max(balances)), "forward_residual_tolerance": 1e-11,
            "nominal_taylor_noise_floor": noise_floor,
            "minimum_observed_e1": np.min([record["e1"] for record in records], axis=0).tolist(),
            "records": records, "orders": orders,
            "taylor_windows": windows, "derivative_passed": bool(derivative_passed),
            "taylor_passed": all(windows.values()), "passed": bool(derivative_passed and all(windows.values())
            and max(final_residuals) <= RESIDUAL_LIMIT and max(balances) <= MASS_LIMIT
            and max(response_parity.max(), constraint_parity) <= EQUIVALENCE_TOLERANCE)}


def snapshot(root, stage_directory, kind):
    folder = root / stage_directory / kind
    path = folder / "fields.npz"
    with np.load(path) as data:
        design, state = data["design"].copy(), data["state"].copy()
    metrics = json.loads((folder / "metrics.json").read_text(encoding="utf-8"))["results"]
    return design, state, metrics, hashlib.sha256(path.read_bytes()).hexdigest()


def checkpoint(path, payload):
    temporary = path.with_suffix(".partial.json")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def audit(root, case_names, direction_names, report_path):
    rng = np.random.default_rng(SEED)
    domain, pipeline = build_fixed_pipe()
    matrix, fixed_rhs = independent_filter(pipeline)
    gray = .12 + .76 * pipeline.initial + rng.uniform(-.025, .025, pipeline.size)
    folders = {0: "stage_01_beta_0", 2: "stage_02_beta_2", 8: "stage_03_beta_8"}
    centers = pipeline.grid.centers_mm[pipeline.free]
    near = np.argmin(np.sum((centers - [9.625, 4.125]) ** 2, axis=1))
    collar = np.argmin(np.sum((centers - [10.375, 5.875]) ** 2, axis=1))
    directions = {"random": rng.uniform(-1, 1, pipeline.size),
                  "junction_cell": np.eye(pipeline.size)[near], "collar_cell": np.eye(pipeline.size)[collar]}
    report = {"snapshot_root": str(root), "seed": SEED, "steps": STEP_SIZES, "responses": RESPONSES,
              "fem": {"degrees_of_freedom": domain.fem.size, "triangles": len(domain.fem.triangles)},
              "thresholds": {"transfer": TRANSFER_TOLERANCE, "port_parity": EQUIVALENCE_TOLERANCE,
                             "nonlinear_residual": RESIDUAL_LIMIT, "relative_mass": MASS_LIMIT,
                             "residual_derivative": RESIDUAL_DERIVATIVE_LIMIT, "scaled_adjoint": ADJOINT_LIMIT,
                             "minimum_taylor_order": TAYLOR_ORDER_LIMIT},
              "cases": {}, "coordinate_directions_mm": {"junction_cell": centers[near].tolist(),
                                                           "collar_cell": centers[collar].tolist()}}
    for name in case_names:
        beta = int(name.removeprefix("beta_")) if name.startswith("beta_") else 8
        stage = Stage(.10, beta, 1e5)
        original, original_state, metrics, original_hash = snapshot(root, folders[beta], "original")
        if name == "gray":
            design, initial, source_hash = gray, original_state, None
        elif name == "original":
            design, initial, source_hash = original, original_state, original_hash
        else:
            design, initial, _, source_hash = snapshot(root, folders[beta], "selected")
        evaluator = FixedPipeEvaluator(domain, pipeline, stage, metrics["power_w_per_m"], metrics["side_fraction"])
        # Evaluate the archived vector first. A copied warm start leaves saved snapshots untouched.
        evaluator.state = State(initial, [], None, np.array([], dtype=int))
        value = evaluator.evaluate(design)
        state = evaluator.state
        reference_power, reference_flux = _power_and_flux(domain, original_state)
        archive_reference_error = max(abs(reference_power / evaluator.power_reference - 1),
                                      abs(reference_flux[2] / reference_flux[0] - evaluator.side_reference))
        report["cases"][name] = case = {"beta": beta, "snapshot_sha256": source_hash,
                                       "power_reference": evaluator.power_reference,
                                       "side_reference": evaluator.side_reference, "directions": {},
                                       "reference_snapshot_sha256": original_hash,
                                       "reference_response_parity": float(archive_reference_error)}
        if name == "gray":
            report["filter_and_transfer"] = filter_transfer_checks(evaluator, matrix, fixed_rhs, design, rng)
        if name == "beta_8" or name == "gray":
            case["state_jacobian"] = state_jacobian_checks(
                evaluator, state, evaluator.transfer.forward(evaluator.density(design).resistance), rng)
        for direction_name in direction_names:
            result = directional_check(evaluator, matrix, design, value, state, directions[direction_name], STEP_SIZES)
            case["directions"][direction_name] = result
            print(f"{name} {direction_name}: derivative={result['derivative_passed']} "
                  f"Taylor={result['taylor_passed']}", flush=True)
            checkpoint(report_path, report)
    report["passed"] = all(case["reference_response_parity"] <= EQUIVALENCE_TOLERANCE
                           and case.get("state_jacobian", {"passed": True})["passed"]
                           and all(item["passed"] for item in case["directions"].values())
                           for case in report["cases"].values()) and report.get(
                               "filter_and_transfer", {"passed": True})["passed"]
    checkpoint(report_path, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshots", type=Path, help="Pilot root containing stage_* snapshot directories")
    parser.add_argument("report", type=Path, help="Machine-readable output JSON")
    parser.add_argument("--cases", nargs="+", default=CASES, choices=CASES)
    parser.add_argument("--directions", nargs="+", default=("random", "junction_cell", "collar_cell"),
                        choices=("random", "junction_cell", "collar_cell"))
    args = parser.parse_args()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    result = audit(args.snapshots, args.cases, args.directions, args.report)
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
