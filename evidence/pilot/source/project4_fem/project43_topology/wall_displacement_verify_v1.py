"""Sol-only independent anchored-filter and fluid derivative verification."""

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import traceback

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from ..fem import State
from .density import DesignGrid
from .domain import build_domain
from .evaluator import Stage
from .gate2 import environment, write_report
from .verify_density import DERIVATIVE_TOLERANCE, STATE_TOLERANCE
from .verify_flow import (
    ADJOINT_LIMIT, EQUIVALENCE_TOLERANCE, FINE_STEPS, MASS_LIMIT, RESIDUAL_DERIVATIVE_LIMIT,
    RESIDUAL_LIMIT, TAYLOR_ORDER_LIMIT, TRANSFER_TOLERANCE, _power_and_flux, _relative_error,
)
from .wall_displacement_density_v1 import WallDisplacementDensityV1
from .wall_displacement_fem_v1 import WallDisplacementEvaluatorV1


SPACING_MM = 0.25
FILTER_RADIUS_MM = 0.5
RESISTANCE_SCALE = 0.004**2 / 0.001002
ALPHA_MAX = 1e5
REFERENCE_PATH = "p43-relaxed-flow.rsOvla/prepare-2404994/sharp_fitted_baseline.json"


def _segment_distance(point, start, end):
    # Every specified wall is axis aligned. Minimize coordinate distance to
    # the closed bounding intervals of the finite segment, including its ends.
    separation = [max(min(a, b) - p, p - max(a, b), 0.0) for p, a, b in zip(point, start, end)]
    return float(np.hypot(*separation))


def _reference_mask(grid, cores=((0.0, 1.0, 20.0, 3.0), (9.5, 2.0, 10.5, 7.0))):
    walls = (((0.0, 0.0), (20.0, 0.0)), ((0.0, 4.0), (9.5, 4.0)),
             ((10.5, 4.0), (20.0, 4.0)), ((9.5, 4.0), (9.5, 7.0)),
             ((10.5, 4.0), (10.5, 7.0)))
    editable = np.zeros(grid.size, dtype=bool)
    wall_count = np.zeros(len(walls), dtype=int)
    for cell, (cx, cy) in enumerate(grid.centers_mm):
        half = SPACING_MM / 2.0
        left, bottom, right, top = cx-half, cy-half, cx+half, cy+half
        if not (6.5 <= left <= right <= 15.5 and 0.0 <= bottom <= top <= 5.5):
            continue
        if any(min(right, x_max) > max(left, x_min) and min(top, y_max) > max(bottom, y_min)
               for x_min, y_min, x_max, y_max in cores):
            continue
        corners = ((left, bottom), (left, top), (right, bottom), (right, top))
        matches = [max(_segment_distance(point, start, end) for point in corners) <= 1.0 for start, end in walls]
        wall_count += matches
        editable[cell] = any(matches)
    x, y = grid.centers_mm.T
    sharp = ((y < 4.0) | ((9.5 < x) & (x < 10.5))).astype(float)
    return editable, sharp, wall_count


class ReferenceDensity:
    """Separate face-by-face finite-volume stencil and fixed-anchor elimination."""

    def __init__(self, grid, *, protected_v2=False):
        self.grid = grid
        if protected_v2:
            self.editable, self.sharp, self.wall_counts = _reference_mask(
                grid, ((0., 0., 20., 3.75), (9.5, 2., 10.5, 7.)))
        else:
            self.editable, self.sharp, self.wall_counts = _reference_mask(grid)
        self.original_area = 41.0
        y_cells, x_cells = grid.shape
        editable_ids = np.flatnonzero(self.editable)
        control_id = np.full(grid.size, -1, dtype=int)
        control_id[editable_ids] = np.arange(len(editable_ids))
        self.fixed_rhs = np.zeros(len(editable_ids))
        diagonal = np.full(len(editable_ids), grid.cell_area_mm2)
        row, column, values = [], [], []
        coefficient = FILTER_RADIUS_MM**2
        for index, cell in enumerate(editable_ids):
            y, x = divmod(cell, x_cells)
            for ny, nx in ((y-1, x), (y+1, x), (y, x-1), (y, x+1)):
                if not (0 <= ny < y_cells and 0 <= nx < x_cells):
                    continue
                neighbor = ny*x_cells + nx
                diagonal[index] += coefficient
                if self.editable[neighbor]:
                    row.append(index)
                    column.append(control_id[neighbor])
                    values.append(-coefficient)
                else:
                    self.fixed_rhs[index] += coefficient * self.sharp[neighbor]
        indices = np.arange(len(editable_ids))
        matrix = coo_matrix((np.r_[values, diagonal], (np.r_[row, indices], np.r_[column, indices])),
                            shape=(len(indices), len(indices))).tocsc()
        self.matrix = matrix
        self.factor = splu(matrix)

    def raw(self, compact):
        raw = self.sharp.copy()
        raw[self.editable] = compact
        return raw

    def filter(self, compact):
        filtered = self.sharp.copy()
        mass = self.grid.cell_area_mm2
        filtered[self.editable] = self.factor.solve(mass * compact + self.fixed_rhs)
        return filtered

    def filter_direction(self, direction):
        action = np.zeros(self.grid.size)
        action[self.editable] = self.factor.solve(self.grid.cell_area_mm2 * direction)
        return action

    def fields(self, compact, stage):
        filtered = self.filter(compact)
        if stage.beta == 0.0:
            projected, projection_slope = filtered.copy(), np.ones_like(filtered)
        else:
            t = stage.beta * (filtered - 0.5)
            denominator = 2.0 * np.tanh(stage.beta * 0.5)
            projected = 0.5 + np.tanh(t) / denominator
            projection_slope = stage.beta / (denominator * np.cosh(t)**2)
        q = stage.q
        resistance = stage.alpha_max * ((1.0 + q) / (q + projected) - 1.0) * q
        resistance_slope = -stage.alpha_max * q * (1.0 + q) / (q + projected)**2
        return filtered, projected, resistance, projection_slope, resistance_slope


def _filter_checks(pipeline, reference, interior, initial, direction, rng):
    editable = reference.editable
    raw = pipeline.raw(interior)
    filtered = pipeline.filter(interior)
    reference_values = reference.filter(interior)
    action = reference.filter_direction(direction)
    covector = rng.normal(size=pipeline.grid.size)
    reverse = pipeline.filter_pullback(covector)
    reference_reverse = (reference.grid.cell_area_mm2
                         * reference.factor.solve(covector[editable], trans="T"))
    rhs = reference.grid.cell_area_mm2 * interior + reference.fixed_rhs
    exact = {
        "mask": bool(np.array_equal(pipeline.free, editable)),
        "initialization": bool(np.array_equal(pipeline.initial, initial)),
        "raw": bool(np.array_equal(raw, reference.raw(interior))),
        "filtered_anchors": bool(np.array_equal(filtered[~editable], reference.sharp[~editable])),
    }
    errors = {
        "independent_anchored_filter": _relative_error(filtered, reference_values),
        "independent_helmholtz_equation": _relative_error(reference.matrix @ filtered[editable], rhs),
        "independent_filter_transpose": _relative_error(reverse, reference_reverse),
        "forward_reverse_pairing": _relative_error(covector @ action, reverse @ direction),
    }
    anchor_only_cotangent = covector.copy()
    anchor_only_cotangent[editable] = 0.0
    exact["zero_anchor_only_transpose"] = not np.any(pipeline.filter_pullback(anchor_only_cotangent))
    steps = (1e-3, 1e-4, 1e-5, 1e-6)
    differences = []
    for step in steps:
        fd = (pipeline.filter(interior + step * direction) - pipeline.filter(interior - step * direction)) / (2*step)
        differences.append({"step": step, "field_derivative_error": _relative_error(fd, action),
                            "covector_pairing_error": _relative_error(covector @ fd, reverse @ direction)})
    return {"passed": (all(exact.values()) and max(errors.values()) <= STATE_TOLERANCE
                       and all(max(record["field_derivative_error"], record["covector_pairing_error"])
                               <= DERIVATIVE_TOLERANCE for record in differences[-2:])),
            "exact": exact, "errors": errors, "finite_differences": differences,
            "total_design_cells": pipeline.grid.size, "editable_cells": pipeline.size,
            "fixed_fluid_cells": int(np.count_nonzero(reference.sharp[~editable])),
            "fixed_solid_cells": int(np.count_nonzero(1.0 - reference.sharp[~editable])),
            "editable_initial_fluid_cells": int(np.count_nonzero(initial)),
            "editable_initial_solid_cells": int(np.count_nonzero(1.0 - initial)),
            "editable_cells_per_wall_with_overlap": reference.wall_counts.tolist()}


class ReferenceFlow:
    """Reference transfer and state solves on the same verified P2/P1 operators."""

    def __init__(self, domain, reference, pde):
        self.domain = domain
        self.density = reference
        self.pde = pde
        self.fem = domain.fem
        self.load = self.fem.load()
        self.free = np.setdiff1d(np.arange(self.fem.size), domain.fixed)
        centroid = self.fem.points[self.fem.triangles].mean(axis=1) * 4.0
        x, y = centroid.T
        grid = reference.grid
        inside = ((grid.x_edges_mm[0] < x) & (x < grid.x_edges_mm[-1])
                  & (grid.y_edges_mm[0] < y) & (y < grid.y_edges_mm[-1]))
        cell_x = ((x[inside] - grid.x_edges_mm[0]) / SPACING_MM).astype(int)
        cell_y = ((y[inside] - grid.y_edges_mm[0]) / SPACING_MM).astype(int)
        self.in_domain = inside
        self.cell = cell_y*grid.shape[1] + cell_x

    def resistance(self, values):
        element = np.zeros(len(self.fem.triangles))
        element[self.in_domain] = values[self.cell] * RESISTANCE_SCALE
        return element[:, None]

    def solve(self, resistance, initial):
        if self.pde == "stokes":
            # Assemble the linear momentum, continuity and drag operator
            # explicitly. Apply the prescribed values by block elimination.
            operator = self.fem.viscous_pressure + self.fem.drag_matrix(resistance)
            free, fixed = self.free, self.domain.fixed
            solution = np.zeros(self.fem.size)
            solution[fixed] = self.domain.prescribed
            rhs = self.load[free] - operator[free][:, fixed] @ solution[fixed]
            solution[free] = splu(operator[free][:, free].tocsc()).solve(rhs)
            error = float(np.linalg.norm((operator @ solution - self.load)[free], ord=np.inf))
            return State(solution, [error], operator, free)
        return self.fem.solve(self.domain.fixed, self.domain.prescribed, self.load,
                              resistance=resistance, initial=initial)

    def continuation(self, resistance):
        if self.pde == "stokes":
            return self.solve(resistance, None)
        state = None
        try:
            for inertia in (0.0, 25.0, 50.0):
                self.fem.inertia = inertia
                state = self.solve(resistance, None if state is None else state.values)
        finally:
            self.fem.inertia = 50.0
        return state

    def residual_direction(self, state, resistance_direction):
        fem = self.fem
        local_velocity = state.values[:2*fem.n_u].reshape(2, fem.n_u).T[fem.velocity_dofs]
        velocity = np.einsum("qa,eai->eqi", fem.basis, local_velocity)
        weighted_resistance = self.resistance(resistance_direction) * fem.weights
        local_force = np.einsum("eq,qa,eqi->eai", weighted_resistance, fem.basis, velocity)
        result = np.zeros(fem.size)
        np.add.at(result, fem.velocity_dofs, local_force[..., 0])
        np.add.at(result, fem.velocity_dofs + fem.n_u, local_force[..., 1])
        return result[state.free]

    def transfer_check(self, evaluator, rng):
        grid = self.density.grid
        values = rng.normal(size=grid.size)
        weights = rng.normal(size=self.fem.weights.shape)
        expected_forward = self.resistance(values)
        expected_reverse = np.zeros(grid.size)
        np.add.at(expected_reverse, self.cell, RESISTANCE_SCALE * weights[self.in_domain].sum(axis=1))
        actual_forward = evaluator.transfer.forward(values)
        actual_reverse = evaluator.transfer.transpose(weights)
        errors = {
            "forward": _relative_error(actual_forward, expected_forward),
            "transpose": _relative_error(actual_reverse, expected_reverse),
            "dot_product": _relative_error(np.sum(actual_forward * weights), values @ actual_reverse),
        }
        return {"passed": max(errors.values()) <= TRANSFER_TOLERANCE, "errors": errors}


def _orders(samples, at_bound):
    if at_bound:
        return [], True
    floors = (5e-11, 5e-11)
    remainders = np.array([sample["response_taylor_remainders"] for sample in samples])
    orders = []
    for first, second in zip(remainders[:-1], remainders[1:]):
        orders.append([float(np.log2(a/b)) if min(a, b) > floor else None
                       for a, b, floor in zip(first, second, floors)])
    passed = True
    for column, floor in enumerate(floors):
        if np.all(remainders[:, column] <= floor):
            continue
        orders_column = [row[column] for row in orders]
        passed &= any(a is not None and b is not None and min(a, b) >= TAYLOR_ORDER_LIMIT
                      for a, b in zip(orders_column[:-1], orders_column[1:]))
    return orders, passed


def _direction_check(evaluator, flow, design, direction, *, at_bound):
    ref, stage, fem = flow.density, evaluator.stage, flow.fem
    reference_power = evaluator.power_reference
    response = evaluator.evaluate(design)
    state = evaluator.state
    density = evaluator.density(design)
    filtered, projected, resistance, projection_slope, resistance_slope = ref.fields(design, stage)
    expected_area = ref.original_area + ref.grid.cell_area_mm2 * projected.sum()
    errors = {
        "raw": _relative_error(density.raw, ref.raw(design)),
        "filtered": _relative_error(density.filtered, filtered),
        "projected": _relative_error(density.projected, projected),
        "resistance": _relative_error(density.resistance, resistance),
        "density_area_mm2": abs(density.area_mm2 - expected_area),
    }
    exact = {"filtered_anchor_values": bool(np.array_equal(density.filtered[~ref.editable], ref.sharp[~ref.editable])),
             "projected_anchor_values": bool(np.array_equal(density.projected[~ref.editable], ref.sharp[~ref.editable]))}
    filter_direction = ref.filter_direction(direction)
    resistance_direction = resistance_slope * projection_slope * filter_direction
    production_action = np.zeros(ref.grid.size)
    editable = ref.editable
    raw_mass = ref.grid.cell_area_mm2
    # The production transpose returns compact controls. Apply it to full
    # test rows via the archived forward/dual identity, while assembling the
    # fixed-state residual in the independent element loop.
    compact_slope = np.asarray(evaluator.pipeline.filter_pullback(projection_slope * resistance_slope))
    errors["filter_to_resistance_pairing"] = _relative_error(resistance_direction.sum(), compact_slope @ direction)
    production_action[editable] = (
        evaluator.pipeline.filter(design + (1e-6 if at_bound else 5e-7)*direction)[editable]
        - evaluator.pipeline.filter(design if at_bound else design - 5e-7*direction)[editable]
    ) / 1e-6
    filter_error = _relative_error(production_action, filter_direction)
    residual_action = flow.residual_direction(state, resistance_direction)
    independent_area_action = raw_mass * np.dot(projection_slope, filter_direction)
    errors["projected_area_pullback"] = _relative_error(
        evaluator.pipeline.area_gradient(density) @ direction, independent_area_action)
    constraint_action = response.jacobian @ direction
    adjoint = np.array([response.gradient @ direction, -0.001 * constraint_action[0]])
    mass_errors, state_residuals = [], []

    def measurements(primal):
        power, flux = _power_and_flux(flow.domain, primal.values)
        state_residuals.append(primal.residual_history[-1])
        mass_errors.append(abs(flux[1:].sum() - flux[0]) / flux[0])
        return np.array([power/reference_power, flux[2]/flux[0]])

    base_resistance = flow.resistance(resistance)
    production_base_resistance = evaluator.transfer.forward(density.resistance)
    base_operator = fem.viscous_pressure + fem.drag_matrix(base_resistance)
    if flow.pde == "stokes":
        base_residual = (base_operator @ state.values - flow.load)[state.free]
        tangent_error = state.tangent - base_operator
        errors["independent_stokes_tangent"] = float(
            abs(tangent_error).max() / max(1.0, abs(base_operator).max()))
        exact["zero_convective_inertia"] = fem.inertia == 0.0
    else:
        base_residual = fem.residual(state.values, base_operator, flow.load, tangent=False)[0][state.free]
    base_value = measurements(state)
    fraction = base_value[1]
    expected_constraints = np.array([(.10-fraction)/.001, (fraction-.90)/.001])
    errors["flow_constraints"] = _relative_error(response.constraints, expected_constraints)
    separate_base = flow.solve(base_resistance, state.values)
    errors["power_and_split"] = float(np.max(np.abs(base_value - [response.objective, response.side_fraction])))
    errors["independent_primal_replay"] = float(np.max(np.abs(base_value - measurements(separate_base))))
    errors["fixed_state_residual_equivalence"] = float(np.max(np.abs(base_residual)))

    def primal(compact):
        independent_resistance = flow.resistance(ref.fields(compact, stage)[2])
        production_resistance = evaluator.transfer.forward(evaluator.density(compact).resistance)
        # Subtract before assembly so the FD is insensitive to cancellation
        # of the design-independent Navier-Stokes and pressure residual.
        drag_delta = fem.drag_matrix(production_resistance - production_base_resistance)
        residual_delta = (drag_delta @ state.values)[state.free]
        return residual_delta, measurements(flow.solve(independent_resistance, state.values))

    count = round(np.log2(FINE_STEPS[0] / FINE_STEPS[-1])) + 1
    steps = FINE_STEPS[0] / 2.0**np.arange(count)
    samples = []
    for step in steps:
        plus_residual, plus_value = primal(design + step * direction)
        if at_bound:
            double_residual, double_value = primal(design + 2.0 * step * direction)
            fd_residual = (4.0*plus_residual - double_residual)/(2.0*step)
            fd_response = (4.0*(plus_value-base_value) - (double_value-base_value))/(2.0*step)
        else:
            minus_residual, minus_value = primal(design - step * direction)
            fd_residual = (plus_residual-minus_residual)/(2.0*step)
            fd_response = (plus_value-minus_value)/(2.0*step)
        response_scale = np.maximum(np.maximum(np.abs(adjoint), np.abs(fd_response)), 1e-4)
        fd_constraints = np.array([-fd_response[1], fd_response[1]]) / .001
        constraint_scale = np.maximum(np.maximum(np.abs(constraint_action), np.abs(fd_constraints)), 0.1)
        samples.append({
            "step": float(step), "residual_derivative_error": _relative_error(fd_residual, residual_action, floor=1e-10),
            "response_errors": (np.abs(fd_response-adjoint)/response_scale).tolist(),
            "constraint_errors": (np.abs(fd_constraints-constraint_action)/constraint_scale).tolist(),
            "response_fd": fd_response.tolist(),
            "response_taylor_remainders": np.abs(plus_value-base_value-step*adjoint).tolist(),
        })
    orders, order_passed = _orders(samples, at_bound)
    finite_difference_passed = sum(item["residual_derivative_error"] <= RESIDUAL_DERIVATIVE_LIMIT
                                   and max(item["response_errors"] + item["constraint_errors"]) <= ADJOINT_LIMIT
                                   for item in samples[-4:]) >= 2
    return {
        "passed": bool(all(exact.values()) and max(errors.values()) <= EQUIVALENCE_TOLERANCE
                       and filter_error <= DERIVATIVE_TOLERANCE and finite_difference_passed and order_passed
                       and max(mass_errors) <= MASS_LIMIT and max(state_residuals) <= RESIDUAL_LIMIT),
        "stencil": "second-order inward" if at_bound else "centered",
        "exact_anchor_checks": exact, "forward_errors": errors,
        "filter_direction_error": filter_error, "responses": ["normalized_power", "side_fraction"],
        "analytic_response_direction": adjoint.tolist(), "finite_differences": samples,
        "maximum_relative_mass_error": float(max(mass_errors)), "maximum_nonlinear_residual": max(state_residuals),
        "fixed_state_residual_direction_inf_norm": float(np.linalg.norm(residual_action, np.inf)),
        "taylor_orders": orders, "taylor_passed": order_passed, "finite_differences_passed": finite_difference_passed,
    }


def _reference_power(path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not payload["completed"]:
        raise ValueError(f"Incomplete fitted reference: {path}")
    if not payload["results"]["fitted_diagnostic_passed"]:
        raise ValueError(f"Fitted reference did not pass the mesh and state screen: {path}")
    levels = payload["results"]["levels"]
    power = float(levels[-1]["mechanical_power_W_per_m"])
    if not np.isfinite(power) or power <= 0.0:
        raise ValueError(f"Invalid fitted reference power: {power}")
    return power, levels


def verify(reference_path, pde="navier-stokes", *, protected_v2=False):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run this numerical verification in a Sol Slurm job")
    reference_power, fitted_levels = _reference_power(reference_path)
    rng = np.random.default_rng(43579)
    grid = DesignGrid(SPACING_MM)
    domain = build_domain(grid, SPACING_MM)
    if protected_v2:
        from .wall_displacement_smooth_v2 import ProtectedDensity

        pipeline = ProtectedDensity(grid)
    else:
        pipeline = WallDisplacementDensityV1(grid)
    ref = ReferenceDensity(grid, protected_v2=protected_v2)
    editable = ref.editable
    initial = ref.sharp[editable]
    interior = 0.14 + 0.72*initial + rng.uniform(-0.04, 0.04, len(initial))
    center_direction = rng.uniform(-0.8, 0.8, len(initial))
    bound_direction = (1.0 - 2.0*initial) * rng.uniform(0.2, 0.8, len(initial))
    filter_check = _filter_checks(pipeline, ref, interior, initial, center_direction, rng)
    flow = ReferenceFlow(domain, ref, pde)
    base_stage = Stage(q=0.01, beta=0.0, alpha_max=ALPHA_MAX)
    seed = flow.continuation(flow.resistance(ref.fields(initial, base_stage)[2]))
    evaluator_options = {"pde": "stokes"} if pde == "stokes" else {}
    evaluator = WallDisplacementEvaluatorV1(domain, pipeline, base_stage, reference_power,
                                            initial_state=seed, **evaluator_options)
    transfer_check = flow.transfer_check(evaluator, rng)
    checks = {}
    for q in ((.10,) if protected_v2 else (.01, .10)):
        for beta in ((8.,) if protected_v2 else (0., 8.)):
            stage = Stage(q=q, beta=beta, alpha_max=ALPHA_MAX)
            for name, design, direction, at_bound in (
                ("interior", interior, center_direction, False),
                ("binary_bounds", initial, bound_direction, True),
            ):
                label = f"q{q:.2f}_beta{beta:g}_{name}"
                try:
                    evaluator = WallDisplacementEvaluatorV1(domain, pipeline, stage, reference_power,
                                                            initial_state=seed, **evaluator_options)
                    record = _direction_check(evaluator, flow, design, direction, at_bound=at_bound)
                except Exception:
                    record = {"passed": False, "error": traceback.format_exc()}
                checks[label] = {"stage": asdict(stage), **record}
                print(f"{label}: {'pass' if record['passed'] else 'fail'}", flush=True)
    return {
        "passed": bool(filter_check["passed"] and transfer_check["passed"]
                       and all(check["passed"] for check in checks.values())),
        "sol_job_id": os.environ["SLURM_JOB_ID"],
        "environment": environment(),
        "pde": pde,
        "convective_inertia_coefficient": 0.0 if pde == "stokes" else 50.0,
        "port_power_uses_original_kinetic_head": True,
        "mesh": {"design_spacing_mm": SPACING_MM, "fem_spacing_mm": SPACING_MM,
                 "triangles": len(domain.fem.triangles), "dofs": domain.fem.size, "reynolds": 50.0},
        "reference_power_W_per_m": reference_power, "reference_source": str(reference_path),
        "reference_fitted_levels": fitted_levels,
        "filter": filter_check, "transfer": transfer_check, "direction_checks": checks,
        "thresholds": {"filter_equivalence": STATE_TOLERANCE, "filter_derivative": DERIVATIVE_TOLERANCE,
                       "transfer": TRANSFER_TOLERANCE, "equivalence": EQUIVALENCE_TOLERANCE,
                       "nonlinear_residual": RESIDUAL_LIMIT, "relative_mass": MASS_LIMIT,
                       "residual_derivative": RESIDUAL_DERIVATIVE_LIMIT, "adjoint": ADJOINT_LIMIT,
                       "minimum_taylor_order": TAYLOR_ORDER_LIMIT},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New report directory")
    parser.add_argument("--pde", choices=("navier-stokes", "stokes"), default="navier-stokes")
    parser.add_argument("--reference", type=Path,
                        default=Path("/scratch") / os.environ.get("USER", "") / REFERENCE_PATH)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical verification in a Sol Slurm job")
    args.output.mkdir(parents=True, exist_ok=False)
    result = verify(args.reference, args.pde)
    completed = write_report(args.output, "wall_displacement_verification_v1", lambda: result)
    raise SystemExit(0 if completed and result["passed"] else 1)


if __name__ == "__main__":
    main()
