"""Verify the topology transfer, state responses and adjoints on Sol."""

import json

import numpy as np

from .density import DensityPipeline, DesignGrid, ReducedDensity
from .domain import build_domain
from .evaluator import Evaluator, Stage, baseline_reference


TRANSFER_TOLERANCE = 1e-12
EQUIVALENCE_TOLERANCE = 1e-10
RESIDUAL_LIMIT = 1e-10
MASS_LIMIT = 1e-9
RESIDUAL_DERIVATIVE_LIMIT = 1e-6
ADJOINT_LIMIT = 1e-4
TAYLOR_ORDER_LIMIT = 1.8
STEPS = (0.02, 0.01, 0.005, 0.0025, 0.000625, 0.00015625, 0.0000390625)
# Sol showed second-order truncation through STEPS[-1] at binary bounds and
# on the beta=0 MMA endpoint. Keep the original samples and add a finer sweep.
FINE_STEPS = STEPS + (0.00001953125, 0.000009765625, 0.0000048828125)


def _relative_error(actual, expected, floor=1.0):
    return float(np.max(np.abs(actual - expected))
                 / max(float(np.max(np.abs(actual))), float(np.max(np.abs(expected))), floor))


def _power_and_flux(domain, state):
    """Independently integrate port traces with five Gauss nodes and outward normals."""
    fem = domain.fem
    viscosity, density, height, reynolds = 0.001002, 998.2, 0.004, 50.0
    bulk_velocity = reynolds * viscosity / (density * height)
    points, weights = np.polynomial.legendre.leggauss(5)
    t = (points + 1.0) / 2.0
    velocity_basis = np.column_stack(((1.0-t)*(1.0-2.0*t), t*(2.0*t-1.0), 4.0*t*(1.0-t)))
    nodal_velocity = state[:2*fem.n_u].reshape(2, fem.n_u).T
    nodal_pressure = state[2*fem.n_u:]
    outflux, power_out = [], []
    for edges, axis, sign in zip(domain.ports, (0, 0, 1), (-1, 1, 1)):
        lengths = np.linalg.norm(fem.points[edges[:, 1]] - fem.points[edges[:, 0]], axis=1)
        line_weights = lengths[:, None] * weights / 2.0
        pressure = nodal_pressure[edges[:, 0], None] * (1.0-t) + nodal_pressure[edges[:, 1], None] * t
        velocity = np.einsum("qa,eai->eqi", velocity_basis, nodal_velocity[edges])
        normal_flux = line_weights * sign * velocity[..., axis]
        head = pressure + reynolds * np.sum(velocity**2, axis=2) / 2.0
        outflux.append(float(normal_flux.sum()))
        power_out.append(float((normal_flux * head).sum()))
    return -viscosity * bulk_velocity**2 * sum(power_out), np.array([-outflux[0], *outflux[1:]])


def _transfer_check(evaluator, rng):
    fem = evaluator.domain.fem
    indices = evaluator.domain.cell_index
    dimensional = rng.standard_normal(evaluator.pipeline.grid.size)
    cotangent = rng.standard_normal(fem.weights.shape)
    actual = evaluator.transfer.forward(dimensional)
    expected = np.zeros(len(indices))
    mask = indices >= 0
    expected[mask] = dimensional[indices[mask]] * (0.004**2 / 0.001002)
    # Independently tally all pairs of FEM element and quadrature index into cells.
    transpose_reference = np.zeros(len(dimensional))
    repeated_indices = np.repeat(indices[mask], fem.weights.shape[1])
    np.add.at(transpose_reference, repeated_indices, cotangent[mask].ravel() * (0.004**2 / 0.001002))
    transpose = evaluator.transfer.transpose(cotangent)
    left, right = float(np.sum(actual * cotangent)), float(dimensional @ transpose)
    errors = {"forward": _relative_error(actual[:, 0], expected),
              "transpose": _relative_error(transpose, transpose_reference),
              "dot_product": abs(left-right) / max(abs(left), abs(right), 1.0)}
    return {"errors": errors, "dot_products": [left, right],
            "passed": max(errors.values()) <= TRANSFER_TOLERANCE}


def _equivalence(full, reduced_evaluator, controls, seed):
    reduced = reduced_evaluator.reduced
    raw = reduced.raw(controls)
    full.state = reduced_evaluator.state = seed
    full_result = full.evaluate(raw)
    reduced_result = reduced_evaluator.evaluate(controls)
    full_density = full.density(raw)
    reduced_density = reduced_evaluator.density(controls)
    errors = {name: _relative_error(getattr(full_density, name), getattr(reduced_density, name))
              for name in ("raw", "filtered", "projected", "resistance")}
    errors.update(
        quadrature_resistance=_relative_error(full.transfer.forward(full_density.resistance),
                                              reduced_evaluator.transfer.forward(reduced_density.resistance)),
        state=_relative_error(full.state.values, reduced_evaluator.state.values),
        normalized_power=abs(full_result.objective - reduced_result.objective),
        side_fraction=abs(full_result.side_fraction - reduced_result.side_fraction),
        area_mm2=abs(full_result.area_mm2 - reduced_result.area_mm2),
        flux=_relative_error(full_result.flux, reduced_result.flux),
        constraints=_relative_error(full_result.constraints, reduced_result.constraints),
        gradient_chain=_relative_error(reduced_result.gradient, reduced.jacobian.T @ full_result.gradient),
        constraints_chain=_relative_error(reduced_result.jacobian, full_result.jacobian @ reduced.jacobian),
    )
    residual = max(full_result.nonlinear_residual, reduced_result.nonlinear_residual)
    mass = max(full_result.mass_imbalance, reduced_result.mass_imbalance)
    return {"max_nonlinear_residual": residual, "max_relative_mass_imbalance": mass,
            "power_w_per_m": full_result.power_w_per_m, "side_fraction": full_result.side_fraction,
            "projected_area_mm2": full_result.area_mm2, "errors": errors,
            "passed": bool(max(errors.values()) <= EQUIVALENCE_TOLERANCE
                           and residual <= RESIDUAL_LIMIT and mass <= MASS_LIMIT)}


def _analytic_residual_direction(evaluator, density, state, raw_direction):
    """Assemble the fixed-state derivative without the production reverse operations."""
    pipeline, stage, fem = evaluator.pipeline, evaluator.stage, evaluator.domain.fem
    filtered_direction = pipeline.filter(raw_direction)
    projection_slope = np.ones_like(filtered_direction)
    if stage.beta > 0.0:
        denominator = np.tanh(stage.beta * pipeline.threshold) + np.tanh(stage.beta * (1.0-pipeline.threshold))
        projection_slope = stage.beta / (denominator * np.cosh(stage.beta * (density.filtered-pipeline.threshold))**2)
    resistance_slope = -stage.alpha_max * stage.q * (stage.q + 1.0) / (stage.q + density.projected)**2
    resistance_direction = (0.004**2 / 0.001002) * resistance_slope * projection_slope * filtered_direction
    indices = evaluator.domain.cell_index
    element_direction = np.zeros(len(indices))
    mask = indices >= 0
    element_direction[mask] = resistance_direction[indices[mask]]

    velocity_at_nodes = state.values[:2*fem.n_u].reshape(2, fem.n_u).T
    velocity_at_quad = np.einsum("qa,eai->eqi", fem.basis, velocity_at_nodes[fem.velocity_dofs])
    weak_drag = np.einsum("e,eq,qa,eqi->eai", element_direction, fem.weights, fem.basis, velocity_at_quad)
    residual_direction = np.zeros(fem.size)
    np.add.at(residual_direction, fem.velocity_dofs, weak_drag[..., 0])
    np.add.at(residual_direction, fem.velocity_dofs + fem.n_u, weak_drag[..., 1])
    return residual_direction[state.free]


def _direction_check(evaluator, design, direction, *, at_bound):
    evaluation = evaluator.evaluate(design)
    base_state = evaluator.state
    domain = evaluator.domain
    fem = domain.fem
    reduced = evaluator.reduced
    raw_direction = direction if reduced is None else reduced.jacobian @ direction
    base_density = evaluator.density(design)
    residual_derivative = _analytic_residual_direction(evaluator, base_density, base_state, raw_direction)
    response_derivative = np.array([evaluation.gradient @ direction,
                                    -0.001 * evaluation.jacobian[0] @ direction,
                                    evaluation.jacobian[2] @ direction])
    constraint_derivative = evaluation.jacobian @ direction
    mass_errors, state_residuals = [], []

    def primal(candidate, *, solve=True):
        density = evaluator.density(candidate)
        resistance = evaluator.transfer.forward(density.resistance)
        operator = fem.viscous_pressure + fem.drag_matrix(resistance)
        residual = fem.residual(base_state.values, operator, evaluator.load, tangent=False)[0][base_state.free]
        state = evaluator.solve(resistance, base_state.values) if solve else base_state
        power, flux = _power_and_flux(domain, state.values)
        mass_errors.append(float(abs(flux[1:].sum() - flux[0]) / flux[0]))
        state_residuals.append(state.residual_history[-1])
        area = 41.0 + evaluator.pipeline.grid.cell_area_mm2 * float(np.sum(density.projected))
        return residual, np.array([power/evaluator.power_reference, flux[2]/flux[0], area])

    base_residual, base_response = primal(design, solve=False)
    response_error = np.abs(base_response - [evaluation.objective, evaluation.side_fraction, evaluation.area_mm2])
    split, area = base_response[1:]
    independent_constraints = np.array([(0.499-split)/0.001, (split-0.501)/0.001, area-83.0])
    constraint_error = _relative_error(evaluation.constraints, independent_constraints)
    records = []
    steps = FINE_STEPS if at_bound or evaluator.stage.beta == 0 else STEPS
    for step in steps:
        plus_residual, plus_response = primal(design + step * direction)
        if at_bound:
            double_residual, double_response = primal(design + 2.0 * step * direction)
            fd_residual = (4.0*(plus_residual-base_residual) - (double_residual-base_residual)) / (2.0*step)
            fd_response = (4.0*(plus_response-base_response) - (double_response-base_response)) / (2.0*step)
        else:
            minus_residual, minus_response = primal(design - step * direction)
            fd_residual = (plus_residual-minus_residual) / (2.0*step)
            fd_response = (plus_response-minus_response) / (2.0*step)
        scale = np.maximum(np.maximum(np.abs(response_derivative), np.abs(fd_response)), 1e-4)
        fd_constraints = np.array([-fd_response[1]/0.001, fd_response[1]/0.001, fd_response[2]])
        constraint_scale = np.maximum(np.maximum(np.abs(constraint_derivative), np.abs(fd_constraints)),
                                      np.array([0.1, 0.1, 1e-4]))
        records.append({"step": step, "response_fd": fd_response.tolist(),
                        "response_errors": (np.abs(response_derivative-fd_response) / scale).tolist(),
                        "constraint_errors": (np.abs(constraint_derivative-fd_constraints) / constraint_scale).tolist(),
                        "residual_error": _relative_error(residual_derivative, fd_residual, floor=1e-10),
                        "response_remainders": np.abs(plus_response-base_response-step*response_derivative).tolist()})

    floors = (5e-11, 5e-11, 1e-12)
    orders = []
    if not at_bound:
        # Use the four largest steps for Taylor order, ahead of the solve/roundoff floor.
        for broad, narrow in zip(records[:3], records[1:4]):
            orders.append([float(np.log(old/new) / np.log(broad["step"]/narrow["step"]))
                           if min(old, new) > floor else None
                           for old, new, floor in zip(broad["response_remainders"], narrow["response_remainders"], floors)])
    # At beta=0 the filtered area is affine, so its exact Taylor remainder
    # vanishes and a logarithmic order estimate would only measure roundoff.
    affine_area_at_floor = evaluator.stage.beta == 0 and all(
        row["response_remainders"][2] <= floors[2] for row in records
    )
    order_passed = at_bound or all(
        sum(value is not None and value >= TAYLOR_ORDER_LIMIT for value in metric) >= 2
        or (index == 2 and affine_area_at_floor)
        for index, metric in enumerate(zip(*orders))
    )
    derivative_passed = sum(row["residual_error"] <= RESIDUAL_DERIVATIVE_LIMIT
                            and max(row["response_errors"] + row["constraint_errors"]) <= ADJOINT_LIMIT
                            for row in records[-4:]) >= 2
    max_mass = max(mass_errors)
    max_residual = max(state_residuals)
    return {"difference": "second-order forward" if at_bound else "centered",
            "responses": ["normalized_power", "side_fraction", "projected_area_mm2"],
            "response_gradient_action": response_derivative.tolist(),
            "residual_derivative_inf_norm": float(np.linalg.norm(residual_derivative, ord=np.inf)),
            "base_measurement_errors": response_error.tolist(), "base_constraint_error": constraint_error,
            "records": records,
            "orders": orders, "affine_area_remainder_at_floor": bool(affine_area_at_floor),
            "max_relative_mass_imbalance": max_mass, "max_nonlinear_residual": max_residual,
            "passed": bool(order_passed and derivative_passed
                           and max(max(response_error), constraint_error) <= EQUIVALENCE_TOLERANCE
                           and max_mass <= MASS_LIMIT and max_residual <= RESIDUAL_LIMIT)}


def run_flow_verification():
    rng = np.random.default_rng(43202)
    grid = DesignGrid(0.5)
    domain = build_domain(grid, 0.5)
    pipeline = DensityPipeline(grid)
    reduced = ReducedDensity(grid)
    alpha_max = 1e5
    power_reference, initial = baseline_reference(domain, pipeline, alpha_max)
    power_independent, flux = _power_and_flux(domain, initial.values)
    baseline_power_error = abs(power_reference - power_independent) / power_reference
    baseline_stage = Stage(0.01, 0.0, alpha_max)
    projected_stage = Stage(0.05, 2.0, alpha_max)
    full_base = Evaluator(domain, pipeline, baseline_stage, power_reference, initial_state=initial)
    reduced_base = Evaluator(domain, pipeline, baseline_stage, power_reference, reduced=reduced, initial_state=initial)
    full = Evaluator(domain, pipeline, projected_stage, power_reference, initial_state=initial)
    reduced_evaluator = Evaluator(domain, pipeline, projected_stage, power_reference, reduced=reduced,
                                  initial_state=initial)

    zeros = np.zeros(reduced.jacobian.shape[1])
    baseline = _equivalence(full_base, reduced_base, zeros, initial)
    baseline["independent_power_relative_error"] = baseline_power_error
    baseline["flux_per_depth_dimensionless"] = flux.tolist()
    baseline["passed"] = bool(baseline["passed"] and baseline_power_error <= TRANSFER_TOLERANCE
                              and abs(baseline["projected_area_mm2"]-83.0) <= EQUIVALENCE_TOLERANCE
                              and np.array_equal(reduced.raw(zeros), grid.base_raw))
    equivalence = [_equivalence(full, reduced_evaluator, controls, initial)
                   for controls in (np.full_like(zeros, 0.4), rng.uniform(0.2, 0.8, len(zeros)),
                                    np.resize([0.2, 0.8], zeros.shape))]
    full_design = 0.16 + 0.68 * grid.base_raw + rng.uniform(-0.03, 0.03, grid.size)
    reduced_design = rng.uniform(0.3, 0.7, len(zeros))
    outlet_cell = np.flatnonzero(np.all(grid.centers_mm == (10.25, 5.75), axis=1)).item()
    outlet_direction = np.zeros(grid.size)
    outlet_direction[outlet_cell] = 1.0
    junction_direction = np.zeros_like(zeros)
    junction_direction[2] = 1.0  # Expansion centered at (11, 4.5) mm, beside the junction.
    transfer = _transfer_check(full, rng)
    derivatives = {
        "full_interior_random": _direction_check(full, full_design, rng.uniform(-0.6, 1.0, grid.size), at_bound=False),
        "full_interior_outlet_cell": _direction_check(full, full_design, outlet_direction, at_bound=False),
        "reduced_interior_random": _direction_check(reduced_evaluator, reduced_design,
                                                    rng.uniform(-0.6, 1.0, len(zeros)), at_bound=False),
        "reduced_interior_junction_control": _direction_check(reduced_evaluator, reduced_design, junction_direction,
                                                              at_bound=False),
        "full_bound_random": _direction_check(full_base, grid.base_raw,
                                               (1.0-2.0*grid.base_raw)*rng.uniform(0.1, 1.0, grid.size), at_bound=True),
        "reduced_bound_random": _direction_check(reduced_base, zeros, rng.uniform(0.1, 1.0, len(zeros)), at_bound=True),
    }
    passed = all(record["passed"] for record in [transfer, baseline, *equivalence, *derivatives.values()])
    return {"passed": bool(passed), "mesh": {"max_axis_spacing_mm": domain.max_spacing_mm,
                                             "triangles": len(domain.fem.triangles), "dofs": domain.fem.size},
            "design": {"spacing_mm": grid.spacing_mm, "cells": grid.size, "reduced_controls": len(zeros)},
            "stages": {"baseline": vars(baseline_stage), "projected": vars(projected_stage)},
            "thresholds": {"transfer": TRANSFER_TOLERANCE, "equivalence": EQUIVALENCE_TOLERANCE,
                           "nonlinear_residual": RESIDUAL_LIMIT, "relative_mass": MASS_LIMIT,
                           "residual_derivative": RESIDUAL_DERIVATIVE_LIMIT, "adjoint": ADJOINT_LIMIT,
                           "minimum_taylor_order": TAYLOR_ORDER_LIMIT},
            "transfer": transfer, "baseline": baseline, "equivalence": equivalence, "derivatives": derivatives}


if __name__ == "__main__":
    report = run_flow_verification()
    print(json.dumps(report, indent=2, allow_nan=False))
    raise SystemExit(0 if report["passed"] else 1)
