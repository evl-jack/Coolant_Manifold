"""Independent Sol verification for the full-density topology and baffle model."""

import argparse
import json
import os
from pathlib import Path
import traceback

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from shapely.geometry import Polygon, box

from .connectivity import ConnectivityGraph, RectangularBaffleClearance, run_baffle_connectivity_verification
from .density import DesignGrid
from .domain import build_domain
from .evaluator import Stage
from .hybrid_density import HybridDensityModel, HybridEvaluator, hybrid_reference
from .verify_flow import (
    ADJOINT_LIMIT, EQUIVALENCE_TOLERANCE, FINE_STEPS, MASS_LIMIT,
    RESIDUAL_DERIVATIVE_LIMIT, RESIDUAL_LIMIT, TAYLOR_ORDER_LIMIT,
    _power_and_flux, _relative_error,
)


HEIGHT_MM = 4.0
LAYER_MM = 0.125
UPPER_HEIGHT_MM = 3.75
RESISTANCE_SCALE = 0.004**2 / 0.001002


def _smooth_indicator(t):
    t = np.clip(t, 0.0, 1.0)
    return t**3 * (10.0 + t * (-15.0 + 6.0 * t))


def _smooth_indicator_derivative(t):
    t = np.clip(t, 0.0, 1.0)
    return 30.0 * (t * (1.0 - t))**2


def _baffle_area(height):
    """Baffle width is 1 mm. Integrate the mask in y using its primitive."""
    if height >= LAYER_MM:
        return height - 0.5 * LAYER_MM
    t = height / LAYER_MM
    return LAYER_MM * t**4 * (2.5 - 3.0 * t + t**2)


class IndependentHybrid:
    """Forward reference using a face-built filter and clipped triangle quadrature."""

    def __init__(self, domain, grid, stage):
        self.domain = domain
        self.grid = grid
        self.stage = stage
        self.fem = domain.fem
        self.load = self.fem.load()
        x, y = grid.centers_mm.T
        self.fixed = (13.0 < x) & (x < 14.0) & (0.0 < y) & (y < 4.0)
        self.design_elements = domain.cell_index >= 0
        self.design_cells = domain.cell_index[self.design_elements]
        element_ids = np.flatnonzero(self.design_elements)
        self.strip_elements = element_ids[self.fixed[self.design_cells]]
        self.strip_vertices = self.fem.points[self.fem.triangles[self.strip_elements]] * HEIGHT_MM
        self.strip_polygons = [Polygon(vertices) for vertices in self.strip_vertices]
        self.filter_factor = self._assemble_filter()
        nodes, weights = np.polynomial.legendre.leggauss(7)
        nodes, weights = (nodes + 1.0) / 2.0, weights / 2.0
        u, v = np.meshgrid(nodes, nodes, indexing="ij")
        self.subtriangle_coordinates = np.stack((1.0-u, u*(1.0-v), u*v), axis=2).reshape(-1, 3)
        self.subtriangle_weights = (u * np.outer(weights, weights)).ravel()
        self.strip_coordinates = np.linalg.inv(
            np.concatenate((np.ones((*self.strip_vertices.shape[:2], 1)), self.strip_vertices), axis=2)
        )

    def _assemble_filter(self):
        cells = np.arange(self.grid.size).reshape(self.grid.shape)
        first = np.r_[cells[:, :-1].ravel(), cells[:-1, :].ravel()]
        last = np.r_[cells[:, 1:].ravel(), cells[1:, :].ravel()]
        row = np.r_[first, last, first, last, cells.ravel()]
        col = np.r_[first, last, last, first, cells.ravel()]
        weights = np.r_[np.full(2 * first.size, 0.5**2),
                        np.full(2 * first.size, -0.5**2),
                        np.full(self.grid.size, self.grid.cell_area_mm2)]
        shape = (self.grid.size, self.grid.size)
        return splu(coo_matrix((weights, (row, col)), shape=shape).tocsc())

    def _filter(self, raw):
        return self.filter_factor.solve(self.grid.cell_area_mm2 * raw)

    def fields(self, design, direction=None):
        raw = np.ones(self.grid.size)
        raw[~self.fixed] = design[:-1]
        filtered = self._filter(raw)
        beta = self.stage.beta
        if beta == 0.0:
            projected = filtered.copy()
            slope = np.ones_like(projected)
        else:
            denominator = 2.0 * np.tanh(0.5 * beta)
            transition = beta * (filtered - 0.5)
            projected = 0.5 + np.tanh(transition) / denominator
            slope = beta / (denominator * np.cosh(transition)**2)
        projected[self.fixed] = 1.0
        slope[self.fixed] = 0.0
        if direction is None:
            return raw, filtered, projected
        raw_direction = np.zeros_like(raw)
        raw_direction[~self.fixed] = direction[:-1]
        return slope * self._filter(raw_direction)

    def _full_grid_mass(self, coefficient):
        dimensional_elements = np.zeros(len(self.fem.triangles))
        dimensional_elements[self.design_elements] = coefficient[self.design_cells]
        return self.fem.drag_matrix((dimensional_elements * RESISTANCE_SCALE)[:, None])

    def _baffle_mass(self, height, derivative=False):
        """Clip each region, triangulate the convex piece and integrate separately."""
        local = np.zeros((len(self.fem.triangles), 15, 15))
        shoulder = max(height - LAYER_MM, 0.0)
        regions = [(box(12.0, shoulder, 15.0, height), True)] if height > 0.0 else []
        if shoulder > 0.0 and not derivative:
            regions.append((box(12.0, 0.0, 15.0, shoulder), False))
        for index, polygon, inv in zip(self.strip_elements, self.strip_polygons, self.strip_coordinates):
            mass = np.zeros((6, 6))
            for region, transition in regions:
                piece = polygon.intersection(region)
                if piece.geom_type != "Polygon" or piece.is_empty:
                    continue
                vertices = np.asarray(piece.exterior.coords)[:-1]
                for left, right in zip(vertices[1:-1], vertices[2:]):
                    triangle = np.stack((vertices[0], left, right))
                    jacobian = abs(np.linalg.det(triangle[1:] - triangle[0])) / HEIGHT_MM**2
                    xy = self.subtriangle_coordinates @ triangle
                    coordinates = np.column_stack((np.ones(len(xy)), xy))
                    ell = coordinates @ inv
                    a, b, c = ell.T
                    p2 = np.column_stack((ell * (2.0 * ell - 1.0), 4.0 * a * b, 4.0 * b * c, 4.0 * c * a))
                    depth = (height - xy[:, 1]) / LAYER_MM
                    if derivative:
                        factor = _smooth_indicator_derivative(depth) * HEIGHT_MM / LAYER_MM
                    else:
                        factor = _smooth_indicator(depth) if transition else np.ones_like(depth)
                    weights = jacobian * self.subtriangle_weights * factor
                    mass += np.einsum("q,qa,qb->ab", weights, p2, p2)
            local[index, :6, :6] = local[index, 6:12, 6:12] = mass
        return self.stage.alpha_max * RESISTANCE_SCALE * self.fem.sparse(local)

    def operators(self, design):
        projected = self.fields(design)[2]
        q = self.stage.q
        resistance = self.stage.alpha_max * q * (1.0 - projected) / (q + projected)
        base_drag = self._full_grid_mass(resistance)
        height = HEIGHT_MM * design[-1]
        baffle_drag = self._baffle_mass(height)
        area = 41.0 + self.grid.cell_area_mm2 * projected.sum() - _baffle_area(height)
        return base_drag, baffle_drag, float(area)

    def derivative(self, design, direction):
        projected = self.fields(design)[2]
        projected_direction = self.fields(design, direction)
        q = self.stage.q
        resistance_direction = (-self.stage.alpha_max * q * (1.0 + q) / (q + projected)**2
                                * projected_direction)
        height = HEIGHT_MM * design[-1]
        height_area_derivative = -HEIGHT_MM * _smooth_indicator(height / LAYER_MM)
        drag = self._full_grid_mass(resistance_direction) + direction[-1] * self._baffle_mass(height, derivative=True)
        area = self.grid.cell_area_mm2 * projected_direction.sum() + direction[-1] * height_area_derivative
        return drag, float(area)

    def solve(self, drag, initial):
        return self.fem.solve(self.domain.fixed, self.domain.prescribed, self.load,
                              initial=initial, added_drag=drag)


def _sparse_error(candidate, expected):
    error = candidate - expected
    scale = max(abs(candidate).max(), abs(expected).max(), 1e-10)
    return float(abs(error).max() / scale)


def _taylor_check(steps, remainders, at_bound):
    floors = (5e-11, 5e-11, 1e-12, 5e-11)
    orders = []
    if at_bound:
        return True, orders, []
    for i in range(len(steps) - 1):
        orders.append([
            float(np.log(old / new) / np.log(steps[i] / steps[i + 1]))
            if min(old, new) > floor else None
            for old, new, floor in zip(remainders[i], remainders[i + 1], floors)
        ])
    at_floor = [bool(np.all(remainders[:, i] <= floor)) for i, floor in enumerate(floors)]

    def consecutive_order(metric):
        return any(left is not None and right is not None
                   and min(left, right) >= TAYLOR_ORDER_LIMIT
                   for left, right in zip(metric[:-1], metric[1:]))

    passed = all(
        at_floor[i] or consecutive_order(metric)
        for i, metric in enumerate(zip(*orders))
    )
    return passed, orders, at_floor


def _direction_check(evaluator, independent, design, direction, *, at_bound, height_cost):
    evaluation = evaluator.evaluate(design)
    state = evaluator.state
    free = state.free
    production_drag, production_height_drag, production_area, production_area_gradient = evaluator.operators(design)
    expected_base, expected_baffle, expected_area = independent.operators(design)
    expected_height_drag = independent._baffle_mass(HEIGHT_MM * design[-1], derivative=True)
    expected_drag_derivative, expected_area_derivative = independent.derivative(design, direction)
    residual_derivative = evaluator.residual_direction(design, state.values, direction)
    if residual_derivative.shape == state.values.shape:
        residual_derivative = residual_derivative[free]
    expected_residual_derivative = (expected_drag_derivative @ state.values)[free]
    density = evaluator.density(design)
    density_expected = independent.fields(design)
    derivative_errors = {
        "assembled_height_drag": _sparse_error(production_height_drag, expected_height_drag),
        "independent_residual_action": _relative_error(residual_derivative, expected_residual_derivative, floor=1e-10),
        "independent_area_action": _relative_error(production_area_gradient @ direction, expected_area_derivative,
                                                   floor=1e-4),
    }
    forward_errors = {
        "assembled_drag": _sparse_error(production_drag, expected_base + expected_baffle),
        "area_mm2": abs(production_area - expected_area),
        "evaluation_area_mm2": abs(evaluation.area_mm2 - expected_area),
        "raw": _relative_error(density.raw, density_expected[0]),
        "filtered": _relative_error(density.filtered, density_expected[1]),
        "projected": _relative_error(density.projected, density_expected[2]),
    }

    mass_errors, state_residuals = [], []

    def measure(candidate, primal, area):
        power, flux = _power_and_flux(evaluator.domain, primal.values)
        mass_errors.append(float(abs(flux[1:].sum() - flux[0]) / flux[0]))
        state_residuals.append(primal.residual_history[-1])
        normalized_power = power / evaluator.power_reference
        return np.array([normalized_power, flux[2] / flux[0], area,
                         normalized_power + height_cost * candidate[-1]])

    def primal(candidate):
        base, baffle, area = independent.operators(candidate)
        # Keep the two matrices separate when differencing to preserve the
        # tiny additions at the zero-height transition.
        response = measure(candidate, independent.solve(base + baffle, state.values), area)
        return base, baffle, response

    base_response = measure(design, state, expected_area)
    audit_state = independent.solve(expected_base + expected_baffle, state.values)
    forward_errors["independent_primal"] = float(
        np.max(np.abs(measure(design, audit_state, expected_area) - base_response)))
    optimizer_values = np.array([evaluation.power_w_per_m / evaluator.power_reference,
                                 evaluation.side_fraction, evaluation.area_mm2, evaluation.objective])
    forward_errors["independent_port_response"] = float(np.max(np.abs(base_response - optimizer_values)))
    split, area = base_response[1:3]
    constraints = np.array([(0.499-split)/0.001, (split-0.501)/0.001, area-83.0])
    forward_errors["independent_constraint"] = _relative_error(evaluation.constraints[:3], constraints)
    objective_derivative = evaluation.gradient @ direction
    response_derivative = np.array([objective_derivative - height_cost * direction[-1],
                                    -0.001 * evaluation.jacobian[0] @ direction,
                                    evaluation.jacobian[2] @ direction, objective_derivative])
    constraint_derivative = evaluation.jacobian[:3] @ direction

    lower = np.zeros(design.size)
    upper = np.r_[np.ones(design.size - 1), UPPER_HEIGHT_MM / HEIGHT_MM]
    outward = direction != 0.0
    inward_limit = np.full(design.size, np.inf)
    grows, shrinks = direction > 0.0, direction < 0.0
    inward_limit[grows] = (upper[grows] - design[grows]) / direction[grows]
    inward_limit[shrinks] = (lower[shrinks] - design[shrinks]) / direction[shrinks]
    if at_bound:
        active_limit = np.min(inward_limit[outward]) / 2.0
    else:
        central_limit = np.minimum(upper[outward] - design[outward], design[outward] - lower[outward])
        active_limit = np.min(central_limit / np.abs(direction[outward]))
    # Fill gaps in the original finite-difference sequence with exact halvings
    # so Taylor order can be checked locally without choosing disconnected steps.
    count = round(np.log2(FINE_STEPS[0] / FINE_STEPS[-1])) + 1
    zero_height_direction = design[-1] == 0.0 and not np.any(direction[:-1])
    if zero_height_direction:
        count += 4
    steps = min(FINE_STEPS[0], 0.8 * active_limit) / (2.0**np.arange(count))
    samples = []
    for step in steps:
        plus_base, plus_baffle, plus_response = primal(design + step * direction)
        if at_bound:
            double_base, double_baffle, double_response = primal(design + 2.0 * step * direction)
            delta = 4.0 * (plus_base - expected_base) - (double_base - expected_base)
            delta += 4.0 * (plus_baffle - expected_baffle) - (double_baffle - expected_baffle)
            fd_response = (4.0 * (plus_response - base_response) - (double_response - base_response)) / (2.0 * step)
        else:
            minus_base, minus_baffle, minus_response = primal(design - step * direction)
            delta = (plus_base - minus_base) + (plus_baffle - minus_baffle)
            fd_response = (plus_response - minus_response) / (2.0 * step)
        fd_residual = (delta @ state.values)[free] / (2.0 * step)
        scale = np.maximum(np.maximum(np.abs(response_derivative), np.abs(fd_response)), 1e-4)
        fd_constraints = np.array([-fd_response[1]/0.001, fd_response[1]/0.001, fd_response[2]])
        constraint_scale = np.maximum(np.maximum(np.abs(constraint_derivative), np.abs(fd_constraints)),
                                      [0.1, 0.1, 1e-4])
        samples.append({
            "step": float(step),
            "response_fd": fd_response.tolist(),
            "response_errors": (np.abs(response_derivative - fd_response) / scale).tolist(),
            "constraint_errors": (np.abs(constraint_derivative - fd_constraints) / constraint_scale).tolist(),
            "residual_error": _relative_error(residual_derivative, fd_residual, floor=1e-10),
            "response_remainders": np.abs(plus_response - base_response - step * response_derivative).tolist(),
        })
    remainders = np.array([sample["response_remainders"] for sample in samples])
    orders_passed, orders, at_floor = _taylor_check(steps, remainders, at_bound)
    finite_differences_passed = sum(
        sample["residual_error"] <= RESIDUAL_DERIVATIVE_LIMIT
        and max(sample["response_errors"] + sample["constraint_errors"]) <= ADJOINT_LIMIT
        for sample in samples[-4:]
    ) >= 2
    exact_zero_height_passed = not zero_height_direction or (
        abs(production_height_drag).max() == 0.0
        and not np.any(residual_derivative)
        and not np.any(response_derivative[:3])
        and response_derivative[3] == height_cost * direction[-1]
    )
    return {
        "height_mm": float(HEIGHT_MM * design[-1]),
        "difference": "second-order inward one-sided" if at_bound else "centered",
        "responses": ["normalized_power", "side_fraction", "area_mm2", "normalized_objective"],
        "independent_assembly_errors": forward_errors,
        "independent_derivative_errors": derivative_errors,
        "response_derivative": response_derivative.tolist(),
        "height_drag_derivative_inf_norm": float(abs(expected_height_drag).max()),
        "fixed_state_residual_derivative_inf_norm": float(np.linalg.norm(residual_derivative, ord=np.inf)),
        "closed_form_area_mm2": expected_area,
        "closed_form_area_derivative": expected_area_derivative,
        "samples": samples,
        "orders": orders,
        "taylor_remainders_at_roundoff_floor": at_floor,
        "exact_zero_height_passed": bool(exact_zero_height_passed) if zero_height_direction else None,
        "maximum_nonlinear_residual": max(state_residuals),
        "maximum_relative_mass_error": max(mass_errors),
        "passed": bool(max(forward_errors.values()) <= EQUIVALENCE_TOLERANCE
                       and max(derivative_errors.values()) <= RESIDUAL_DERIVATIVE_LIMIT
                       and finite_differences_passed and orders_passed and exact_zero_height_passed
                       and max(state_residuals) <= RESIDUAL_LIMIT and max(mass_errors) <= MASS_LIMIT),
    }


def _graph_chain_check(domain, model, stage, reference, seed, independent,
                       interior, topology_direction, height_direction):
    clearance = RectangularBaffleClearance()
    graph = ConnectivityGraph(model.grid, baffle=clearance)
    evaluator = HybridEvaluator(domain, model, stage, reference, initial_state=seed, graph=graph)
    design = interior.copy()
    design[-1] = (0.75 - clearance.margin_mm - 0.5 * clearance.transition_mm) / HEIGHT_MM
    result = evaluator.evaluate(design)
    projected = independent.fields(design)[2]
    direct = graph.evaluate(projected, baffle_height_mm=HEIGHT_MM * design[-1])
    forward_error = _relative_error(result.constraints[3:], direct.constraints / graph.minimum_conductance)

    def responses(candidate):
        projected = independent.fields(candidate)[2]
        return (graph.evaluate(projected, baffle_height_mm=HEIGHT_MM * candidate[-1]).constraints
                / graph.minimum_conductance)

    errors = {}
    for name, direction in (("height", height_direction), ("topology", topology_direction),
                            ("mixed", 0.4 * height_direction + topology_direction)):
        action = result.jacobian[3:] @ direction
        rows = []
        for step in (1e-4, 1e-5, 1e-6, 1e-7):
            finite_difference = (responses(design + step * direction) - responses(design - step * direction)) / (2 * step)
            rows.append({"step": step, "fd": finite_difference.tolist(),
                         "error": _relative_error(action, finite_difference, floor=0.1)})
        errors[name] = {"derivative": action.tolist(), "samples": rows,
                        "passed": sum(row["error"] <= ADJOINT_LIMIT for row in rows) >= 2}

    covector = np.sin(np.arange(model.grid.size) * 1.25)
    forward_action = independent.fields(design, topology_direction)
    reverse_covector = model.projected_pullback(evaluator.density(design), covector)
    dot_error = _relative_error(covector @ forward_action, reverse_covector @ topology_direction[:-1])
    height_norm = np.linalg.norm(result.jacobian[3:, -1], ord=np.inf)
    return {"passed": bool(forward_error <= EQUIVALENCE_TOLERANCE and dot_error <= EQUIVALENCE_TOLERANCE
                           and height_norm > 1e-8 and all(record["passed"] for record in errors.values())),
            "filter_projection_transpose_dot_error": dot_error,
            "graph_constraints_scaled_forward_error": forward_error,
            "graph_height_jacobian_inf_norm": float(height_norm), "directional_checks": errors}


def verify_hybrid():
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run hybrid verification in a Sol Slurm job")
    grid = DesignGrid(0.5)
    domain = build_domain(grid, 0.5)
    model = HybridDensityModel(grid)
    reference, seed = hybrid_reference(domain, model, 1e5)
    stages = {"initial": Stage(0.01, 0.0, 1e5), "projected": Stage(0.05, 4.0, 3e5)}
    references = {name: IndependentHybrid(domain, grid, stage) for name, stage in stages.items()}
    reserved_area = domain.fem.weights[references["initial"].strip_elements].sum() * HEIGHT_MM**2
    bounds = model.initial.copy()
    interior = bounds.copy()
    rng = np.random.default_rng(43265)
    interior[:-1] = 0.17 + 0.66 * bounds[:-1] + rng.uniform(-0.04, 0.04, bounds.size - 1)
    height_direction = np.zeros_like(bounds)
    height_direction[-1] = 1.0
    topology_direction = np.r_[rng.uniform(-0.6, 0.6, bounds.size - 1), 0.0]
    topology_inward = np.r_[(1.0 - 2.0 * bounds[:-1]) * rng.uniform(0.2, 0.7, bounds.size - 1), 0.0]
    tests = [
        ("zero_height", "initial", bounds, 0.0, height_direction, True, 0.0),
        ("near_zero_height", "projected", interior, 0.045, height_direction, False, 0.0),
        ("transition_starts_at_wall", "projected", interior, 0.125, height_direction, False, 0.1),
        ("tip_crosses_cell_face", "projected", interior, 0.5, height_direction, False, 0.0),
        ("transition_base_crosses_cell_face", "projected", interior, 0.625, height_direction, False, 0.1),
        ("interior_height", "initial", interior, 1.313, height_direction, False, 0.0),
        ("near_upper_height", "projected", interior, 3.74, height_direction, False, 0.0),
        ("upper_height_bound", "projected", interior, 3.75, -height_direction, True, 0.1),
        ("topology_binary_bounds", "initial", bounds, 0.0, topology_inward, True, 0.0),
        ("projected_topology_direction", "projected", interior, 1.313, topology_direction, False, 0.1),
        ("coupled_topology_and_height_direction", "projected", interior, 1.313,
         topology_direction + 0.4 * height_direction, False, 0.1),
    ]
    outputs = {}
    for name, stage, raw, height, direction, one_sided, height_cost in tests:
        design = raw.copy()
        design[-1] = height / HEIGHT_MM
        evaluator = HybridEvaluator(domain, model, stages[stage], reference, initial_state=seed, height_cost=height_cost)
        try:
            result = _direction_check(evaluator, references[stage], design, direction,
                                      at_bound=one_sided, height_cost=height_cost)
        except Exception:
            result = {"passed": False, "error": traceback.format_exc()}
        outputs[name] = {"stage": stage, "height_cost": height_cost, **result}
        print(f"{name}: {'pass' if result['passed'] else 'fail'}", flush=True)
    topology_layout_passed = (np.array_equal(model.free, ~references["initial"].fixed)
                              and np.array_equal(bounds[:-1], grid.base_raw[model.free])
                              and bounds[-1] == 0.0 and bounds.size == 225
                              and abs(reserved_area - 4.0) <= EQUIVALENCE_TOLERANCE)
    reference_power, reference_flux = _power_and_flux(domain, seed.values)
    reference_error = abs(reference_power / reference - 1.0)
    try:
        graph_geometry = run_baffle_connectivity_verification()
    except Exception:
        graph_geometry = {"passed": False, "error": traceback.format_exc()}
    try:
        graph_chain = _graph_chain_check(domain, model, stages["projected"], reference, seed,
                                         references["projected"], interior, topology_direction, height_direction)
    except Exception:
        graph_chain = {"passed": False, "error": traceback.format_exc()}
    print(f"graph geometry: {'pass' if graph_geometry['passed'] else 'fail'}", flush=True)
    print(f"graph chain: {'pass' if graph_chain['passed'] else 'fail'}", flush=True)
    return {
        "passed": bool(topology_layout_passed and reference_error <= EQUIVALENCE_TOLERANCE
                       and graph_geometry["passed"] and graph_chain["passed"]
                       and all(record["passed"] for record in outputs.values())),
        "sol_job_id": os.environ["SLURM_JOB_ID"],
        "fixed_fluid_strip_mm2": float(reserved_area),
        "topology_layout_passed": bool(topology_layout_passed),
        "topology_decision_variables": int(bounds.size - 1),
        "original_reference_power_W_per_m": float(reference_power),
        "original_reference_side_fraction": float(reference_flux[2] / reference_flux[0]),
        "original_reference_power_relative_error": reference_error,
        "limits": {"primal_equivalence": EQUIVALENCE_TOLERANCE, "residual_derivative": RESIDUAL_DERIVATIVE_LIMIT,
                   "total_sensitivity": ADJOINT_LIMIT, "state_residual": RESIDUAL_LIMIT,
                   "relative_mass": MASS_LIMIT, "taylor_order": TAYLOR_ORDER_LIMIT},
        "directions": outputs,
        "baffle_graph_geometry": graph_geometry,
        "baffle_graph_chain": graph_chain,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New JSON report path")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run hybrid verification in a Sol Slurm job")
    with args.output.open("x", encoding="utf-8") as report:
        result = verify_hybrid()
        json.dump(result, report, indent=2, allow_nan=False)
        report.write("\n")
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
