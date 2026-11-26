"""Independent density map, filtering, and differentiation checks for Sol."""

import json

import numpy as np

from .density import DensityPipeline, DesignGrid, ReducedDensity, brinkman, project


SPACINGS_MM = (0.5, 0.25, 0.125)
STATE_TOLERANCE = 1e-11
DERIVATIVE_TOLERANCE = 1e-6
MINIMUM_TAYLOR_ORDER = 1.8


def _scaled_error(candidate, reference):
    error = np.max(np.abs(np.asarray(candidate) - reference))
    return float(error / max(1., np.max(np.abs(reference))))


def _stencil_residual(grid, filtered, raw, radius):
    # Assemble the reference balance by adding equal and opposite face fluxes.
    field = filtered.reshape(grid.shape)
    residual = (filtered - raw).reshape(grid.shape) * grid.cell_area_mm2
    for axis in (0, 1):
        differences = radius ** 2 * np.diff(field, axis=axis)
        lower = [slice(None)] * 2
        upper = lower.copy()
        lower[axis] = slice(None, -1)
        upper[axis] = slice(1, None)
        residual[tuple(lower)] -= differences
        residual[tuple(upper)] += differences
    return float(np.max(np.abs(residual)))


def _basis_reference(grid):
    reference = np.zeros((grid.size, 20))
    for row, (x, y) in enumerate(grid.centers_mm):
        initially_fluid = y < 4 or 9.5 < x < 10.5
        heights, radius, offset, sign = ([1., 3.], 2.5, 10, -1.) if initially_fluid else ([4.5, 5.5], 2., 0, 1.)
        distances = (np.hypot(x - cx, y - cy) / radius for cy in heights for cx in (7., 9., 11., 13., 15.))
        group = np.array([(1. - r) ** 4 * (4. * r + 1.) if r < 1. else 0. for r in distances])
        reference[row, offset:offset + 10] = sign * group / max(1., group.sum())
    return reference


def _spatial_coverage(grid, jacobian, base):
    expansion = jacobian[:, :10].sum(axis=1)
    contraction = -jacobian[:, 10:].sum(axis=1)
    solid = base == 0.
    fluid = ~solid
    return {
        "array_order": "y outer, x inner, matching shape_ny_nx",
        "x_edges_mm": grid.x_edges_mm.tolist(),
        "y_edges_mm": grid.y_edges_mm.tolist(),
        "expansion_max_per_cell": expansion.tolist(),
        "contraction_max_per_cell": contraction.tolist(),
        "max_attainable_density_change_per_cell": np.abs(jacobian).sum(axis=1).tolist(),
        "expandable_solid_cells": int(np.count_nonzero(expansion[solid])),
        "contractible_fluid_cells": int(np.count_nonzero(contraction[fluid])),
        "eligible_solid_cells": int(np.count_nonzero(solid)),
        "eligible_fluid_cells": int(np.count_nonzero(fluid)),
    }


def _verify_grid(spacing, rng):
    grid = DesignGrid(spacing)
    reduced = ReducedDensity(grid)
    pipeline = DensityPipeline(grid)
    x, y = grid.centers_mm.T
    base = ((y < 4.) | ((y < 6.) & (x > 9.5) & (x < 10.5))).astype(float)
    baseline_exact = bool(np.array_equal(grid.base_raw, base) and np.array_equal(reduced.raw(np.zeros(20)), base))
    design_fluid = float(grid.cell_area_mm2 * base.sum())
    design_solid = float(grid.cell_area_mm2 * (1. - base).sum())
    reference_basis = _basis_reference(grid)
    basis_error = _scaled_error(reduced.jacobian, reference_basis)
    singular_values = np.linalg.svd(reduced.jacobian, compute_uv=False)
    rank = int(np.linalg.matrix_rank(reduced.jacobian))
    lower = base + np.minimum(reduced.jacobian, 0.).sum(axis=1)
    upper = base + np.maximum(reduced.jacobian, 0.).sum(axis=1)
    bounds_error = float(max(0., -lower.min(), upper.max() - 1.))
    control_direction = rng.uniform(0.1, 0.9, 20)
    steps = (1e-3, 1e-4, 1e-5)
    basis_fd_errors = [_scaled_error((reduced.raw(step * control_direction) - base) / step,
                                     reduced.jacobian @ control_direction) for step in steps]

    raw = rng.uniform(0., 1., grid.size)
    filtered = pipeline.filter(raw)
    reference_angle = np.pi * (x - 6.) / 10.
    mode = 0.5 + 0.25 * np.cos(reference_angle) * np.cos(np.pi * y / 6.)
    eigenvalue = (np.pi / 10.) ** 2 + (np.pi / 6.) ** 2
    mode_exact = 0.5 + (mode - 0.5) / (1. + pipeline.filter_radius_mm ** 2 * eigenvalue)
    mode_error = float(np.max(np.abs(pipeline.filter(mode) - mode_exact)))
    baseline = pipeline.evaluate(base, q=0.01, beta=0., alpha_max=1e6)
    test_covector = rng.normal(size=grid.size)
    transpose_error = _scaled_error(test_covector @ filtered, pipeline.filter_pullback(test_covector) @ raw)
    constant_error = _scaled_error(pipeline.filter(np.full(grid.size, 0.37)), 0.37)
    conservation_error = _scaled_error(grid.cell_area_mm2 * filtered.sum(), grid.cell_area_mm2 * raw.sum())
    filtered_extrema_error = float(max(0., filtered.max() - raw.max(), raw.min() - filtered.min()))
    stencil_error = _stencil_residual(grid, filtered, raw, pipeline.filter_radius_mm)
    beta0_exact = bool(np.array_equal(baseline.filtered, baseline.projected)
                       and np.array_equal(baseline.projection_derivative, np.ones(grid.size)))

    reduced_full_error = 0.
    for controls in (np.zeros(20), np.ones(20), rng.uniform(0., 1., 20)):
        candidate = reduced.raw(controls)
        full_raw = reference_basis @ controls + base
        full = pipeline.evaluate(full_raw, q=0.1, beta=4., alpha_max=1e6)
        mapped = pipeline.evaluate(candidate, q=0.1, beta=4., alpha_max=1e6)
        reduced_full_error = max(reduced_full_error, _scaled_error(candidate, full_raw),
                                 _scaled_error(mapped.filtered, full.filtered),
                                 _scaled_error(mapped.projected, full.projected),
                                 _scaled_error(mapped.resistance, full.resistance),
                                 _scaled_error(mapped.area_mm2, full.area_mm2))

    state_checks = {
        "basis_reference_error": basis_error,
        "box_bounds_error": bounds_error,
        "filter_constant_error": constant_error,
        "filter_conservation_error": conservation_error,
        "filter_extrema_error": filtered_extrema_error,
        "filter_stencil_residual": stencil_error,
        "filter_transpose_error": transpose_error,
        "reduced_full_density_equivalence_error": reduced_full_error,
        "baseline_area_error_mm2": abs(baseline.area_mm2 - 83.),
    }
    passed = (baseline_exact and beta0_exact and design_fluid == 42. and design_solid == 18.
              and rank == 20 and max(state_checks.values()) <= STATE_TOLERANCE
              and max(basis_fd_errors) <= DERIVATIVE_TOLERANCE)
    condition_number = float(singular_values[0] / singular_values[-1]) if singular_values[-1] > 0 else None
    return {"spacing_mm": spacing, "shape_ny_nx": list(grid.shape), "cells": grid.size,
            "baseline_exact": baseline_exact, "beta_zero_identity_exact": beta0_exact,
            "design_fluid_area_mm2": design_fluid, "design_solid_area_mm2": design_solid,
            "fixed_fluid_area_mm2": 41., "basis_rank": rank, "basis_singular_values": singular_values.tolist(),
            "basis_condition_number": condition_number,
            "basis_column_norms": np.linalg.norm(reduced.jacobian, axis=0).tolist(),
            "spatial_coverage": _spatial_coverage(grid, reduced.jacobian, base),
            "basis_one_sided_steps": steps, "basis_one_sided_derivative_errors": basis_fd_errors,
            "physical_neumann_mode_error": mode_error, "errors": state_checks, "passed": bool(passed)}


def _verify_derivatives(rng):
    grid = DesignGrid(0.5)
    pipeline = DensityPipeline(grid)
    values = np.array([0.005, 0.22, 0.5, 0.78, 0.995])
    direction = rng.uniform(0.1, 0.5, grid.size)
    raw = rng.uniform(0.2, 0.8, grid.size)
    weights = rng.uniform(0.1, 1., grid.size)
    weights /= weights.sum()
    alpha_max = 1e6
    perturbations = (1e-3, 1e-4, 1e-5, 1e-6)
    taylor_steps = np.array([0.01, 0.005, 0.0025])
    records = []
    for q, beta in ((0.01, 0.), (0.05, 2.), (0.10, 8.)):
        def response(z):
            density = pipeline.evaluate(z, q, beta, alpha_max)
            return np.array([weights @ density.resistance / alpha_max, density.area_mm2])

        projected, projection_derivative = project(values, beta, pipeline.threshold)
        _, resistance_derivative = brinkman(values, q, alpha_max)
        state = pipeline.evaluate(raw, q, beta, alpha_max)
        baseline = pipeline.evaluate(grid.base_raw, q, beta, alpha_max)
        bound_direction = (1. - 2. * grid.base_raw) * direction
        base_derivative = np.array([pipeline.pullback(baseline, weights / alpha_max) @ bound_direction,
                                    pipeline.area_gradient(baseline) @ bound_direction])
        interior_derivative = np.array([pipeline.pullback(state, weights / alpha_max) @ direction,
                                        pipeline.area_gradient(state) @ direction])
        base_response = response(grid.base_raw)
        interior_response = response(raw)
        errors = []
        for step in perturbations:
            fp = project(values + step, beta, pipeline.threshold)[0]
            fm = project(values - step, beta, pipeline.threshold)[0]
            ap = brinkman(values + step, q, alpha_max)[0]
            am = brinkman(values - step, q, alpha_max)[0]
            central = (response(raw + step * direction) - response(raw - step * direction)) / (2 * step)
            forward_step = response(grid.base_raw + step * bound_direction)
            forward_double = response(grid.base_raw + 2 * step * bound_direction)
            one_sided = (4 * forward_step - forward_double - 3 * base_response) / (2 * step)
            errors.append({"step": step,
                           "projection": _scaled_error(projection_derivative, (fp - fm) / (2 * step)),
                           "brinkman": _scaled_error(resistance_derivative, (ap - am) / (2 * step)),
                           "interior_chain": _scaled_error(interior_derivative, central),
                           "bound_chain": _scaled_error(base_derivative, one_sided)})
        remainders = np.array([abs((response(raw + step * direction) - interior_response
                                   - step * interior_derivative)[0]) for step in taylor_steps])
        order = np.log2(remainders[:-1] / remainders[1:])
        endpoints, endpoint_gradient = brinkman(np.array([0., 1.]), q, alpha_max)
        endpoint_error = _scaled_error(endpoints, np.array([alpha_max, 0.]))
        passed = (max(errors[-1][name] for name in ("projection", "brinkman", "interior_chain", "bound_chain"))
                  <= DERIVATIVE_TOLERANCE and min(order) >= MINIMUM_TAYLOR_ORDER
                  and np.all(np.isfinite(projected)) and endpoint_error <= STATE_TOLERANCE
                  and endpoint_gradient[-1] < 0)
        records.append({"q": q, "beta": beta, "errors": errors, "taylor_steps": taylor_steps.tolist(),
                        "taylor_remainders": remainders.tolist(), "taylor_orders": order.tolist(),
                        "brinkman_endpoint_error": endpoint_error,
                        "fully_fluid_resistance_derivative": float(endpoint_gradient[-1]), "passed": bool(passed)})
    return records


def run_density_verification():
    rng = np.random.default_rng(43201)
    grids = [_verify_grid(spacing, rng) for spacing in SPACINGS_MM]
    derivatives = _verify_derivatives(rng)
    mode_errors = np.array([grid["physical_neumann_mode_error"] for grid in grids])
    orders = np.log2(mode_errors[:-1] / mode_errors[1:])
    passed = all(item["passed"] for item in grids + derivatives) and min(orders) >= MINIMUM_TAYLOR_ORDER
    return {"passed": bool(passed), "grids": grids, "derivatives": derivatives,
            "physical_filter_refinement_orders": orders.tolist(),
            "limits": {"state_error": STATE_TOLERANCE, "derivative_error": DERIVATIVE_TOLERANCE,
                       "minimum_taylor_and_filter_order": MINIMUM_TAYLOR_ORDER}}


if __name__ == "__main__":
    report = run_density_verification()
    print(json.dumps(report, indent=2, allow_nan=False))
    raise SystemExit(0 if report["passed"] else 1)
