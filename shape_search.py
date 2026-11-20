"""Search connected pipe shapes under the manifold's volume and outlet constraints."""

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.ndimage import label
from scipy.optimize import LinearConstraint, minimize

from .case import CASE
from .optimize import flow_metrics, outlet_violation, setup
from .solve import solve
from mesh import manifold_grid_axes


NAMES = ("upstream_pinch", "downstream_pinch", "branch_flare")
BOUNDS = np.array([(0., 1.), (0., 1.), (0., .5)])
MAIN_X = np.array(CASE["main_channel"]["x"]) * 1e3
MAIN_Y = np.array(CASE["main_channel"]["y"]) * 1e3
BRANCH_X = np.array(CASE["side_branch"]["x"]) * 1e3
ROOT_X = np.array(CASE["branch_root"]["x"]) * 1e3
ROOT_TOP = CASE["branch_root"]["y"][1] * 1e3
SIDE_TOP = CASE["side_branch"]["y"][1] * 1e3
STRAIGHT_START = SIDE_TOP - CASE["fixed_side_run_m"] * 1e3
MAIN_RUN = CASE["fixed_main_run_m"] * 1e3
UPSTREAM = (MAIN_X[0] + MAIN_RUN, ROOT_X[0])
DOWNSTREAM = (ROOT_X[1], MAIN_X[1] - MAIN_RUN)
AREA_GRADIENT = np.r_[-np.diff(UPSTREAM), -np.diff(DOWNSTREAM), STRAIGHT_START - ROOT_TOP]


def bump(coordinates, start, stop):
    position = (coordinates - start) / (stop - start)
    return np.where((0 < position) & (position < 1), np.sin(np.pi * position) ** 2, 0)


def shape_density(controls, spacing):
    """Cell coverage from a single open interval through each passage."""
    x, y = (edges * 1e3 for edges in manifold_grid_axes(spacing))
    xc, yc = (x[:-1] + x[1:]) / 2, (y[:-1] + y[1:]) / 2
    upstream, downstream, flare = controls
    pinch = upstream * bump(xc, *UPSTREAM) + downstream * bump(xc, *DOWNSTREAM)
    widening = flare * bump(yc, ROOT_TOP, STRAIGHT_START)
    widening += np.where((MAIN_Y[1] <= yc) & (yc < ROOT_TOP), BRANCH_X[0] - ROOT_X[0], 0)

    def coverage(edges, lower, upper):
        return np.clip(np.minimum(edges[1:], upper) - np.maximum(edges[:-1], lower), 0,
                       np.diff(edges)) / np.diff(edges)

    main = coverage(y, MAIN_Y[0] + pinch[:, None], MAIN_Y[1] - pinch[:, None])
    branch = coverage(x, BRANCH_X[0] - widening[:, None], BRANCH_X[1] + widening[:, None]).T
    branch[:, yc < MAIN_Y[1]] = 0
    return main + branch


def solve_walls(controls, spacing, iterations=1600):
    mask = shape_density(controls, spacing) >= .5
    if label(mask)[1] != 1:
        raise RuntimeError("Binary pipe geometry is disconnected")
    run = solve(iterations, model="explicit", include_side_outlet=True, spacing_m=spacing,
                fluid_mask=mask, convergence_tolerance=1e-4, verbose=False)
    return mask, run


def search(max_steps=5, iterations=1600):
    target, tolerance = CASE["side_to_main_ratio"], CASE["outlet_balance_tolerance"]
    spacing = CASE["spacing_m"]
    _, areas, _, _, initial, _, _, editable = setup()
    areas = areas.reshape(initial.shape) * 1e6
    budget = CASE["design_volume_fraction"] * areas[editable].sum()
    fixed_fluid = (areas[~editable] * initial[~editable]).sum()
    area_cap = fixed_fluid + budget
    reference_area = (areas * initial).sum()
    sampled_gradient = [((areas * shape_density(control, spacing)).sum() - reference_area) * 2
                        for control in .5 * np.eye(3)]
    area_constraint = LinearConstraint([AREA_GRADIENT, sampled_gradient], -np.inf,
                                       area_cap - reference_area)
    initial_controls = np.array([.40, .40, .10])
    records, cache = [], {}

    def evaluate(controls, state=None):
        key = tuple(np.round(controls, 9))
        if key in cache:
            return cache[key]
        density = shape_density(controls, spacing)
        # Verify the search geometry on the same grid used by Model B.
        fixed_ok = np.allclose(density[~editable], initial[~editable], atol=1e-10)
        volume_ok = fixed_ok and (areas[editable] * density[editable]).sum() <= budget + 1e-6
        run = solve(iterations, model="B", design_density=density, alpha_max=1e6,
                    convergence_tolerance=1e-4, initial_state=state, verbose=False)
        metrics = flow_metrics(run, target)
        metrics.update(zip(NAMES, controls))
        metrics.update(fluid_area_mm2=float((areas * density).sum()),
                       analytic_area_mm2=float(reference_area + AREA_GRADIENT @ controls),
                       volume_ok=volume_ok, fixed_geometry_ok=fixed_ok,
                       converged=run["converged"], iterations=run["iterations"])
        records.append(metrics)
        cache[key] = metrics, run
        return cache[key]

    def sensitivities(controls, h):
        _, state = evaluate(controls)
        jacobian = []
        for axis, (lower, upper) in enumerate(BOUNDS):
            left, right = controls.copy(), controls.copy()
            left[axis] = max(lower, controls[axis] - h)
            right[axis] = min(upper, controls[axis] + h)
            low, _ = evaluate(left, state)
            high, _ = evaluate(right, state)
            if not (low["converged"] and high["converged"]):
                raise RuntimeError(f"Incomplete flow convergence for {NAMES[axis]} sensitivity")
            jacobian.append([(high[key] - low[key]) / (right[axis] - left[axis])
                             for key in ("side_split", "outlet_balance_error", "pressure_drop_pa")])
        return np.array(jacobian).T

    def metrics_at(controls):
        metrics, _ = evaluate(controls)
        if not metrics["converged"]:
            raise RuntimeError("A search candidate did not converge")
        return metrics

    def split_objective(controls):
        error = metrics_at(controls)["outlet_balance_error"]
        return .5 * error ** 2

    def split_gradient(controls):
        return metrics_at(controls)["outlet_balance_error"] * sensitivities(controls, .08)[1]

    start = metrics_at(initial_controls)
    coarse_grad = sensitivities(initial_controls, .08)
    fine_grad = sensitivities(initial_controls, .04)
    print(f"Starting side share {start['side_split']:.2%}, target {target / (1 + target):.2%}")
    for name, g1, g2 in zip(NAMES, coarse_grad[0], fine_grad[0]):
        print(f"{name}: {100*g1:+.2f} / {100*g2:+.2f} percentage points of side share per mm")

    try:
        result = minimize(split_objective, initial_controls, jac=split_gradient,
                          method="SLSQP", bounds=BOUNDS, constraints=[area_constraint],
                          options={"maxiter": max_steps, "ftol": 1e-6})
        print(f"Flow-balance search: {result.message}")
    except RuntimeError as error:
        print(f"Flow-balance search stopped: {error}")

    def allowed(row):
        return row["converged"] and row["volume_ok"] and row["analytic_area_mm2"] <= area_cap + 1e-6

    def rank(row):
        return outlet_violation(row, tolerance), row["pressure_drop_pa"]

    best = min(filter(allowed, records), key=rank)
    if outlet_violation(best, tolerance) == 0:
        def balance(controls):
            error = metrics_at(controls)["outlet_balance_error"]
            return np.array([tolerance - error, tolerance + error])

        minimize(lambda controls: metrics_at(controls)["pressure_drop_pa"] / 1000,
                 [best[name] for name in NAMES], jac=lambda controls: sensitivities(controls, .08)[2] / 1000,
                 method="SLSQP", bounds=BOUNDS,
                 constraints=[area_constraint, {"type": "ineq", "fun": balance}],
                 options={"maxiter": max_steps})
        best = min(filter(allowed, records), key=rank)

    controls = [best[name] for name in NAMES]
    print(f"Best Model B side/main ratio {best['side_to_main_ratio']:.4f}")
    wall_checks = []
    for wall_spacing in (.0005, .00025):
        _, run = solve_walls(controls, wall_spacing, iterations)
        wall = flow_metrics(run, target)
        wall.update(spacing_mm=wall_spacing * 1e3, converged=run["converged"],
                    mass_residual=run["global_mass_residual"], wall_leakage=float(np.max(np.abs(
                        [flux for flux, face in zip(run["boundary_flux"], run["boundary"]) if face["tag"] == "wall"]))),
                    fluid_area_mm2=float(run["areas"].sum() * 1e6))
        wall["feasible"] = (wall["converged"] and wall["mass_residual"] <= 1e-4
                            and wall["wall_leakage"] == 0 and wall["fluid_area_mm2"] <= area_cap + 1e-6
                            and outlet_violation(wall, tolerance) == 0)
        wall_checks.append(wall)
        print(f"Explicit walls at {wall['spacing_mm']:g} mm: ratio {wall['side_to_main_ratio']:.4f}, "
              f"converged {wall['converged']}, constraints met {wall['feasible']}")

    out = Path(__file__).resolve().parent / "results"
    out.mkdir(exist_ok=True)
    gradients = [dict(control=name, area_mm2_per_mm=area, side_share_per_mm_h08=g1,
                      side_share_per_mm_h04=g2, pressure_pa_per_mm_h04=p)
                 for name, area, g1, g2, p in zip(NAMES, AREA_GRADIENT, coarse_grad[0], fine_grad[0], fine_grad[2])]
    for name, rows in (("shape_gradients", gradients), ("shape_candidates", records),
                       ("shape_best", [best]), ("shape_wall_checks", wall_checks)):
        with (out / f"{name}.csv").open("w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
    print(f"Shape results saved in {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=5)
    args = parser.parse_args()
    search(args.steps)
