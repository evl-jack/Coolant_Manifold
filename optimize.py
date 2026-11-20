"""Small, constrained 2D density optimization for the T-manifold case."""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.ndimage import label

try:
    from .solve import density_geometry, solve
    from .case import CASE
except ImportError:
    from solve import density_geometry, solve
    from case import CASE

from mesh import in_region


def flow_metrics(result, side_to_main_ratio):
    """Measure pressure drop and net boundary mass flow per unit depth."""
    def port_pressure(tag):
        faces = [face for face in result["boundary"] if face["tag"] == tag]
        return np.average(result["p"][[face["cell"] for face in faces]],
                          weights=[face["length"] for face in faces])

    dp = port_pressure("inlet") - port_pressure("main_outlet")
    flow = {tag: sum(f for f, face in zip(result["boundary_flux"], result["boundary"])
                     if face["tag"] == tag) for tag in ("inlet", "main_outlet", "side_outlet")}
    inlet_flow = -flow["inlet"]
    main, side = flow["main_outlet"], flow["side_outlet"]
    return {"pressure_drop_pa": float(dp), "side_split": float(side / inlet_flow),
            "inlet_mass_flow": float(inlet_flow), "side_mass_flow": float(side),
            "main_mass_flow": float(main), "side_to_main_ratio": float(side / main) if main > 0 else np.nan,
            "outlet_balance_error": float((side - side_to_main_ratio * main) / inlet_flow)}


def outlet_violation(metrics, tolerance):
    """Measure the equality mismatch beyond tolerance and reject net backflow."""
    inlet = metrics["inlet_mass_flow"]
    return max(abs(metrics["outlet_balance_error"]) - tolerance,
               -metrics["side_mass_flow"] / inlet, -metrics["main_mass_flow"] / inlet, 0.0)


def setup(spacing=None):
    spacing = CASE["spacing_m"] if spacing is None else spacing
    centers, areas, internal, boundary, initial, x, y = density_geometry(spacing)
    shape = (len(x), len(y))
    initial = initial.reshape(shape)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    # Fix the port sections, the branch root, and the flow spines.
    branch = CASE["side_branch"]
    main_center = np.mean(CASE["main_channel"]["y"])
    main_spine = (np.abs(yy - main_center) <= spacing / 2 + 1e-12) & (initial == 1)
    branch_column = x[(x >= branch["x"][0]) & (x <= branch["x"][1])][0]
    branch_spine = (xx == branch_column) & (yy >= main_center) & (initial == 1)
    left, right = CASE["main_channel"]["x"]
    run = CASE["fixed_main_run_m"]
    end_sections = (xx < left + run) | (xx > right - run)
    root = in_region(xx, yy, CASE["branch_root"])
    side_outlet_section = yy >= branch["y"][1] - CASE["fixed_side_run_m"]
    design = ~(end_sections | root | main_spine | branch_spine | side_outlet_section)
    return centers, areas, internal, boundary, initial, x, y, design


def expand(coarse, x, y, initial, design):
    cx = np.linspace(x[design.any(axis=1)].min(), x[design.any(axis=1)].max(), coarse.shape[0])
    cy = np.linspace(y[design.any(axis=0)].min(), y[design.any(axis=0)].max(), coarse.shape[1])
    along_y = np.array([np.interp(y, cy, row) for row in coarse])
    field = np.column_stack([np.interp(x, cx, along_y[:, j]) for j in range(len(y))])
    result = initial.copy()
    result[design] = np.clip(initial[design] + field[design], 0, 1)
    return result


def project_volume(coarse, target, x, y, initial, design, areas):
    """Shift design densities uniformly until their area-weighted volume matches."""
    def shifted(offset):
        return expand(np.clip(coarse + offset, -1, 1), x, y, initial, design)
    if abs(np.average(shifted(0)[design], weights=areas[design]) - target) < 1e-10:
        return coarse
    lo, hi = -2., 2.
    for _ in range(32):
        mid = (lo + hi) / 2
        if np.average(shifted(mid)[design], weights=areas[design]) > target:
            hi = mid
        else:
            lo = mid
    return np.clip(coarse + (lo + hi) / 2, -1, 1)


def connected(density, x, y, boundary):
    active = density >= .5
    components, _ = label(active, structure=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]))
    inlet_ids = set()
    ny = active.shape[1]
    for face in boundary:
        if face["tag"] != "inlet":
            continue
        cell = face["cell"]
        i, j = divmod(cell, ny)
        if components[i, j] > 0:
            inlet_ids.add(components[i, j])
    outlets = []
    for face in boundary:
        if face["tag"] in ("main_outlet", "side_outlet"):
            i, j = divmod(face["cell"], ny)
            outlets.append(components[i, j])
    return bool(inlet_ids) and all(cid > 0 and cid in inlet_ids for cid in outlets)


def optimize(steps=5, model="B", variables=(4, 3), iterations=1000,
             sensitivity_step=.05, move_limit=.2, side_to_main_ratio=None,
             volume_fraction=None, alpha_max=1e6):
    target = CASE["side_to_main_ratio"] if side_to_main_ratio is None else side_to_main_ratio
    tolerance = CASE["outlet_balance_tolerance"]
    if steps < 1 or not np.isfinite(target) or target <= 0:
        raise ValueError("steps and side_to_main_ratio must be positive")
    _, areas, _, boundary, initial, x, y, design = setup()
    areas = areas.reshape(initial.shape)
    if volume_fraction is None:
        volume_fraction = CASE["design_volume_fraction"]
    if not 0 <= volume_fraction <= 1:
        raise ValueError("volume_fraction must be between 0 and 1")
    target_volume = volume_fraction
    # Start from the original T, then project to the requested fluid budget.
    coarse = np.zeros(tuple(variables))
    coarse = project_volume(coarse, target_volume, x, y, initial, design, areas)

    def evaluate(values, initial_state=None, required=False):
        density = expand(values, x, y, initial, design)
        if not connected(density, x, y, boundary):
            return None
        result = solve(iterations, model=model, design_density=density, alpha_max=alpha_max,
                       convergence_tolerance=5e-4, initial_state=initial_state)
        if not result["converged"]:
            if required:
                raise RuntimeError(f"Base RANS solve did not converge in {result['iterations']} iterations")
            return None
        return flow_metrics(result, target), result

    def rank(metrics):
        return outlet_violation(metrics, tolerance), metrics["pressure_drop_pa"]

    base = evaluate(coarse, required=True)
    if base is None:
        raise RuntimeError("Starting density disconnects an inlet from an outlet")

    history = []
    for step in range(steps):
        metrics, base_result = base
        violation, _ = rank(metrics)
        # Match the outlet ratio before reducing pressure drop.
        key = "outlet_balance_error" if violation else "pressure_drop_pa"
        direction = -1 if violation and metrics[key] < 0 else 1
        score = direction * metrics[key]
        grad = np.zeros_like(coarse)
        for index in np.ndindex(coarse.shape):
            plus, minus = coarse.copy(), coarse.copy()
            plus[index] = min(1., plus[index] + sensitivity_step)
            minus[index] = max(-1., minus[index] - sensitivity_step)
            plus = project_volume(plus, target_volume, x, y, initial, design, areas)
            minus = project_volume(minus, target_volume, x, y, initial, design, areas)
            fp, fm = evaluate(plus, base_result), evaluate(minus, base_result)
            if fp and fm and plus[index] != minus[index]:
                grad[index] = direction * (fp[0][key] - fm[0][key]) / (plus[index] - minus[index])
            elif fp and plus[index] != coarse[index]:
                grad[index] = (direction * fp[0][key] - score) / (plus[index] - coarse[index])
            elif fm and minus[index] != coarse[index]:
                grad[index] = (score - direction * fm[0][key]) / (coarse[index] - minus[index])
        scale = max(np.max(np.abs(grad)), 1e-30)
        accepted, step_size = False, move_limit
        for _ in range(5):
            trial = np.clip(coarse - step_size * grad / scale, -1, 1)
            trial = project_volume(trial, target_volume, x, y, initial, design, areas)
            candidate = evaluate(trial, base_result)
            if candidate and rank(candidate[0]) < rank(metrics):
                coarse, base, accepted = trial, candidate, True
                break
            step_size *= .5
        next_metrics = base[0]
        next_violation, _ = rank(next_metrics)
        history.append({"step": step, "model": model, "variables": "x".join(map(str, coarse.shape)),
                        "pressure_drop_pa": metrics["pressure_drop_pa"], "side_split": metrics["side_split"],
                        "target_side_to_main_ratio": target, "outlet_balance_tolerance": tolerance,
                        "design_volume_fraction": target_volume,
                        "gradient_norm": float(np.linalg.norm(grad)),
                        "next_pressure_drop_pa": next_metrics["pressure_drop_pa"],
                        "next_side_split": next_metrics["side_split"],
                        "side_to_main_ratio": next_metrics["side_to_main_ratio"],
                        "outlet_balance_error": next_metrics["outlet_balance_error"],
                        "outlet_balance_violation": next_violation,
                        "side_mass_flow_kg_s_m": next_metrics["side_mass_flow"],
                        "main_mass_flow_kg_s_m": next_metrics["main_mass_flow"],
                        "feasible": next_violation == 0,
                        "accepted": accepted})
        print(history[-1])

    final_density = base[1]["density"].reshape(initial.shape)
    out = Path(__file__).resolve().parent / "results"
    out.mkdir(exist_ok=True)
    with (out / "optimization.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=history[0].keys())
        writer.writeheader()
        writer.writerows(history)
    np.save(out / "optimized_density.npy", final_density)
    status = "met" if history[-1]["feasible"] else "not met"
    print(f"Outlet ratio {status}: side/main = {base[0]['side_to_main_ratio']:.4f} (target {target:.4f})")
    print(f"Density and history saved under {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--model", choices=("A", "B", "C"), default="B")
    parser.add_argument("--volume-fraction", type=float,
                        help="fluid fraction in editable cells (default: %(default)s from case.py)",
                        default=CASE["design_volume_fraction"])
    parser.add_argument("--move-limit", type=float, default=.2)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--side-to-main-ratio", type=float,
                        help="target side outlet mass flow divided by main outlet mass flow (default: 0.75)")
    args = parser.parse_args()
    optimize(steps=args.steps, model=args.model, volume_fraction=args.volume_fraction,
             move_limit=args.move_limit, iterations=args.iterations, side_to_main_ratio=args.side_to_main_ratio)
