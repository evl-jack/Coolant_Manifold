"""Check and screen the Project 4 fillet formulation before solving the design problem."""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize, nnls

from .fillet import ALPHA, BASELINE, CONTROLS, RESPONSES, TOLERANCE, Problem


FEASIBILITY_TOL = 1e-6
ACTIVE_TOL = 1e-5
ITERATION_LIMIT = 35
GRADIENT_STEPS = (.04, .02, .01, .005, .0025, .00125, .000625)


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2)+"\n", encoding="utf-8")


def verify_map(problem):
    mapping = problem.map
    y = np.linspace(mapping.root-.01, mapping.root+mapping.radius+.01, 1001)
    _, clearance, _ = mapping.parameters(y)
    signs = np.array([np.zeros(CONTROLS), np.ones(CONTROLS), np.arange(CONTROLS) % 2])
    offsets = np.array([mapping.offset(y, d) for d in signs])
    fixed = (y <= mapping.root) | (y >= mapping.root+mapping.radius)
    concentrations = [mapping.density(d) for d in (BASELINE, *signs)]
    gamma, slope = concentrations[0]
    density_rank = np.linalg.matrix_rank(slope)
    direction = np.linspace(-1, 1, CONTROLS)
    h = 1e-7
    low, high = (mapping.density(BASELINE+sign*h*direction)[0] for sign in (-1, 1))
    derivative_error = float(np.max(abs((high-low)[mapping.index]/(2*h)-slope@direction)))

    mirrored = problem.fem.quad_points.copy()
    mirrored[..., 0] = 22/(2*mapping.main)-mirrored[..., 0]
    mirror = type(mapping)(problem.geometry)
    mirror.prepare(mirrored)
    reflection = float(max(np.max(abs(field-mirror.density(d)[0])) for d, (field, _) in zip((BASELINE, *signs), concentrations)))
    invariant = all(np.array_equal(field[~mapping.edit], gamma[~mapping.edit]) and
                    np.min(field) >= 0 and np.max(field) <= 1 for field, _ in concentrations)
    min_width = float(2*np.min(mapping.stem+clearance+offsets))
    # Check continuity of the tapered displacement at both fixed joins using angular coordinates.
    t = np.array([0., 1e-5, 1-1e-5, 1.])
    y_join = mapping.root+mapping.radius*(1-np.cos(.5*np.pi*t))
    join_slope = max(abs(mapping.offset(y_join, np.ones(CONTROLS))[[1, 2]])/1e-5)
    passed = (density_rank == CONTROLS and invariant and reflection < 2e-9 and derivative_error < 2e-6
              and np.max(abs(offsets[:, fixed])) < 1e-12 and min_width >= 2*mapping.stem-1e-12
              and join_slope < 1e-4)
    return dict(passed=bool(passed), density_tangent_rank=int(density_rank), reflected_density_error=reflection,
                derivative_absolute_error=derivative_error, fixed_regions_exact=bool(invariant),
                max_join_angular_slope=float(join_slope), minimum_sampled_width_mm=min_width,
                connectivity="single positive-width branch interval attached to a continuous fixed bore")


def verify_gradients(problem, baseline=BASELINE, directions=None, response_names=RESPONSES, steps=GRADIENT_STEPS):
    start = perf_counter()
    base, jacobian = problem.evaluate(baseline, gradients=True)
    base_state = problem.diagnostics(baseline)
    checks = []
    if directions is None:
        rng = np.random.default_rng(41)
        directions = (rng.normal(size=len(baseline)), rng.normal(size=len(baseline)), np.ones(len(baseline)))
    for direction in directions:
        d = direction/np.linalg.norm(direction)
        adjoint = jacobian@d
        records = []
        prior_remainder = None
        perturbations_pass = True
        for h in steps:
            values = []
            for sign in (-1, 1):
                coordinates = baseline+sign*h*d
                values.append(problem.evaluate(coordinates))
                diagnostics = problem.diagnostics(coordinates)
                perturbations_pass &= diagnostics["mass_imbalance"] < 1e-9 and diagnostics["residual"] <= 1e-11
            low, high = values
            finite = (high-low)/(2*h)
            absolute = abs(finite-adjoint)
            remainder = abs(high-base-h*adjoint)
            order = None if prior_remainder is None else np.log2(np.maximum(prior_remainder, 1e-20)/np.maximum(remainder, 1e-20))
            records.append(dict(step=h, finite_difference=finite.tolist(), absolute_error=absolute.tolist(),
                                relative_error=(absolute/np.maximum(abs(adjoint), 1e-12)).tolist(),
                                taylor_remainder=remainder.tolist(), taylor_order=None if order is None else order.tolist()))
            prior_remainder = remainder
        orders = np.array([r["taylor_order"] for r in records[1:]])
        order_ok = np.any((orders > 1.7) & (orders < 2.3), axis=0) | (prior_remainder < 1e-12)
        passed = np.all(absolute <= np.maximum(1e-7, .002*abs(adjoint))) and np.all(order_ok) and perturbations_pass
        checks.append(dict(direction=d.tolist(), adjoint=adjoint.tolist(), steps=records,
                           taylor_passed=order_ok.tolist(), perturbations_converged=bool(perturbations_pass), passed=bool(passed)))
    return dict(passed=all(ch["passed"] for ch in checks) and base_state["mass_imbalance"] < 1e-9,
                responses=response_names, baseline=base.tolist(), baseline_state=base_state, checks=checks,
                wall_seconds=perf_counter()-start)


class Constraints:
    def __init__(self, target, baseline):
        self.target = target
        self.objective_scale = baseline[0]
        self.pressure_cap = 1.05*baseline[2]
        self.area_cap = baseline[3]
        if min(self.objective_scale, self.pressure_cap, self.area_cap) <= 0:
            raise ValueError("Baseline normalization requires positive pressure loss, dissipation and solid area")

    def values(self, r):
        return np.array([(r[1]-self.target)/TOLERANCE-1., (self.target-r[1])/TOLERANCE-1.,
                         r[2]/self.pressure_cap-1., r[3]/self.area_cap-1.])

    def jacobian(self, j):
        return np.stack((j[1]/TOLERANCE, -j[1]/TOLERANCE, j[2]/self.pressure_cap, j[3]/self.area_cap))


def kkt(x, objective_gradient, residual, jacobian):
    values = np.r_[residual, -x, x-1]
    gradients = np.vstack((jacobian, -np.eye(len(x)), np.eye(len(x))))
    names = ["side upper", "side lower", "pressure", "area"] + [f"lower {i+1}" for i in range(len(x))] + [f"upper {i+1}" for i in range(len(x))]
    active = values >= -ACTIVE_TOL
    multipliers = np.zeros(len(values))
    if active.any():
        multipliers[active] = nnls(gradients[active].T, -objective_gradient)[0]
    return dict(primal_violation=float(max(0, values.max())), dual_violation=float(max(0, -multipliers.min())),
                stationarity=float(np.linalg.norm(objective_gradient+multipliers@gradients, ord=np.inf)),
                complementarity=float(np.max(abs(multipliers*values))),
                active_multipliers={name: float(value) for name, value, keep in zip(names, multipliers, active) if keep})


def search(problem, limits, output, optimize_design=False):
    rows = []
    last = None
    for design in (BASELINE, np.zeros(CONTROLS), np.ones(CONTROLS),
                   np.linspace(0, 1, CONTROLS), np.linspace(1, 0, CONTROLS),
                   .5+.45*np.sin(2*np.pi*np.linspace(0, 1, CONTROLS)),
                   .5-.45*np.sin(2*np.pi*np.linspace(0, 1, CONTROLS))):
        response = problem.evaluate(design)
        rows.append(dict(design=design.tolist(), **dict(zip(RESPONSES, response.tolist())),
                         scaled_constraints=limits.values(response).tolist(), **problem.diagnostics(design)))
    print(f"Re {problem.reynolds:g} screen side share range "
          f"{min(r['side_share'] for r in rows):.4%} to {max(r['side_share'] for r in rows):.4%}", flush=True)
    # Separate constrained range probes from the target-band optimization.
    side = rows[0]["side_share"]
    orientation = 1. if side >= limits.target else -1.
    start = min((r for r in rows if max(r["scaled_constraints"][2:]) <= FEASIBILITY_TOL),
                key=lambda row: orientation*row["side_share"], default=rows[0])
    history = []

    def trial(x):
        nonlocal last
        if last is None or not np.array_equal(x, last[0]):
            last = (x.copy(), *problem.evaluate(x, gradients=True))
        return last[1:]

    def constraint(x, first):
        return -limits.values(trial(x)[0])[first:]

    def constraint_jac(x, first):
        return -limits.jacobian(trial(x)[1])[first:]

    def log(x):
        r = trial(x)[0]
        history.append(dict(iteration=len(history)+1, design=x.tolist(), values=r.tolist(), constraints=limits.values(r).tolist()))

    probe = minimize(lambda x: orientation*trial(x)[0][1], np.array(start["design"]),
                     jac=lambda x: orientation*trial(x)[1][1], bounds=[(0, 1)]*CONTROLS,
                     constraints=[dict(type="ineq", fun=lambda x: constraint(x, 2), jac=lambda x: constraint_jac(x, 2))],
                     method="SLSQP", callback=log, options=dict(maxiter=ITERATION_LIMIT, ftol=1e-9))
    result, _ = trial(probe.x)
    base_jac = problem.evaluate(BASELINE, True)[1]
    cosine = float(base_jac[0]@base_jac[2]/(np.linalg.norm(base_jac[0])*np.linalg.norm(base_jac[2])))
    correlated = np.corrcoef([r["dissipation"] for r in rows], [r["pressure_loss"] for r in rows])[0, 1]
    feasible = [r for r in rows if max(r["scaled_constraints"]) <= FEASIBILITY_TOL]
    if max(limits.values(result)) <= FEASIBILITY_TOL:
        feasible.append(dict(design=probe.x.tolist(), dissipation=result[0]))
    if not feasible and max(limits.values(result)[2:]) <= FEASIBILITY_TOL:
        # If the range probe crossed the target, search the segment back to the baseline.
        for t in np.linspace(.1, .9, 9):
            x = t*probe.x+(1-t)*BASELINE
            r = trial(x)[0]
            if max(limits.values(r)) <= FEASIBILITY_TOL:
                feasible.append(dict(design=x.tolist(), dissipation=r[0]))
    print(f"Re {problem.reynolds:g} constrained probe side share {result[1]:.4%}; {probe.message}", flush=True)
    report = dict(reynolds=problem.reynolds, samples=rows, pressure_dissipation_sample_correlation=float(correlated),
                  pressure_dissipation_gradient_cosine=cosine, observed_feasible=bool(feasible),
                  probe=dict(values=dict(zip(RESPONSES, result.tolist())), design=probe.x.tolist(),
                             constraints=limits.values(result).tolist(), termination=str(probe.message), history=history.copy()))
    if feasible and optimize_design:
        history.clear()
        seed = min(feasible, key=lambda r: r["dissipation"])
        optimum = minimize(lambda x: trial(x)[0][0]/limits.objective_scale, np.array(seed["design"]),
                           jac=lambda x: trial(x)[1][0]/limits.objective_scale, bounds=[(0, 1)]*CONTROLS,
                           constraints=[dict(type="ineq", fun=lambda x: constraint(x, 0), jac=lambda x: constraint_jac(x, 0))],
                           method="SLSQP", callback=log, options=dict(maxiter=ITERATION_LIMIT, ftol=1e-9))
        r, jac = trial(optimum.x)
        stationarity = kkt(optimum.x, jac[0]/limits.objective_scale, limits.values(r), limits.jacobian(jac))
        report["optimum"] = dict(design=optimum.x.tolist(), values=dict(zip(RESPONSES, r.tolist())),
                                 scaled_constraints=limits.values(r).tolist(), kkt=stationarity,
                                 termination=str(optimum.message), converged=bool(optimum.success), history=history)
    write(output/f"search-re-{problem.reynolds:g}.json", report)
    return report


def run(output, geometry, side_main_target):
    output.mkdir(parents=True, exist_ok=True)
    ns, stokes = Problem(geometry), Problem(geometry, reynolds=0.)
    mapping = verify_map(ns)
    write(output/"map.json", mapping)
    print(json.dumps(mapping, indent=2), flush=True)
    if not mapping["passed"]:
        return False
    leakage = []
    for alpha in (5000., ALPHA, 80000.):
        model = ns if alpha == ALPHA else Problem(geometry, alpha=alpha)
        leakage.append(dict(alpha=alpha, **dict(zip(RESPONSES, model.evaluate(BASELINE).tolist())),
                            **model.diagnostics(BASELINE)))
    write(output/"penalty.json", leakage)
    ns_gradients, stokes_gradients = verify_gradients(ns), verify_gradients(stokes)
    write(output/"gradients-re-50.json", ns_gradients)
    write(output/"gradients-re-0.json", stokes_gradients)
    print(f"Adjoint checks: Re 50 {ns_gradients['passed']}; Stokes {stokes_gradients['passed']}", flush=True)
    if not (ns_gradients["passed"] and stokes_gradients["passed"]):
        return False
    constraints = Constraints(side_main_target/(1+side_main_target), ns.evaluate(BASELINE))
    write(output/"limits.json", dict(**vars(constraints), side_main_target=side_main_target, side_tolerance=TOLERANCE,
                                     alpha_max=ALPHA, response_names=RESPONSES, active_tolerance=ACTIVE_TOL,
                                     scaled_feasibility_tolerance=FEASIBILITY_TOL))
    for problem in (ns, stokes):
        search(problem, constraints, output, optimize_design=True)
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("project4_fem/output/fillet-study"))
    parser.add_argument("--side-main-target", type=float, default=.25)
    args = parser.parse_args()
    if not np.isfinite(args.side_main_target) or args.side_main_target <= 0:
        parser.error("The side/main target must be positive and finite")
    raise SystemExit(0 if run(args.output, (1., .625, .5, 1.5), args.side_main_target) else 1)
