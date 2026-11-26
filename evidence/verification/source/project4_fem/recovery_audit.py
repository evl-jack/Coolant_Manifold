"""Replay the approved fillet model and assess an unapproved transition on paper."""

import argparse
import csv
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import quad
from scipy.sparse import coo_matrix

from .fillet import BASELINE, CONTROLS, RESPONSES, FilletMap, Problem


GEOMETRY = (1., .625, .5, 1.5)
TRANSITION_KNOTS = np.r_[np.zeros(4), 1/3, 2/3, np.ones(4)]
TRANSITION_LOCATIONS = 6.5 + 1.5*np.array([np.mean(TRANSITION_KNOTS[i+1:i+4]) for i in range(6)])


def transition(y):
    t = np.clip((y-6.5)/1.5, 0., 1.)
    baseline = .5+.125*t*t*(3-2*t)
    authority = 16*t*t*(1-t)**2*(baseline-.5)
    return baseline, authority


def outlet_metric(problem):
    """Exact boundary integral of the square of the main-outlet normal velocity."""
    fem, edges = problem.fem, problem.ports[1]
    length = np.linalg.norm(fem.points[edges[:, 1]]-fem.points[edges[:, 0]], axis=1)
    # P2 line mass matrix in endpoint, endpoint, midpoint order.
    coefficients = np.array([[4., -1., 2.], [-1., 4., 2.], [2., 2., 16.]])/30
    values = length[:, None, None]*coefficients
    rows = np.broadcast_to(edges[:, :, None], values.shape)
    cols = np.broadcast_to(edges[:, None, :], values.shape)
    matrix = coo_matrix((values.ravel(), (rows.ravel(), cols.ravel())), shape=(fem.size, fem.size)).tocsr()
    return matrix, float(length.sum())


def measured(problem, design, metric):
    measured_values = problem.evaluate(design)
    record = dict(zip(RESPONSES, measured_values.tolist()), **problem.diagnostics(design))
    gamma, slope, resistance, state = problem.cache
    w = state.values
    fem = problem.fem
    u, grad_u = fem.velocity_fields(w)
    pressure = np.einsum("qa,ea->eq", fem.p_basis, w[2*fem.n_u:][fem.triangles])
    continuity_work = -2*np.sum(fem.weights*pressure*np.trace(grad_u, axis1=-2, axis2=-1))
    advective_work = problem.reynolds*np.einsum("eq,eqi,eqj,eqij->", fem.weights, u, u, grad_u)
    operator = fem.viscous_pressure+fem.drag_matrix(resistance)
    residual = fem.residual(w, operator, problem.load, tangent=False)[0]
    reaction_work = w[problem.fixed]@residual[problem.fixed]
    k, length = metric
    q = problem.flux[1]@w
    weighted_velocity = k@w
    uniformity = length*(w@weighted_velocity)/q**2-1.
    y = np.linspace(problem.map.root, problem.map.root+problem.map.radius, 801)
    width = 2*(problem.map.stem+problem.map.parameters(y)[1]+problem.map.offset(y, design))
    record.update(minimum_width_mm=float(width.min()), outlet_nonuniformity=float(uniformity),
                  inlet_reaction_work=float(reaction_work), advective_work=float(advective_work),
                  discrete_pressure_work=float(continuity_work),
                  energy_balance_residual=float(reaction_work-record["dissipation"]-advective_work-continuity_work),
                  converged=bool(record["residual"] <= 1e-11 and record["mass_imbalance"] < 1e-9))
    return record


def uniformity_gradient(problem, design, metric):
    problem.evaluate(design)
    gamma, slope, _, state = problem.cache
    w = state.values
    k, length = metric
    q = problem.flux[1]@w
    ku = k@w
    partial = 2*length*(ku/q**2-(w@ku)*problem.flux[1]/q**3)
    derivative = np.zeros((*gamma.shape, CONTROLS))
    idx = problem.map.index
    derivative[idx] = -2*problem.alpha*(1-gamma[idx])[:, None]*slope
    return problem.fem.drag_adjoint(state, partial[None, :], derivative)[0]


def derivatives(problem, metric):
    start = perf_counter()
    _, physical = problem.evaluate(BASELINE, gradients=True)
    uniform = uniformity_gradient(problem, BASELINE, metric)
    gradient = np.vstack((physical, uniform))
    names = (*RESPONSES, "outlet_nonuniformity")
    base = measured(problem, BASELINE, metric)
    rng = np.random.default_rng(20)
    checks = []
    for raw in (np.ones(CONTROLS), rng.normal(size=CONTROLS)):
        d = raw/np.linalg.norm(raw)
        exact = float(uniform@d)
        for step in (.005, .00125):
            low, high = [measured(problem, BASELINE+sign*step*d, metric)[names[-1]] for sign in (-1, 1)]
            finite = (high-low)/(2*step)
            checks.append(dict(step=step, adjoint=exact, finite_difference=finite, absolute_error=abs(finite-exact)))
    normalized = {name: dict(gradient=row.tolist(), l2=float(np.linalg.norm(row)),
                             relative_l2=float(np.linalg.norm(row)/abs(base[name])),
                             cosine_with_dissipation=float(row@gradient[0]/(np.linalg.norm(row)*np.linalg.norm(gradient[0]))))
                  for name, row in zip(names, gradient)}
    return normalized, checks, perf_counter()-start


def theoretical_authority(mapping, design, start_share):
    """Calibrated branch-resistance proxy for a proposed shape, not a PDE calculation."""
    def fillet(y):
        clearance = mapping.parameters(np.asarray(y))[1]
        return .5+clearance+mapping.offset(np.asarray(y), design)

    fixed = quad(lambda y: (2*fillet(y))**-3, 5., 6.5)[0]+(11.-8.)/(1.25**3)
    def resistance(amplitude):
        return fixed+quad(lambda y: (2*sum(a*b for a,b in zip(transition(y), (1., amplitude))))**-3,
                          6.5, 8.)[0]

    base, ratio = resistance(0.), start_share/(1-start_share)
    cases = []
    for dial in np.linspace(-1., 1., 9):
        factor = resistance(dial)/base
        cases.append(dict(dial=float(dial), resistance_ratio=factor, side_share=ratio/(ratio+factor)))
    extreme = (fixed+1.5)/base
    cap = 2*quad(lambda y: transition(y)[1], 6.5, 8.)[0]
    result = dict(start_share=start_share, branch_resistance_ratio_required_at_upper_band=ratio/(.205/.795),
                  branch_resistance_ratio_required_at_target=ratio/.25, cases=cases,
                  extreme_min_width_resistance_ratio=extreme, extreme_min_width_side_share=ratio/(ratio+extreme),
                  area_removed_at_full_narrowing_mm2=cap,
                  max_half_width_excursion_mm=float(np.max(transition(np.linspace(6.5, 8., 1001))[1])))
    return result


def audit(archive, output):
    output.mkdir(parents=True, exist_ok=True)
    report = {"geometry_mm": GEOMETRY, "target": [.195, .205], "transition_greville_y_mm": TRANSITION_LOCATIONS.tolist()}
    all_rows = []
    for re in (0., 50.):
        archived = json.loads((archive/f"search-re-{re:g}.json").read_text())
        problem = Problem(GEOMETRY, reynolds=re)
        metric = outlet_metric(problem)
        rows = []
        start = perf_counter()
        for i, item in enumerate([*archived["samples"], archived["probe"]]):
            design = np.array(item["design"])
            row = measured(problem, design, metric)
            reference = item.get("values", item)
            delta = max(abs(row[k]-reference[k]) for k in RESPONSES)
            if delta > 5e-7 or not row["converged"]:
                raise RuntimeError(f"Re {re:g} archive replay differs by {delta:.3e} in case {i}")
            row.update(case="baseline" if i == 0 else "range probe" if i == 7 else f"archived sample {i}",
                       reynolds=re, archive_difference=float(delta))
            rows.append(row)
        replay_seconds = perf_counter()-start
        amplitude_points = []
        for value in np.linspace(-1., 1., 5):
            d = np.full(CONTROLS, .5*(1+value))
            row = measured(problem, d, metric)
            y = np.linspace(5., 6.5, 1001)
            amplitude_points.append(dict(dial=float(value), amplitude_mm=float(np.sign(value)*np.max(abs(problem.map.offset(y, d)))),
                                         **row))
        gradients, checks, derivative_seconds = derivatives(problem, metric)
        refined = Problem(GEOMETRY, reynolds=re, refine=2)
        start = perf_counter()
        fine = measured(refined, BASELINE, outlet_metric(refined))
        fine_seconds = perf_counter()-start
        candidates = {name: float(np.corrcoef([r[name] for r in rows[:7]], [r["dissipation"] for r in rows[:7]])[0, 1])
                      for name in ("pressure_loss", "outlet_nonuniformity")}
        report[f"re_{re:g}"] = dict(archived=rows, amplitude=amplitude_points, gradients=gradients, uniformity_checks=checks,
                                    correlations=candidates, fine_baseline=fine, replay_seconds=replay_seconds,
                                    derivative_seconds=derivative_seconds, fine_solve_seconds=fine_seconds,
                                    mesh_triangles=len(problem.fem.triangles), fine_triangles=len(refined.fem.triangles))
        all_rows.extend(rows)
        print(f"Re {re:g}: reproduced {len(rows)} states and evaluated outlet uniformity on two grids", flush=True)
    fitted = json.loads((archive/"explicit-probe.json").read_text())
    probe = fitted["pressure-and-area-feasibility-probe"]
    report["transition_proxy"] = theoretical_authority(FilletMap(GEOMETRY), np.array(probe["design"]),
                                                       probe["meshes"][-1]["responses"][-1]["side_share"])
    (output/"audit.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    with (output/"replayed-cases.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.archive, args.output)
