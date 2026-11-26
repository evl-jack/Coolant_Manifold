"""Run the approved 26-control FEM pilot without a design optimizer."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import quad_vec
from scipy.interpolate import BSpline
from scipy.optimize import brentq

from .fillet import ALPHA, CONTROLS, RESPONSES, FilletMap, Problem
from .fillet_study import FEASIBILITY_TOL, GRADIENT_STEPS, verify_gradients, write
from .recovery_audit import GEOMETRY, TRANSITION_KNOTS, TRANSITION_LOCATIONS, transition


COUNT = CONTROLS+6
REFERENCE = np.full(COUNT, .5)
LOWER, UPPER = .195, .205
PRESSURE_CAP = 95.31983921402909  # Archived limit, displayed as 95.3198.
TRANSITION_ROWS = 18
PILOT_RESPONSES = (*RESPONSES, "side_to_main", "physical_fluid_area_change_mm2")


class TransitionMap(FilletMap):
    """Two separate mirrored supports with the same compact density law."""

    def __init__(self, geometry):
        super().__init__(geometry)
        self.transition_spline = BSpline(TRANSITION_KNOTS, np.eye(6), 3)

    def basis(self, y):
        s, clearance, fillet_amplitude = self.parameters(y)
        t = np.clip((y-6.5)/1.5, 0., 1.)
        branch, transition_amplitude = transition(y)
        in_fillet = (5. < y) & (y < 6.5)
        in_transition = (6.5 <= y) & (y < 8.)
        width = np.where(y < 6.5, self.stem+clearance, branch)
        amplitude = np.where(in_fillet, fillet_amplitude, np.where(in_transition, transition_amplitude, 0.))
        basis = np.concatenate((self.spline(s)*in_fillet[..., None],
                                self.transition_spline(t)*in_transition[..., None]), axis=-1)
        return width, amplitude, basis

    def offset(self, y, design=None):
        _, amplitude, basis = self.basis(np.asarray(y))
        return 1.14*amplitude if design is None else amplitude*(basis@(2*design-1))

    def prepare(self, points):
        x, y = 2*points[..., 0]-1, 2*points[..., 1]+3
        width, amplitude, basis = self.basis(y)
        distance = abs(x-10.)-width
        self.index = np.nonzero(amplitude > 0)
        self.shape = y.shape
        self.edit = (amplitude > 0) & (abs(distance) <= 1.14*amplitude)
        self.origin = .5-distance[self.index]/(.28*amplitude[self.index])
        self.design_basis = basis[self.index]/.14

    def physical_area_gradient(self):
        def integrand(y):
            _, amplitude, basis = self.basis(np.asarray(y))
            return 4*amplitude*basis
        curve_knots = 5+1.5*(1-np.cos(.5*np.pi*np.unique(self.spline.t)))
        knots = np.r_[curve_knots[1:-1], 6.5, 7., 7.5]
        return quad_vec(integrand, 5., 8., points=knots, epsabs=1e-10)[0]


class PilotProblem(Problem):
    def __init__(self, reynolds, refine=1):
        super().__init__(GEOMETRY, reynolds, refine, ALPHA, mapping=TransitionMap(GEOMETRY),
                         transition_rows=TRANSITION_ROWS)
        self.physical_gradient = self.map.physical_area_gradient()

    def evaluate(self, design, gradients=False):
        result = super().evaluate(design, gradients)
        values, jac = result if gradients else (result, None)
        f = values[1]
        values = np.r_[values, f/(1-f), self.physical_gradient@(design-REFERENCE)]
        return (values, np.vstack((jac, jac[1]/(1-f)**2, self.physical_gradient))) if gradients else values


def verify_map(problem):
    shape = problem.map
    points = problem.fem.quad_points
    y = np.linspace(5., 8., 4001)
    base_width, _, _ = shape.basis(y)
    extremes = [REFERENCE, np.zeros(COUNT), np.ones(COUNT), np.arange(COUNT) % 2]
    minimum = min(2*np.min(base_width+shape.offset(y, d)) for d in extremes)
    gamma, slope = shape.density(REFERENCE)
    rank = np.linalg.matrix_rank(slope)
    singular = np.linalg.svd(slope, compute_uv=False)
    mirrored = points.copy()
    mirrored[..., 0] = 11-mirrored[..., 0]
    reflected = TransitionMap(GEOMETRY)
    reflected.prepare(mirrored)
    mirror_error, fixed_error, min_density, max_density = 0., 0., 1., 0.
    for d in extremes:
        density, _ = shape.density(d)
        mirror_error = max(mirror_error, np.max(abs(density-reflected.density(d)[0])))
        fixed_error = max(fixed_error, np.max(abs(density[~shape.edit]-gamma[~shape.edit])))
        min_density, max_density = min(min_density, density.min()), max(max_density, density.max())
    rng = np.random.default_rng(30)
    errors = []
    for direction in (rng.normal(size=COUNT), np.r_[np.zeros(CONTROLS), np.ones(6)]):
        direction /= np.linalg.norm(direction)
        h = 1e-7
        low, high = [shape.density(REFERENCE+sign*h*direction)[0][shape.index] for sign in (-1, 1)]
        errors.append(float(np.max(abs((high-low)/(2*h)-slope@direction))))
    joins = np.array([5., 6.5, 8.])
    move = np.ones(COUNT)
    join = max(abs(shape.offset(joins, move)))
    eps = 1e-5
    probes = np.array([6.5-eps, 6.5+eps, 8.-eps, 8.+eps])
    tangent = max(abs(shape.offset(probes, move)))/eps
    # Density polynomial has zero derivative at both transition limits.
    z = np.array([0., 1.])
    ramp_endpoint_slopes = (6*z*(1-z)).tolist()
    port_mm = [float(2*np.linalg.norm(problem.fem.points[p[:, 1]]-problem.fem.points[p[:, 0]], axis=1).sum())
               for p in problem.ports]
    port_planes = [float(2*problem.fem.points[p[0, 0], dim]+origin)
                   for p, dim, origin in zip(problem.ports, (0, 0, 1), (-1., -1., 3.))]
    passed = (rank == COUNT and max(errors) < 2e-6 and mirror_error < 2e-9 and fixed_error == 0.
              and min_density >= 0. and max_density <= 1. and minimum >= 1.-1e-12 and join < 1e-12
              and tangent < 1e-4 and np.allclose(port_mm, (2., 2., 1.25), atol=1e-12, rtol=0)
              and np.allclose(port_planes, (-1., 21., 11.), atol=1e-12, rtol=0))
    return dict(passed=bool(passed), design_rank=int(rank), tangent_smallest_singular_value=float(singular[-1]),
                tangent_largest_singular_value=float(singular[0]), density_derivative_errors=errors,
                mirrored_density_error=float(mirror_error), fixed_region_error=float(fixed_error),
                density_bounds=[float(min_density), float(max_density)], ramp_endpoint_slopes=ramp_endpoint_slopes,
                minimum_width_mm=float(minimum), fixed_join_displacement_mm=float(join),
                max_endpoint_tangent_defect=float(tangent), port_widths_mm=port_mm, port_planes_mm=port_planes,
                connectivity="A positive-width single interval at each branch station joins the continuous main bore.",
                transition_greville_y_mm=TRANSITION_LOCATIONS.tolist())


def diffuse_solid(problem, design):
    gamma = problem.map.density(design)[0]
    return float(np.sum(problem.area_weights*(1-gamma)))


def compensate(problem, design, group):
    """Find the smallest affine control shift meeting both pilot area measurements."""
    baseline = diffuse_solid(problem, REFERENCE)

    def shifted(t):
        candidate = design.copy()
        scale = (1-design[group]) if t >= 0. else design[group]
        candidate[group] += t*scale
        return candidate

    def minimum_gain(t):
        candidate = shifted(t)
        return min(baseline-diffuse_solid(problem, candidate), problem.physical_gradient@(candidate-REFERENCE))

    if minimum_gain(1.) < 0. or minimum_gain(-1.) > 0.:
        raise ValueError("The chosen compensation direction cannot bracket zero net area")
    t = brentq(minimum_gain, -1., 1., xtol=1e-12)
    return shifted(min(t+1e-10, 1.))


def designs(problem, archive):
    yield "baseline", "width sweep", 0., REFERENCE.copy()
    for dial in (-1., -.75, -.5, -.25, .5, 1.):
        d = REFERENCE.copy()
        d[CONTROLS:] = .5+.5*dial
        name = "narrow" if dial < 0 else "widen"
        yield f"{name} {abs(dial):.0%}", "width sweep", dial, d
    curve = REFERENCE.copy()
    curve[CONTROLS:] += np.array([-.30, -.27, -.15, .15, .27, .30])
    yield "curvature redistribution", "curvature", None, compensate(problem, curve, slice(CONTROLS, None))
    d = REFERENCE.copy()
    d[:CONTROLS] = archive
    yield "archived fillet probe", "archive comparison", None, d.copy()
    for label, seed in (("baseline fillet", REFERENCE), ("archived fillet", d)):
        for dial in (-.5, -1.):
            expanded = seed.copy()
            expanded[CONTROLS:] = .5+.5*dial
            name = f"{label} + {abs(dial):.0%} narrowing, compensated"
            yield name, "combined", dial, compensate(problem, expanded, slice(0, CONTROLS))


def evaluate_case(problem, baseline_solid, name, mode, dial, design):
    y = np.linspace(5., 8., 4001)
    offsets = problem.map.offset(y, design)
    width = problem.map.basis(y)[0]+offsets
    diffuse_change = baseline_solid-diffuse_solid(problem, design)
    physical_change = float(problem.physical_gradient@(design-REFERENCE))
    row = dict(case=name, mode=mode, transition_dial=dial, reynolds=problem.reynolds,
               minimum_width_mm=float(2*width.min()), diffuse_fluid_area_change_mm2=diffuse_change,
               physical_fluid_area_change_mm2=physical_change,
               area_discrepancy_mm2=diffuse_change-physical_change, design=design.tolist())
    residuals = dict(diffuse_area=float(-diffuse_change/baseline_solid),
                     physical_area=float(-physical_change/baseline_solid),
                     width=float((1.-row["minimum_width_mm"])))
    try:
        values = problem.evaluate(design)
        diagnostic = problem.diagnostics(design)
        if not np.isfinite(values).all():
            raise FloatingPointError("Nonfinite hydraulic response")
        row.update(dict(zip(PILOT_RESPONSES, values.tolist())), **diagnostic)
        f, pressure = values[1:3]
        residuals.update(flow_upper=float((f-UPPER)/.005), flow_lower=float((LOWER-f)/.005),
                         pressure=float(pressure/PRESSURE_CAP-1.))
        converged = row["residual"] <= 1e-11 and row["mass_imbalance"] < 1e-9
        row.update(converged=converged, failure=None if converged else "Newton or mass tolerance missed",
                   area_feasible=max(residuals[k] for k in ("diffuse_area", "physical_area", "width")) <= FEASIBILITY_TOL,
                   pressure_feasible=residuals["pressure"] <= FEASIBILITY_TOL,
                   all_constraints_feasible=converged and max(residuals.values()) <= FEASIBILITY_TOL)
    except (RuntimeError, np.linalg.LinAlgError, FloatingPointError) as error:
        row.update(converged=False, failure=str(error), area_feasible=None, pressure_feasible=None,
                   all_constraints_feasible=False)
    row["scaled_constraint_residuals"] = residuals
    return row


def run_model(reynolds, archive, output, steps):
    start = perf_counter()
    problem = PilotProblem(reynolds)
    mapping = verify_map(problem)
    print(f"Re {reynolds:g} map: rank {mapping['design_rank']}, pass {mapping['passed']}", flush=True)
    write(output/f"map-re-{reynolds:g}.json", mapping)
    if not mapping["passed"]:
        return False
    data = json.loads((archive/"search-re-50.json").read_text())
    probe = np.array(data["probe"]["design"])
    baseline_solid = diffuse_solid(problem, REFERENCE)
    rows = []
    for name, mode, dial, d in designs(problem, probe):
        rows.append(evaluate_case(problem, baseline_solid, name, mode, dial, d))
        latest = rows[-1]
        print(f"Re {reynolds:g} {name}: side {latest.get('side_share', float('nan')):.5%}", flush=True)
    eligible = [r for r in rows if r["converged"] and r["area_feasible"] and r["pressure_feasible"]]
    saved = min(eligible, key=lambda r: r["side_share"], default=rows[0])
    if reynolds == 50. and saved["converged"]:
        problem.evaluate(np.array(saved["design"]))
        fem, (_, _, _, state) = problem.fem, problem.cache
        np.savez_compressed(output/"velocity-re-50.npz", velocity_points=fem.velocity_points,
                            velocity_dofs=fem.velocity_dofs, state=state.values, case=saved["case"])
    write(output/f"cases-re-{reynolds:g}.json", dict(reynolds=reynolds, baseline_solid_area_mm2=baseline_solid,
                                                   triangles=len(problem.fem.triangles), rows=rows,
                                                   seconds=perf_counter()-start))
    rng = np.random.default_rng(41)
    directions = (rng.normal(size=COUNT), rng.normal(size=COUNT), np.ones(COUNT),
                  np.r_[np.zeros(CONTROLS), np.ones(6)])
    checks = verify_gradients(problem, REFERENCE, directions, PILOT_RESPONSES, steps=steps)
    write(output/f"gradients-re-{reynolds:g}.json", checks)
    print(f"Re {reynolds:g}: gradients pass {checks['passed']}", flush=True)
    return checks["passed"]


def run(archive, output, serial=False, fine_steps=False):
    output.mkdir(parents=True, exist_ok=True)
    steps = (*GRADIENT_STEPS, .0003125, .00015625) if fine_steps else GRADIENT_STEPS
    write(output/"settings.json", dict(target=.20, flow_band=[LOWER, UPPER], pressure_cap=PRESSURE_CAP,
                                       original_controls=CONTROLS, transition_controls=6, minimum_width_mm=1.,
                                       alpha=ALPHA, area_rule="baseline-relative fluid area on the complete paired support",
                                       gradient_steps=steps, transition_rows=TRANSITION_ROWS))
    if serial:
        results = [run_model(re, archive, output, steps) for re in (0., 50.)]
    else:
        with ProcessPoolExecutor(max_workers=2) as pool:
            results = [pool.submit(run_model, re, archive, output, steps) for re in (0., 50.)]
            results = [job.result() for job in results]
    paths = [output/f"cases-re-{re}.json" for re in ("0", "50")]
    if all(path.exists() for path in paths):
        fields = (*[k for k in ("case", "mode", "transition_dial", "reynolds", "side_share", "side_to_main",
                                 "pressure_loss", "dissipation", "diffuse_fluid_area_change_mm2",
                                 "physical_fluid_area_change_mm2", "area_discrepancy_mm2", "minimum_width_mm",
                                 "solid_rms_speed", "converged", "mass_imbalance", "area_feasible", "pressure_feasible",
                                 "all_constraints_feasible", "failure")], "scaled_constraint_residuals")
        with (output/"pilot-cases.csv").open("w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for path in paths:
                writer.writerows(json.loads(path.read_text())["rows"])
    return all(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--serial", action="store_true")
    parser.add_argument("--fine-gradient-steps", action="store_true")
    args = parser.parse_args()
    raise SystemExit(0 if run(args.archive, args.output, args.serial, args.fine_gradient_steps) else 1)
