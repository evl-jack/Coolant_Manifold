"""Matched Brinkman and boundary-fitted checks for two archived pilot shapes."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from .fem import TaylorHood
from .fillet import ALPHA, pressure_weights
from .fillet_study import FEASIBILITY_TOL, write
from .junction import boundary_data, junction_mesh
from .transition_pilot import GEOMETRY, LOWER, PRESSURE_CAP, REFERENCE, TRANSITION_ROWS, UPPER, TransitionMap


LEVELS = (1, 2, 3)
METHODS = ("Brinkman", "Fitted")
MODELS = (0., 50.)


def archived_shapes(pilot):
    cases = json.loads((pilot/"cases-re-50.json").read_text())["rows"]
    baseline = next(row for row in cases if row["case"] == "baseline")
    candidate = min((row for row in cases if row["converged"] and row["area_feasible"]
                     and row["pressure_feasible"]), key=lambda row: row["side_share"])
    return {"Baseline": baseline, "Candidate": candidate}


def measure(method, level, name, case):
    start = perf_counter()
    design = np.array(case["design"])
    mapping = TransitionMap(GEOMETRY)
    immersed = method == "Brinkman"
    offset = mapping.offset if immersed else lambda y: mapping.offset(y, design)
    fem = TaylorHood(*junction_mesh(GEOMETRY, level, offset, TRANSITION_ROWS))
    fixed, values, ports, flux = boundary_data(fem)
    reference_pressure, outlet_pressure = pressure_weights(fem, ports)
    scale = 2*GEOMETRY[0]
    dimensions = (0, 0, 1)
    origins = (-1., -1., 3.)
    widths = [scale*np.linalg.norm(fem.points[p[:, 1]]-fem.points[p[:, 0]], axis=1).sum() for p in ports]
    planes = [scale*fem.points[p[:, :2], k]+origin for p, k, origin in zip(ports, dimensions, origins)]
    port_ok = (np.allclose(widths, (2., 2., 1.25), rtol=0., atol=1e-12)
               and all(np.allclose(actual, expected, rtol=0., atol=1e-12)
                       for actual, expected in zip(planes, (-1., 21., 11.))))
    if immersed:
        mapping.prepare(fem.quad_points)
        density, _ = mapping.density(design)
        resistance = ALPHA*(1-density)**2
        area = -scale**2*np.sum(fem.weights*mapping.edit*(1-density))
        solid = density < .05
    else:
        resistance, solid = 0., None
        area = scale**2*fem.weights.sum()

    y = np.linspace(5., 8., 4001)
    width = float(2*np.min(mapping.basis(y)[0]+mapping.offset(y, design)))
    physical_area = float(mapping.physical_area_gradient()@(design-REFERENCE))
    load = fem.load()
    initial = None
    rows = []
    for re in MODELS:
        row = dict(method=method, level=level, shape=name, reynolds=re, triangles=len(fem.triangles),
                   dofs=fem.size, discrete_area_coordinate_mm2=float(area), physical_fluid_area_change_mm2=physical_area,
                   minimum_width_mm=width, fixed_ports_passed=bool(port_ok))
        try:
            fem.inertia = re
            state = fem.solve(fixed, values, load, resistance=resistance, initial=initial)
            initial = state.values
            inlet, main, side = flux@initial
            velocity = initial[:2*fem.n_u]
            viscous = velocity@(fem.viscous_pressure[:2*fem.n_u, :2*fem.n_u]@velocity)
            speed2 = np.sum(fem.velocity_fields(initial)[0]**2, axis=-1)
            drag = np.sum(fem.weights*resistance*speed2)
            leakage = None if solid is None else np.sqrt(np.average(speed2[solid], weights=fem.weights[solid]))
            share = side/inlet
            pressure = reference_pressure@initial-initial@(outlet_pressure@initial)/(main+side)
            imbalance = abs(inlet-main-side)/inlet
            residual = state.residual_history[-1]
            passed = port_ok and abs(inlet-1.) <= 1e-12 and residual <= 1e-11 and imbalance < 1e-9
            row.update(side_share=float(share), side_to_main=float(side/main), pressure_loss=float(pressure),
                       viscous_dissipation=float(viscous), brinkman_dissipation=float(drag),
                       dissipation=float(viscous+drag), inlet_flow=float(inlet), mass_imbalance=float(imbalance),
                       residual=residual, newton_updates=len(state.residual_history)-1,
                       solid_rms_speed=None if leakage is None else float(leakage), converged=bool(passed),
                       failure=None if passed else "Residual, mass, or port check failed",
                       flow_passed=bool(max((share-UPPER)/.005, (LOWER-share)/.005) <= FEASIBILITY_TOL),
                       pressure_passed=bool(pressure/PRESSURE_CAP-1 <= FEASIBILITY_TOL))
        except (RuntimeError, FloatingPointError, np.linalg.LinAlgError) as error:
            row.update(converged=False, failure=str(error))
        rows.append(row)
        print(f"{method} L{level} {name} Re {re:g} pass={row['converged']} "
              f"f={row.get('side_share', float('nan')):.6%}", flush=True)
    return rows, perf_counter()-start


def analyze(rows, area_scale):
    index = {(r["method"], r["level"], r["shape"], r["reynolds"]): r for r in rows}
    comparisons = {}
    for re in MODELS:
        grid, wall, gains = {}, {}, {}
        for method in METHODS:
            pairs = [tuple(index[method, level, name, re] for name in ("Baseline", "Candidate")) for level in LEVELS]
            for baseline, candidate in pairs:
                base_area = baseline["discrete_area_coordinate_mm2"]
                for row in (baseline, candidate):
                    area = row["discrete_area_coordinate_mm2"]-base_area
                    exact = row["physical_fluid_area_change_mm2"]
                    # Use the archived diffuse strip normalization for both wall treatments.
                    area_tol = FEASIBILITY_TOL*area_scale
                    row.update(discrete_fluid_area_change_mm2=area, area_discrepancy_mm2=area-exact,
                               discrete_area_passed=area >= -area_tol, physical_area_passed=exact >= -area_tol,
                               width_passed=row["minimum_width_mm"] >= 1.-1e-12,
                               area_agreement_passed=abs(area-exact) <= .002)
            gains[method] = [100*(base["side_share"]-candidate["side_share"])
                             if base["converged"] and candidate["converged"] else None for base, candidate in pairs]
            grid[method] = {fine["shape"]: difference(mid, fine) for mid, fine in zip(*pairs[-2:])}
        for name in ("Baseline", "Candidate"):
            wall[name] = difference(index["Brinkman", 3, name, re], index["Fitted", 3, name, re])
        resolved = dict(grid_changes=grid, wall_changes=wall, gains_pp=gains, gain_resolved=False)
        if all(gain is not None for method in METHODS for gain in gains[method][-2:]):
            grid_gain = max(abs(g[2]-g[1]) for g in gains.values())
            wall_gain = abs(gains["Brinkman"][2]-gains["Fitted"][2])
            grid_ok = all(change["side_passed"] for group in grid.values() for change in group.values())
            resolved.update(gain_grid_pp=grid_gain, gain_wall_pp=wall_gain, gain_envelope_pp=grid_gain+wall_gain,
                            gain_resolved=grid_ok and min(g[2] for g in gains.values()) > grid_gain+wall_gain)
        comparisons[f"{re:g}"] = resolved
    return comparisons


def difference(before, after):
    if not before["converged"] or not after["converged"]:
        return {"converged": False, "side_passed": False}
    share = 100*(after["side_share"]-before["side_share"])
    pressure = 100*(after["pressure_loss"]-before["pressure_loss"])/abs(after["pressure_loss"])
    energy = 100*(after["dissipation"]-before["dissipation"])/abs(after["dissipation"])
    leakage = None
    if after["solid_rms_speed"] is not None:
        leakage = 100*(after["solid_rms_speed"]-before["solid_rms_speed"])/after["solid_rms_speed"]
    return dict(converged=True, side_pp=share, side_passed=abs(share) <= .10,
                pressure_relative_pct=pressure, pressure_passed=abs(pressure) <= .5,
                dissipation_relative_pct=energy, dissipation_passed=abs(energy) <= .5,
                leakage_relative_pct=leakage, leakage_passed=None if leakage is None else abs(leakage) <= 10.)


def run(pilot, output, workers):
    output.mkdir(parents=True, exist_ok=True)
    shapes = archived_shapes(pilot)
    arguments = [(m, level, name, case) for level in LEVELS for m in METHODS for name, case in shapes.items()]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = [pool.submit(measure, *args) for args in arguments]
        groups = [job.result() for job in results]
    rows = [row for group, _ in groups for row in group]
    area_scale = json.loads((pilot/"cases-re-50.json").read_text())["baseline_solid_area_mm2"]
    comparisons = analyze(rows, area_scale)
    write(output/"matched-verification.json", dict(shapes={n: r["design"] for n, r in shapes.items()},
          archived_candidate=shapes["Candidate"]["case"], rows=rows, comparisons=comparisons,
          worker_seconds=[secs for _, secs in groups]))
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (output/"matched-verification.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    run(args.pilot, args.output, args.workers)
