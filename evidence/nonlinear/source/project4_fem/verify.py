"""Analytical and numerical checks for the Project 4 finite-element solver."""

import argparse
import json
from math import factorial
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import quad

from .fem import LOCAL_EDGES, TaylorHood, quadratic_basis, rectangle, triangle_rule


def error(actual, exact, weights):
    if actual.ndim == weights.ndim + 1:
        weights = weights[..., None]
    return float(np.sqrt(np.sum(weights * (actual - exact)**2) / np.sum(weights * exact**2)))


def fields_error(fem, state, velocity, pressure):
    u, _ = fem.velocity_fields(state.values)
    p = np.einsum("qa,ea->eq", fem.p_basis, state.values[2*fem.n_u:][fem.triangles])
    reference_p = pressure(fem.quad_points)
    return dict(velocity_l2=error(u, velocity(fem.quad_points), fem.weights),
                pressure_l2=error(p, reference_p, fem.weights))


def element_check():
    unit = np.eye(3)
    nodes = np.vstack((unit, [(unit[i] + unit[j])/2 for i, j in LOCAL_EDGES]))
    basis, derivative = quadratic_basis(nodes)
    lambdas, weights = triangle_rule()
    x, y = lambdas[:, 1:].T
    moments = [abs(weights @ (x**i * y**j) - factorial(i)*factorial(j)/factorial(i+j+2))
               for i in range(8) for j in range(8-i)]
    return dict(interpolation_defect=float(np.max(abs(basis-np.eye(6)))),
                partition_defect=float(np.max(abs(basis.sum(axis=1)-1))),
                gradient_sum_defect=float(np.max(abs(derivative.sum(axis=1)))),
                moment_error=float(max(moments)))


def layered_profile(alpha, gradient):
    """Exact symmetric channel with a clear core and two Brinkman layers."""
    half_core, half_solid = .25, .25
    k = np.sqrt(alpha)
    offset = gradient / alpha
    B = -gradient * half_core / k
    A = -(offset + B*np.sinh(k*half_solid))/np.cosh(k*half_solid)

    def stream(y):
        distance = abs(y-.5)
        z = k*np.maximum(distance-half_core, 0.)
        core = offset+A+.5*gradient*(half_core**2-distance**2)
        solid = offset+A*np.cosh(z)+B*np.sinh(z)
        return np.where(distance <= half_core, core, solid)

    def shear(y):
        distance = abs(y-.5)
        z = k*np.maximum(distance-half_core, 0.)
        slope = np.where(distance <= half_core, -gradient*distance, k*(A*np.sinh(z)+B*np.cosh(z)))
        return np.sign(y-.5)*slope

    return stream, shear


def channel(nx, ny, *, alpha=0., angle=0., inertia=0., layered=False):
    length, gradient = 2., 12.
    points, triangles = rectangle(nx, ny, length)
    c, s = np.cos(angle), np.sin(angle)
    rotation = np.array([[c, -s], [s, c]])
    fem = TaylorHood(points @ rotation.T, triangles, inertia=inertia)

    def stream(y):
        if not alpha:
            return .5 * gradient * y * (1-y)
        return gradient / alpha * (1-np.cosh(np.sqrt(alpha)*(y-.5))/np.cosh(np.sqrt(alpha)*.5))

    def shear(y):
        if not alpha:
            return gradient * (.5-y)
        return -gradient / np.sqrt(alpha) * np.sinh(np.sqrt(alpha)*(y-.5))/np.cosh(np.sqrt(alpha)*.5)

    if layered:
        stream, shear = layered_profile(alpha, gradient)

    def velocity(position):
        u = stream((position @ rotation)[..., 1])
        return np.stack((u, np.zeros_like(u)), axis=-1) @ rotation.T

    def pressure(position):
        return gradient * (length-(position @ rotation)[..., 0])

    def traction(position):
        tau = shear((position @ rotation)[..., 1])
        return np.stack((np.zeros_like(tau), tau), axis=-1) @ rotation.T

    x, y = (fem.velocity_points @ rotation).T
    fixed, values = fem.velocity_boundary(np.flatnonzero(np.isclose(x, 0) | np.isclose(y, 0) | np.isclose(y, 1)), velocity)
    edges = fem.boundary_edges
    outlet = edges[np.isclose(x[edges[:, 2]], length)]
    inlet = edges[np.isclose(x[edges[:, 2]], 0)]
    rhs = fem.load(traction=(outlet, traction))
    resistance = alpha * (abs((fem.quad_points@rotation)[..., 1]-.5) > .25) if layered else alpha
    state = fem.solve(fixed, values, rhs, resistance=resistance)
    incoming = fem.flux_functional(inlet, rotation[:, 0]) @ state.values
    outgoing = fem.flux_functional(outlet, rotation[:, 0]) @ state.values
    exact_flow = gradient/12 if not alpha else gradient/alpha*(1-2*np.tanh(np.sqrt(alpha)/2)/np.sqrt(alpha))
    result = dict(**fields_error(fem, state, velocity, pressure), exact_flow=exact_flow,
                  computed_flow=float(outgoing), relative_flux_error=float(abs(outgoing-exact_flow)/exact_flow),
                  relative_mass_imbalance=float(abs(outgoing-incoming)/incoming))
    if layered:
        exact_flow = quad(stream, 0., 1., points=(.25, .75))[0]
        solid_outlet = outlet[abs(y[outlet[:, 2]]-.5) > .25]
        solid_flow = fem.flux_functional(solid_outlet, rotation[:, 0]) @ state.values
        exact_solid_flow = 2*quad(stream, 0., .25)[0]
        result.update(alpha=alpha, exact_flow=exact_flow, relative_flux_error=float(abs(outgoing-exact_flow)/exact_flow),
                      solid_flow=float(solid_flow), exact_solid_flow=exact_solid_flow,
                      solid_flow_fraction=float(solid_flow/outgoing), clear_channel_exact_flow=gradient*.5**3/12)
    return fem, state, fixed, values, rhs, outlet, inlet, result


def manufactured(n, inertia=0.):
    fem = TaylorHood(*rectangle(n, n), inertia=inertia)

    def components(position):
        x, y = position[..., 0], position[..., 1]
        return [(t*t*(1-t)**2, 2*t-6*t*t+4*t**3, 2-12*t+12*t*t, -12+24*t) for t in (x, y)]

    def velocity(position):
        (f, df, _, _), (g, dg, _, _) = components(position)
        return 8*np.stack((f*dg, -df*g), axis=-1)

    def pressure(position):
        return position[..., 0] + .5*position[..., 1]

    def force(position):
        (f, df, ddf, dddf), (g, dg, ddg, dddg) = components(position)
        laplacian = 8*np.stack((ddf*dg+f*dddg, -dddf*g-df*ddg), axis=-1)
        ux, uy = np.moveaxis(velocity(position), -1, 0)
        advection = 8*np.stack((ux*df*dg+uy*f*ddg, -ux*ddf*g-uy*df*dg), axis=-1)
        return -laplacian + np.array([1., .5]) + inertia*advection

    fixed, values = fem.velocity_boundary(np.unique(fem.boundary_edges), velocity, (0, 0.))
    rhs = fem.load(force=force)
    state = fem.solve(fixed, values, rhs)
    result = dict(n=n, **fields_error(fem, state, velocity, pressure), residuals=state.residual_history)
    if inertia:
        direction = np.zeros(fem.size)
        direction[state.free] = np.random.default_rng(4).normal(size=len(state.free))
        direction /= np.linalg.norm(direction)
        step = 1e-5
        low, high = [fem.residual(state.values+sign*step*direction, fem.viscous_pressure, rhs, False)[0]
                     for sign in (-1, 1)]
        finite = (high[state.free]-low[state.free])/(2*step)
        exact = (state.tangent@direction)[state.free]
        result["tangent_relative_error"] = float(np.linalg.norm(finite-exact)/np.linalg.norm(exact))
    return result


def adjoint_check():
    fem, _, fixed, values, rhs, outlet, inlet, _ = channel(8, 6, inertia=15.)
    centers = np.array(np.meshgrid(np.linspace(.16, 1.84, 6), np.linspace(.15, .85, 4))).reshape(2, -1).T
    distances = fem.quad_points[..., None, :] - centers
    basis = np.exp(-np.sum(distances**2, axis=-1)/.045)
    basis /= np.maximum(1., basis.sum(axis=-1, keepdims=True))
    design = np.linspace(.33, .67, len(centers))

    def analysis(d):
        gamma = .28 + basis@d
        return fem.solve(fixed, values, rhs, resistance=25*(1-gamma)**2), gamma

    def pressure_integral(edges):
        row = np.zeros(fem.size)
        length = np.linalg.norm(fem.points[edges[:, 1]]-fem.points[edges[:, 0]], axis=1)
        np.add.at(row, 2*fem.n_u+edges[:, :2], length[:, None]*np.ones((1, 2))/2)
        return row

    upper = outlet[fem.velocity_points[outlet[:, 2], 1] > .5]
    metrics = np.stack((pressure_integral(inlet)-pressure_integral(outlet),
                        fem.flux_functional(upper, (1, 0))))
    state, gamma = analysis(design)
    gradient = fem.drag_adjoint(state, metrics, -50*(1-gamma[..., None])*basis)
    results = []
    rng = np.random.default_rng(7)
    for direction in (np.ones(len(design)), rng.normal(size=len(design))):
        direction /= np.linalg.norm(direction)
        exact = gradient@direction
        response = metrics@state.values
        previous_step, previous_remainder = None, None
        for step in (4e-2, 2e-2, 1e-2, 5e-3, 2.5e-3):
            low, high = [metrics@analysis(design+sign*step*direction)[0].values for sign in (-1, 1)]
            finite = (high-low)/(2*step)
            remainder = abs(high-response-step*exact)
            record = dict(step=step, adjoint=exact.tolist(), finite_difference=finite.tolist(),
                          relative_errors=(abs(finite-exact)/np.maximum(abs(exact), 1e-10)).tolist(),
                          taylor_remainders=remainder.tolist())
            if previous_step is not None:
                record["taylor_orders"] = (np.log(previous_remainder/remainder)/np.log(previous_step/step)).tolist()
            results.append(record)
            previous_step, previous_remainder = step, remainder
    return dict(design_variables=len(design), metrics=["inlet-minus-outlet pressure", "upper-half outlet flow"], data=results)


def verify(output):
    start = perf_counter()
    element = element_check()
    straight = channel(8, 6)[-1]
    rotated = channel(8, 6, angle=.61)[-1]
    refinement = [manufactured(n) for n in (4, 8, 16)]
    orders = {key: list(np.log2([a[key]/b[key] for a, b in zip(refinement, refinement[1:])]))
              for key in ("velocity_l2", "pressure_l2")}
    nonlinear = manufactured(12, inertia=35.)
    brinkman = channel(16, 16, alpha=25.)[-1]
    leakage = [channel(8, 64, alpha=alpha, layered=True)[-1] for alpha in (100., 400., 1600.)]
    adjoint = adjoint_check()
    limits = dict(elements=1e-12, poiseuille=1e-8, rotation=1e-8, mass=1e-10,
                  velocity_order=2.6, pressure_order=1.7, nonlinear_velocity=0.01,
                  nonlinear_tangent=1e-7, brinkman_flow=5e-4, layered_flow=5e-4, layered_solid_flow=1e-3,
                  adjoint=2e-4, taylor_order=1.8)
    smallest_steps = adjoint["data"][4::5]
    passed = dict(element=max(element.values()) < limits["elements"],
                  poiseuille=max(straight[key] for key in ("velocity_l2", "pressure_l2", "relative_flux_error")) < limits["poiseuille"],
                  rotation=max(rotated[key] for key in ("velocity_l2", "pressure_l2", "relative_flux_error")) < limits["rotation"],
                  mass=max(case["relative_mass_imbalance"] for case in (straight, rotated, brinkman)) < limits["mass"],
                  refinement=orders["velocity_l2"][-1] > limits["velocity_order"] and orders["pressure_l2"][-1] > limits["pressure_order"],
                  nonlinear=nonlinear["velocity_l2"] < limits["nonlinear_velocity"] and nonlinear["tangent_relative_error"] < limits["nonlinear_tangent"],
                  brinkman=brinkman["relative_flux_error"] < limits["brinkman_flow"],
                  penalty_leakage=all(row["relative_flux_error"] < limits["layered_flow"] and
                                      abs(row["solid_flow"]/row["exact_solid_flow"]-1) < limits["layered_solid_flow"]
                                      for row in leakage) and
                  all(a["solid_flow_fraction"] > b["solid_flow_fraction"] for a, b in zip(leakage, leakage[1:])),
                  adjoint=all(max(row["relative_errors"]) < limits["adjoint"] and min(row["taylor_orders"]) > limits["taylor_order"]
                              for row in smallest_steps))
    passed = {name: bool(result) for name, result in passed.items()}
    report = dict(passed=passed, thresholds=limits, element=element, poiseuille=straight, rotated=rotated,
                  manufactured=refinement, orders=orders, nonlinear=nonlinear, brinkman=brinkman, layered=leakage,
                  adjoint=adjoint, wall_seconds=perf_counter()-start)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(dict(passed=passed, orders=orders, wall_seconds=report["wall_seconds"]), indent=2), flush=True)
    return all(passed.values())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("project4_fem/output/verification.json"))
    args = parser.parse_args()
    raise SystemExit(0 if verify(args.output) else 1)
