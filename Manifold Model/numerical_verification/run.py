"""Run the six blind numerical verification cases from the supplied inputs."""

import argparse
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy

from fem import (Assembly, Problem, RectangleGrid, RectangleMesh, basis, direct,
                 line_rule, newton, vector_constant)


MU, DEPTH, RHO = .001, .005, 1000.
S, V, P_STAR = .10, .001, .001

# Set the study sequences and acceptance criteria before running either check.
NEWTON_SEEDS = (317, 911)
NEWTON_AMPLITUDES = (.5, .25, .125, .0625, .03125, .015625)
NEWTON_ERROR_FLOOR = 1e-10
QUADRATIC_BAND = (1.8, 2.2)
TANGENT_ADJOINT_TOL = 1e-8
CENTERED_DIFFERENCE_TOL = 1e-5
CENTERED_DIFFERENCE_EPSILON = 1e-4
TAYLOR_DEPTH_EPSILONS = (1e-2, 1e-3, 1e-4)


# Polynomial factors for the smooth Stokes manufactured field and its forcing.


def f(t):
    return t**2*(1-t)**2


def df(t):
    return 2*t - 6*t**2 + 4*t**3


def ddf(t):
    return 2 - 12*t + 12*t**2


def dddf(t):
    return -12 + 24*t


def linear_force(x, y):
    """Assemble the supplied Stokes forcing in normalized coordinates."""
    x, y = x/S, y/S
    a = np.array([f(x)*df(y), -df(x)*f(y)])
    b = np.array([ddf(x)*df(y) + f(x)*dddf(y), -dddf(x)*f(y) - df(x)*ddf(y)])
    return .01*np.array([y-.5, x-.5]) - .0001*b + .48*a


def rotated_force(x, y):
    """Pull back coordinates by the inverse rotation, then rotate the force."""
    fx, fy = linear_force(y, -x)
    return np.array([-fy, fx])


def nonlinear_force(x, y):
    """Evaluate the supplied forcing with convection and plate drag included."""
    x, y = x/S, y/S
    return np.array([.02*x**3 + .0098 + .48*x**2, .02*x**2*y - .01 - .96*x*y])


def nonlinear_boundary(x, y):
    """Values on the four supplied edges in normalized coordinates."""
    u, v = np.zeros_like(x), np.zeros_like(y)
    right, bottom, top = x == 1, y == 0, y == 1
    u[right], v[right] = V, -2*V*y[right]
    u[bottom] = V*x[bottom]**2
    u[top], v[top] = V*x[top]**2, -2*V*x[top]
    return np.array([u, v])


def exact_linear(x, y):
    """Return exact Stokes velocity, pressure, and gradient for error evaluation."""
    x, y = x/S, y/S
    return (np.stack((V*f(x)*df(y), -V*df(x)*f(y)), axis=-1), P_STAR*(x-.5)*(y-.5),
            V/S*np.stack((np.stack((df(x)*df(y), f(x)*ddf(y)), axis=-1),
                          np.stack((-ddf(x)*f(y), -df(x)*df(y)), axis=-1)), axis=-2))


def exact_nonlinear(x, y):
    """Return the representable nonlinear fields and velocity gradient."""
    x, y = x/S, y/S
    return (np.stack((V*x**2, -2*V*x*y), axis=-1), P_STAR*(x-y),
            V/S*np.stack((np.stack((2*x, np.zeros_like(x)), axis=-1),
                          np.stack((-2*y, -2*x), axis=-1)), axis=-2))


def pressure_on_rule(rule, full):
    """Interpolate P1 pressure at the integration points, skipping the gauge."""
    grid = rule.grid
    return full[2*grid.nv:-1][grid.p_tri] @ rule.p.T


def integrate(rule, values):
    """Integrate element quadrature samples over the in-plane domain."""
    return float(rule.mesh.det*np.sum(values @ rule.weights))


def field_errors(mesh, full, exact):
    """Use a separate quadrature rule to integrate the squared field errors."""
    rule = Assembly(mesh, order=8)
    velocity, gradient = rule.velocity_on_rule(full)
    pressure = pressure_on_rule(rule, full)
    x, y = np.moveaxis(mesh.physical_points(rule.points), -1, 0)
    u, p, grad = exact(x, y)
    return [float(np.sqrt(integrate(rule, values))) for values in
            (np.sum((velocity-u)**2, axis=-1), (pressure-p)**2,
             np.sum((gradient-grad)**2, axis=(-2, -1)))]


def mesh_metadata(problem):
    """Report resolution and degrees of freedom before and after reduction."""
    mesh, grid = problem.assembly.mesh, problem.assembly.grid
    return {"subdivisions": [grid.nx, grid.ny], "h_max": float(mesh.h_max),
            "triangles": len(grid.ids), "velocity_dofs": 2*grid.nv,
            "pressure_dofs": grid.np, "reduced_dofs_with_gauge": len(problem.rhs)}


def strip_parameters():
    """Return the active-Brinkman strip and its independent physical inputs."""
    return dict(length=.10, height=.020, depth=DEPTH, mu=MU, density=RHO,
                driving=1., gamma=.5, alpha_max=480., q_alpha=2.)


def strip_solve(grid, m, local_terms):
    """Reassemble the linear periodic strip for the supplied parameter set."""
    plate = local_terms.plate_drag_and_partials(m["mu"], m["depth"])
    brinkman = local_terms.brinkman_coefficient_and_partials(m["gamma"], m["alpha_max"], m["q_alpha"])
    rule = Assembly(RectangleMesh(grid, m["length"], m["height"]))
    force = lambda x, y: vector_constant(x, y, m["driving"])
    drag = plate["coefficient"] + brinkman["alpha"]
    problem = Problem(rule, m["mu"], drag, force, periodic=True)
    reduced, error = direct(problem.matrix, problem.rhs)
    return problem, reduced, error, plate, brinkman, force


def axial_peak(grid, full):
    """Find the maximum of the P2 field using nodal and stationary values."""
    reference = np.array([(0, 0), (1, 0), (0, 1), (.5, 0), (.5, .5), (0, .5)])
    x, y = reference.T
    vandermonde = np.column_stack((np.ones(6), x, y, x*x, x*y, y*y))
    nodes = full[:grid.nv]
    coeff = np.linalg.solve(vandermonde, nodes[grid.v_tri].T).T
    peak = nodes.max()
    edges = [(reference[i], reference[j]-reference[i]) for i, j in ((0, 1), (1, 2), (2, 0))]
    for a, bx, by, cxx, cxy, cyy in coeff:
        b = np.array([bx, by])
        hess = np.array([[2*cxx, cxy], [cxy, 2*cyy]])
        candidates = []
        if np.linalg.det(hess) != 0:
            point = np.linalg.solve(hess, -b)
            if np.min(point) >= 0 and np.sum(point) <= 1:
                candidates.append(point)
        for origin, direction in edges:
            curvature = direction @ hess @ direction
            if curvature < 0:
                t = -direction @ (b + hess @ origin)/curvature
                if 0 < t < 1:
                    candidates.append(origin + t*direction)
        for point in candidates:
            peak = max(peak, a + b @ point + point @ hess @ point/2)
    return float(peak)


def strip_boundaries(rule, full, mu, depth):
    """Integrate outward end fluxes, wall shear magnitude, and wall traction."""
    mesh, grid = rule.mesh, rule.grid
    velocity = full[:2*grid.nv].reshape(2, -1).T[grid.v_tri]
    ref = np.array([(0, 0), (1, 0), (0, 1)])
    t, w = line_rule(6)
    left = right = shear = wall_force = 0.
    for i, j in ((0, 1), (1, 2), (2, 0)):
        _, n, grad_ref = basis(ref[i] + t[:, None]*(ref[j]-ref[i]))
        coords = grid.p[grid.p_tri[:, [i, j]]]
        for axis, side in ((0, 0), (0, 1), (1, 0), (1, 1)):
            ids = np.flatnonzero(np.all(coords[..., axis] == side, axis=1))
            normal_sign = 2*side - 1
            if axis == 0:
                # End flux uses the normal velocity and the physical face area.
                flux = depth*(mesh.height/grid.ny)*normal_sign*np.sum((velocity[ids, :, 0] @ n.T) @ w)
                if side == 0:
                    left += flux
                else:
                    right += flux
            else:
                # On each side wall, the shear stress gives the axial traction.
                grad = grad_ref[None, ...] @ mesh.inverse[grid.kind[ids], None, ...]
                strain = np.einsum("eqi,ei->eq", grad[..., 1], velocity[ids, :, 0])
                strain += np.einsum("eqi,ei->eq", grad[..., 0], velocity[ids, :, 1])
                tau = mu*strain
                edge_length = mesh.length/grid.nx
                shear += edge_length*np.sum(np.abs(tau) @ w)
                wall_force += depth*edge_length*normal_sign*np.sum(tau @ w)
    return float(left), float(right), float(shear/(2*mesh.length)), float(wall_force)


def strip_responses(problem, full, m):
    """Return the domain-averaged section flow and the driving power."""
    rule = problem.assembly
    velocity, _ = rule.velocity_on_rule(full)
    q = m["depth"]/m["length"]*integrate(rule, velocity[..., 0])
    return np.array([q, m["driving"]*m["length"]*q])


def measure_strip(solution, m):
    """Collect the strip observables and signed mass, force, and power budgets."""
    problem, state, backward, plate, brinkman, _ = solution
    rule, depth = problem.assembly, m["depth"]
    full = problem.recover(state)
    velocity, grad = rule.velocity_on_rule(full)
    q, input_power = strip_responses(problem, full, m)
    left, right, shear, wall = strip_boundaries(rule, full, m["mu"], depth)
    volume_u = depth*integrate(rule, velocity[..., 0])
    volume_speed_sq = depth*integrate(rule, np.sum(velocity**2, axis=-1))
    strain = (grad + grad.swapaxes(-1, -2))/2
    viscous = depth*integrate(rule, 2*m["mu"]*np.sum(strain**2, axis=(-2, -1)))
    k, alpha = plate["coefficient"], brinkman["alpha"]
    # Positive entries drive the flow and negative entries oppose it.
    force = [m["driving"]*rule.mesh.area*depth, wall, -k*volume_u, -alpha*volume_u]
    power = [input_power, -viscous, -k*volume_speed_sq, -alpha*volume_speed_sq]
    return {"mesh": mesh_metadata(problem), "backward": backward, "gamma": m["gamma"],
            "plate": k, "brinkman": alpha, "flow": q, "peak": axial_peak(rule.grid, full),
            "shear": shear, "input_power": input_power, "flux": [left, right],
            "mass_residual": left+right, "force": force, "force_residual": sum(force),
            "power": power, "dissipation": -sum(power[1:]), "power_residual": sum(power),
            "max_abs_transverse": float(np.max(np.abs(full[rule.grid.nv:2*rule.grid.nv]))),
            "mean_pressure": integrate(rule, pressure_on_rule(rule, full))/rule.mesh.area}


def strip_gradients(grid, m, solution, local_terms, solve_tangent):
    """Compare discrete tangent and adjoint derivatives with fresh solves."""
    from core.sensitivities.implicit import adjoint_directional_derivative, solve_adjoint

    problem, state, _, plate, brinkman, force = solution
    rule, full = problem.assembly, problem.recover(state)
    response = strip_responses(problem, full, m)

    # Differentiate Q = (b/L) integral(u_x) by integrating each P2 basis
    # function. Map that full-space derivative to the free state coordinates.
    element_flow_weights = (m["depth"]/m["length"])*rule.mesh.det*(rule.weights @ rule.n)
    full_flow_partial = np.zeros(rule.grid.size)
    full_flow_partial[:rule.grid.nv] = np.bincount(
        rule.grid.v_tri.ravel(),
        weights=np.tile(element_flow_weights, len(rule.grid.v_tri)),
        minlength=rule.grid.nv,
    )
    flow_state_partial = problem.transpose @ full_flow_partial
    power_state_partial = m["driving"]*m["length"]*flow_state_partial
    jacobian = problem.matrix  # The strip residual is linear in the reduced state.
    flow_adjoint = solve_adjoint(jacobian, flow_state_partial)
    power_adjoint = solve_adjoint(jacobian, power_state_partial)

    drag = plate["coefficient"] + brinkman["alpha"]
    drag_partials = {"mu": plate["d_coefficient_d_viscosity"], "depth": plate["d_coefficient_d_depth"]}
    drag_partials.update({key: brinkman["d_alpha_d_" + key] for key in ("gamma", "alpha_max", "q_alpha")})
    rows = []
    for key, m0 in m.items():
        # Differentiate the assembled residual while holding the state fixed.
        log_l = 1/m0 if key == "length" else 0.
        log_h = 1/m0 if key == "height" else 0.
        log_b = 1/m0 if key == "depth" else 0.
        dforce = lambda x, y: vector_constant(x, y, float(key == "driving"))
        direction = float(key == "mu"), drag_partials.get(key, 0.), log_l, log_h, dforce
        da, df = rule.linear(m["mu"], drag, force, direction)
        partial = problem.transpose @ (da @ full - df)
        tangent = solve_tangent(problem.matrix, partial)
        dq_state = strip_responses(problem, problem.prolong @ tangent, m)[0]
        dq = dq_state + response[0]*(log_b + log_h)
        dp = (m["driving"]*m["length"]*dq + response[1]*log_l
              + response[0]*m["length"]*float(key == "driving"))
        derivative = np.array([dq, dp])

        # The integration measure contributes L and H, while Q's 1/L
        # cancels the direct L contribution. Power also varies through G*L.
        direct_flow = response[0]*(log_b + log_h)
        direct_power = (m["driving"]*m["length"]*direct_flow
                        + response[1]*log_l
                        + response[0]*m["length"]*float(key == "driving"))
        adjoint_derivative = np.array([
            adjoint_directional_derivative(direct_flow, flow_adjoint, partial),
            adjoint_directional_derivative(direct_power, power_adjoint, partial),
        ])
        # Normalize by each response scale, including when the derivative is zero.
        response_scale = m0/np.abs(response)
        tangent_adjoint_difference = np.abs(adjoint_derivative - derivative)*response_scale

        print(f"Checking {key} perturbations", flush=True)
        for epsilon in (1e-2, 1e-3, 1e-4, 1e-5):
            # Re-solve each side of the perturbation on the same connectivity.
            values, errors = [], []
            for sign in (-1, 1):
                perturbed = dict(m, **{key: m0*(1+sign*epsilon)})
                p, z, backward, _, _, _ = strip_solve(grid, perturbed, local_terms)
                values.append(strip_responses(p, p.recover(z), perturbed))
                errors.append(backward)
            delta = m0*epsilon
            fd = (values[1] - values[0])/(2*delta)
            signed_steps = np.array([-delta, delta])[:, None]
            remainder = np.abs(np.asarray(values) - response - signed_steps*derivative)
            adjoint_remainder = np.abs(
                np.asarray(values) - response - signed_steps*adjoint_derivative
            )
            rows.append({"parameter": key, "baseline": m0, "epsilon": epsilon,
                         "derivative": derivative.tolist(), "centered": fd.tolist(),
                         "remainder_minus_plus": remainder.tolist(), "backward_minus_plus": errors,
                         "adjoint": adjoint_derivative.tolist(),
                         "adjoint_remainder_minus_plus": adjoint_remainder.tolist(),
                         "scaled_tangent_adjoint_difference": tangent_adjoint_difference.tolist(),
                         "scaled_adjoint_centered_difference":
                             (np.abs(adjoint_derivative - fd)*response_scale).tolist()})
    return rows


def rotation_comparison(base, rotated, z_base, z_rot):
    """Map the rotated field back and compare nodal values and integrals."""
    original, turned = base.recover(z_base), rotated.recover(z_rot)
    rule, rr = base.assembly, rotated.assembly
    velocity, grad = rule.velocity_on_rule(original)
    velocity_r, grad_r = rr.velocity_on_rule(turned)
    velocity_back = velocity_r @ rr.mesh.rotation
    pressure, pressure_r = pressure_on_rule(rule, original), pressure_on_rule(rr, turned)
    dv, dp = velocity_back-velocity, pressure_r-pressure
    x, y = np.moveaxis(rule.mesh.physical_points(rule.points), -1, 0)
    xr, yr = np.moveaxis(rr.mesh.physical_points(rr.points), -1, 0)
    strain, strain_r = [(g+g.swapaxes(-1, -2))/2 for g in (grad, grad_r)]
    return {"mesh": mesh_metadata(rotated),
            "max_nodal_velocity": float(np.linalg.norm(
                turned[:2*rule.grid.nv].reshape(2, -1).T @ rr.mesh.rotation
                - original[:2*rule.grid.nv].reshape(2, -1).T, axis=1).max()),
            "max_nodal_pressure": float(np.max(np.abs(turned[2*rule.grid.nv:-1]-original[2*rule.grid.nv:-1]))),
            "velocity_l2": float(np.sqrt(integrate(rule, np.sum(dv**2, axis=-1)))),
            "pressure_l2": float(np.sqrt(integrate(rule, dp**2))),
            "area_difference": rr.mesh.area-rule.mesh.area,
            "rms_speed_difference": (np.sqrt(integrate(rr, np.sum(velocity_r**2, axis=-1))/rr.mesh.area)
                                     - np.sqrt(integrate(rule, np.sum(velocity**2, axis=-1))/rule.mesh.area)),
            "force_power_difference": DEPTH*(integrate(rr, np.einsum("aeq,eqa->eq", rotated_force(xr, yr), velocity_r))
                                               - integrate(rule, np.einsum("aeq,eqa->eq", linear_force(x, y), velocity))),
            "viscous_power_difference": 2*MU*DEPTH*(integrate(rr, np.sum(strain_r**2, axis=(-2,-1)))
                                                    - integrate(rule, np.sum(strain**2, axis=(-2,-1)))),
            "plate_power_difference": 480*DEPTH*(integrate(rr, np.sum(velocity_r**2, axis=-1))
                                                  - integrate(rule, np.sum(velocity**2, axis=-1)))}


def jacobian_check(problem, state):
    """Compare the assembled Jacobian action with centered residual differences."""
    direction = np.random.default_rng(317).standard_normal(state.size)
    direction *= np.linalg.norm(state)/np.linalg.norm(direction)
    action = problem.jacobian(state) @ direction
    norm = np.linalg.norm(action)
    rows = []
    for epsilon in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        fd = (problem.residual(state+epsilon*direction)-problem.residual(state-epsilon*direction))/(2*epsilon)
        error = np.linalg.norm(fd-action)
        rows.append([epsilon, float(error), float(error/norm)])
    return {"rows": rows, "action_norm": float(norm), "direction_norm": float(np.linalg.norm(direction))}


def newton_order_study(problem, reference):
    """Measure local order across controlled starts with one Newton step each."""
    grid = problem.assembly.grid
    velocity_end = 2*grid.nv
    free_velocity = np.unique(problem.prolong[:velocity_end].nonzero()[1])
    reference_field = problem.recover(reference)

    def state_error(state):
        # Nondimensionalize the physical fields and omit the gauge multiplier.
        difference = problem.recover(state) - reference_field
        scaled = np.concatenate((difference[:velocity_end]/V,
                                 difference[velocity_end:-1]/P_STAR))
        return float(np.linalg.norm(scaled)/np.sqrt(scaled.size))

    studies = []
    for seed in NEWTON_SEEDS:
        # Perturb only free velocity unknowns to retain the lift and pressure gauge.
        direction = np.zeros_like(reference)
        direction[free_velocity] = np.random.default_rng(seed).standard_normal(free_velocity.size)
        direction *= V*np.sqrt(free_velocity.size)/np.linalg.norm(direction)
        rows = []
        for amplitude in NEWTON_AMPLITUDES:
            start = reference + amplitude*direction
            updated, diagnostics = newton(problem, initial_state=start, max_steps=1)
            before, after = state_error(start), state_error(updated)
            order = None
            full_step = diagnostics["steps"] == [1.0]
            if rows:
                previous = rows[-1]
                if (full_step and previous["full_step"]
                        and previous["after"] > after > NEWTON_ERROR_FLOOR
                        and previous["before"] > before > after):
                    order = float(np.log(previous["after"]/after)
                                  / np.log(previous["before"]/before))
            rows.append({"amplitude": amplitude, "before": before, "after": after,
                         "full_step": full_step, "diagnostics": diagnostics, "order": order})

        # Require two adjacent usable intervals in the quadratic band.
        in_band = [row["order"] is not None and QUADRATIC_BAND[0] <= row["order"] <= QUADRATIC_BAND[1]
                   for row in rows]
        studies.append({"seed": seed, "rows": rows,
                        "passed": any(a and b for a, b in zip(in_band, in_band[1:]))})
    return {"order_band": QUADRATIC_BAND, "state_error_floor": NEWTON_ERROR_FLOOR,
            "required_consecutive_intervals_per_seed": 2, "studies": studies,
            "passed": all(study["passed"] for study in studies)}


def adjoint_check(rows):
    """Apply the declared tolerances to both responses and all nine inputs."""
    selected = [row for row in rows if row["epsilon"] == CENTERED_DIFFERENCE_EPSILON]
    tangent_difference = max(max(row["scaled_tangent_adjoint_difference"]) for row in rows)
    centered_difference = max(max(row["scaled_adjoint_centered_difference"]) for row in selected)

    # Depth is a curved response with measurable Taylor remainders on both sides.
    depth = {row["epsilon"]: row for row in rows if row["parameter"] == "depth"}
    remainder = np.array([depth[epsilon]["adjoint_remainder_minus_plus"]
                          for epsilon in TAYLOR_DEPTH_EPSILONS])
    orders = []
    for i, (coarse, fine) in enumerate(zip(TAYLOR_DEPTH_EPSILONS, TAYLOR_DEPTH_EPSILONS[1:])):
        if np.all(remainder[i:i+2] > 0):
            orders.append((np.log(remainder[i]/remainder[i+1])/np.log(coarse/fine)).tolist())
    taylor_passed = len(orders) == 2 and bool(np.all(
        (QUADRATIC_BAND[0] <= np.asarray(orders)) & (np.asarray(orders) <= QUADRATIC_BAND[1])
    ))
    return {"tangent_adjoint_tolerance": TANGENT_ADJOINT_TOL,
            "centered_difference_epsilon": CENTERED_DIFFERENCE_EPSILON,
            "centered_difference_tolerance": CENTERED_DIFFERENCE_TOL,
            "max_scaled_tangent_adjoint_difference": tangent_difference,
            "max_scaled_adjoint_centered_difference": centered_difference,
            "depth_taylor_epsilons": TAYLOR_DEPTH_EPSILONS,
            "depth_taylor_orders_by_sign_and_response": orders,
            "taylor_order_band": QUADRATIC_BAND,
            "passed": bool(tangent_difference <= TANGENT_ADJOINT_TOL
                           and centered_difference <= CENTERED_DIFFERENCE_TOL and taylor_passed)}


def main():
    """Run the benchmark suite and write one machine-readable results file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.model_root.resolve()))
    from core.sensitivities import local_terms
    from core.sensitivities.implicit import solve_tangent

    result = {"environment": [platform.python_version(), np.__version__, scipy.__version__]}
    # Rungs 1 and 3 share the two straight-strip solves.
    grid, m = RectangleGrid(16, 96), strip_parameters()
    print("Solving the two periodic strips", flush=True)
    active = strip_solve(grid, m, local_terms)
    hollow = dict(m, gamma=1.)
    result["strip"] = [measure_strip(strip_solve(grid, hollow, local_terms), hollow), measure_strip(active, m)]

    result["linear"] = []
    # Refine the Stokes square for Rung 2 and rotate the middle mesh for Rung 4.
    for count in (8, 16, 32):
        print(f"Solving the linear manufactured case on {count} x {count}", flush=True)
        mesh = RectangleMesh(RectangleGrid(count, count), S, S)
        base = Problem(Assembly(mesh), MU, 480, linear_force)
        state, backward = direct(base.matrix, base.rhs)
        result["linear"].append({"mesh": mesh_metadata(base), "backward": backward,
                                  "errors": field_errors(mesh, base.recover(state), exact_linear)})
        if count == 16:
            rotated_mesh = RectangleMesh(mesh.grid, S, S, rotated=True)
            rotated = Problem(Assembly(rotated_mesh), MU, 480, rotated_force)
            rstate, rbackward = direct(rotated.matrix, rotated.rhs)
            result["rotation"] = rotation_comparison(base, rotated, state, rstate)
            result["rotation"]["backward"] = [backward, rbackward]

    print("Solving the nonlinear manufactured case", flush=True)
    mesh = RectangleMesh(RectangleGrid(16, 16), S, S)
    problem = Problem(Assembly(mesh), MU, 480, nonlinear_force, boundary=nonlinear_boundary, density=RHO)
    state, diagnostics = newton(problem)
    result["nonlinear"] = {"mesh": mesh_metadata(problem), "errors": field_errors(mesh, problem.recover(state), exact_nonlinear),
                           "newton": diagnostics, "jacobian": jacobian_check(problem, state)}
    if diagnostics["converged"]:
        result["newton_order"] = newton_order_study(problem, state)
    else:
        result["newton_order"] = {"passed": False, "reason": "Reference Newton solve did not converge."}
    result["gradients"] = strip_gradients(grid, m, active, local_terms, solve_tangent)
    result["adjoint_check"] = adjoint_check(result["gradients"])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0 if result["newton_order"]["passed"] and result["adjoint_check"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
