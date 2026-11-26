"""Check the Re=50 laminar density adjoint for an opening of both branch shoulders."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import expit

from .project43_baseline import Geometry
from .project43_brinkman import create_problem
from .project43_power import port_power


def run(wall_report):
    candidate = wall_report["best_feasible"]
    width = candidate["shoulder_width_mm"]
    throat = candidate["open_height_mm"]
    g = Geometry()
    spacing = .25
    fem, base, fixed, values, ports, flux, _ = create_problem(g, spacing)
    x, y = fem.quad_points.transpose(2, 0, 1) * g.height_mm
    shoulders = (y > 4) & (y < 6)
    baffle = (x > 13) & (x < 14)
    mu, rho, height, re = .001002, 998.2, .004, 50.
    speed = re*mu/(rho*height)
    power_scale = mu*speed**2
    alpha = 1.e6*height**2/mu
    transition = spacing/2
    load = fem.load()

    def resistance(shoulder):
        fraction = expit((.5+shoulder-abs(x-10))/transition)
        gamma = np.where(shoulders, fraction, base)
        gamma = np.where(baffle & (y < 4), expit((throat-y)/transition), gamma)
        derivative = np.where(shoulders, fraction*(1-fraction)/transition, 0.)
        return alpha*(1-gamma)**2, -2*alpha*(1-gamma)*derivative

    def responses(state):
        inlet, main, side = flux @ state.values
        if min(inlet, main, side) <= 0 or abs(main+side-inlet)/inlet > 1.e-9:
            raise ValueError("Outlet flux check failed")
        power, grad = port_power(fem, ports, state.values, re)
        share = side/inlet
        return (np.array([sum(power)*power_scale, share]),
                np.stack([grad*power_scale, (flux[2]-share*flux[0])/inlet]))

    drag, derivative = resistance(width)
    state = None
    for inertia in (0., re/2, re):
        fem.inertia = inertia
        state = fem.solve(fixed, values, load, resistance=drag,
                          initial=None if state is None else state.values)
    baseline, partials = responses(state)
    adjoint = fem.drag_adjoint(state, partials, derivative[:,:,None])[:,0]
    checks = []
    for step in (.002, .001):
        pair = []
        for sign in (-1, 1):
            shifted, _ = resistance(width+sign*step)
            trial = fem.solve(fixed, values, load, resistance=shifted, initial=state.values)
            pair.append(responses(trial)[0])
        centered = (pair[1]-pair[0])/(2*step)
        checks.append(dict(step_mm=step, centered=centered.tolist(),
                           relative_difference=(abs(adjoint-centered)/np.maximum(abs(adjoint),1e-16)).tolist()))
    passed = all(max(check["relative_difference"]) < 1e-3 for check in checks)
    return dict(re=re, spacing_mm=spacing, transition_scale_mm=transition,
                alpha_max_kg_m3_s=1e6, shoulder_width_mm=width, throat_mm=throat,
                quantities=["port_power_W_per_m", "side_share_fraction"], derivative_per="mm of shoulder width",
                explicit_wall=[candidate["pumping_power_W_per_m"], candidate["side_share"]],
                density_baseline=baseline.tolist(), adjoint=adjoint.tolist(), checks=checks, passed=passed)


if __name__ == "__main__":
    result = run(json.loads(Path(sys.argv[1]).read_text()))
    with Path(sys.argv[2]).open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
    print(json.dumps(result))
    sys.exit(0 if result["passed"] else 2)
