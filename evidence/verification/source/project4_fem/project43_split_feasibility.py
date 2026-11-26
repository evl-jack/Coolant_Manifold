"""Two-control, fixed-envelope 50/50 split feasibility test on the frozen mesh."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import minimize
from .project43_baseline import Geometry
from .project43_brinkman import create_problem


def run(out):
    out.mkdir()
    g = Geometry()
    fem, base, fixed, values, _, flux, _ = create_problem(g, .5)
    xy = fem.points[fem.triangles].mean(axis=1) * g.height_mm
    x, y = xy.T
    # 1 mm long throttle inside the downstream main passage; open two adjacent
    # 1 mm shoulders of the branch over its first 2 mm, below the fixed exit.
    throttle = (x > 13) & (x < 14) & (y < 4)
    shoulders = (abs(x-10) > .5) & (abs(x-10) < 1.5) & (y > 4) & (y < 6)
    directions = np.stack((-throttle.astype(float), shoulders.astype(float)), axis=-1)[:, None, :]
    area_gradient = np.einsum('eq,eld->d', fem.weights, directions)*g.height_mm**2
    np.testing.assert_allclose(area_gradient, [-4, 4])
    alpha = 1e6 * .004**2 / .001002
    fem.inertia = 50.
    load = fem.load()
    cache = {}
    history = []
    best = None
    def evaluate(design):
        nonlocal best
        key = tuple(design)
        if key in cache:
            return cache[key]
        gamma = base + (directions @ design)
        resistance = alpha*(1-gamma)**2
        try:
            state = fem.solve(fixed, values, load, resistance=resistance)
        except RuntimeError:
            initial = None
            for inertia in (0., 12.5, 25., 37.5, 50.):
                fem.inertia = inertia
                state = fem.solve(fixed, values, load, resistance=resistance, initial=initial)
                initial = state.values
        inlet, main, side = flux @ state.values
        imbalance = abs(main+side-inlet)/inlet
        if min(inlet, main, side) <= 0 or imbalance > 1e-9:
            raise RuntimeError('Signed port flow or mass balance failed')
        share = side/inlet
        partial = (flux[2]-share*flux[0])/inlet
        gradient = fem.drag_adjoint(state, [partial], -2*alpha*(1-gamma[:,:,None])*directions)[0]
        area = 83. + area_gradient @ design
        record = dict(design=list(key), fluid_area_mm2=area, side_share=share, main_share=main/inlet,
                      share_gradient=gradient.tolist(), residual=state.residual_history[-1], imbalance=imbalance)
        history.append(record)
        if area <= 83.+1e-9 and (best is None or abs(share-.5)<abs(best['side_share']-.5)):
            best = record
            np.savez_compressed(out/'best.npz', vertices_mm=fem.points*g.height_mm,
                                triangles=fem.triangles, gamma=gamma[:,0], state=state.values)
        (out/'history.json').write_text(json.dumps(history, indent=2))
        value = (.5*(share-.5)**2, (share-.5)*gradient)
        cache[key] = value
        return value
    evaluate(np.zeros(2))
    result = minimize(evaluate, [.08,.04], jac=True, method='SLSQP', bounds=[(0,1)]*2,
                      constraints=[dict(type='ineq', fun=lambda z: -area_gradient@z,
                                        jac=lambda z: -area_gradient)],
                      options=dict(ftol=1e-13, maxiter=35))
    report = dict(optimizer_success=bool(result.success), message=str(result.message), best=best,
                  target_met=bool(abs(best['side_share']-.5)<=.001), split_tolerance_fraction=.001,
                  interpretation='Relaxed density feasibility within two fixed control regions.')
    (out/'report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report))


if __name__ == '__main__':
    run(Path(sys.argv[1]))
