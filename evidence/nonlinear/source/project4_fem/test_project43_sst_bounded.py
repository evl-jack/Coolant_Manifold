"""Fixed-checkpoint verification executed before any nonlinear sweeps."""

import numpy as np


def gate_a(iteration, state, dt, report):
    raw, target = iteration.residual(state, dt)
    report["original"] = iteration.norms(raw)
    report["target"] = iteration.norms(target)
    report["cf"] = float(iteration.model.cf_at(state))
    report["dt"] = dt
    np.testing.assert_allclose(report["target"][2:], [1.003e-7, 2.859e-5], rtol=.005)
    report["scalars"] = []
    for j in range(2):
        _, _, audit = iteration.scalar(state, j, dt)
        report["scalars"].append(audit)
        np.testing.assert_allclose(audit["relaxation"], .7, atol=1e-14)
        expected = [1.44645e-8, 6.24443][j]
        # Compare the undamped proposal and let the audit retain both minima.
        np.testing.assert_allclose(audit["full_step_minimum"], expected, rtol=.01)
        np.testing.assert_allclose(audit["positivity_nodes"], [169, 14][j], atol=2, rtol=0)
    raw_after, target_after = iteration.residual(state, dt)
    np.testing.assert_array_equal(raw_after, raw)
    np.testing.assert_array_equal(target_after, target)
    report["steady_residual_unchanged"] = True
    trial = state.copy()
    for j in range(2):
        ids, candidate, _ = iteration.scalar(trial, j, dt)
        trial[ids] = candidate
    sequential = iteration.norms(iteration.residual(trial, dt)[1])
    report["sequential_target"] = sequential
    if sequential[3] > .9*report["target"][3]:
        raise ValueError("Sequential trial reduced target omega residual by less than 10%")
