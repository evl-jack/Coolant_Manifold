"""Read-only coupled Newton probe at the archived sweep-200 state."""

import argparse
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.sparse import block_diag, csr_matrix, diags
from scipy.sparse.linalg import splu

from .project43_sst_bounded import BETA, CHECKPOINT_HASH, Iteration, require, segregated


def relative_error(actual, approx):
    return float(np.linalg.norm(actual - approx) /
                 max(np.linalg.norm(actual), np.finfo(float).tiny))


def probe(model, state, initial, free, scale, dt, report):
    iteration = Iteration(model, free, scale)
    fixed = np.setdiff1d(np.arange(len(state)), free)
    scalar_ids = np.r_[iteration.ids[2], iteration.ids[3]]
    q = state[scalar_ids]
    require(np.all(q > 0) and np.array_equal(state[fixed], initial[fixed]),
            "Invalid checkpoint scalar or prescribed trace")

    def target(z):
        return iteration.residual(z, dt)[1]

    def merit(f):
        return float(.5 * (f @ f))

    start = perf_counter()
    steady = target(state)
    initial_blocks = iteration.norms(steady)
    report["initial_blocks"] = initial_blocks
    np.testing.assert_allclose(initial_blocks,
                               [3.28865e-14, 2.89787e-17, 1.002976e-7, 2.859111e-5],
                               rtol=.005, atol=0)
    f = steady[free]/scale[free]
    initial_merit = merit(f)
    graphs = [iteration.operator(state, j, dt)[1] for j in range(2)]
    coupled = model.jacobian(state) + BETA*block_diag(
        [csr_matrix((model.fem.size, model.fem.size)), *graphs], format="csr")
    jacobian = (diags(1/scale[free]) @ coupled[free][:, free]).tocsc()
    report["assembly_seconds"] = perf_counter() - start
    start = perf_counter()
    lu = splu(jacobian)
    report["factorization_seconds"] = perf_counter() - start
    start = perf_counter()
    delta = np.zeros_like(state)
    delta[free] = lu.solve(-f)
    report["solve_seconds"] = perf_counter() - start
    linear = relative_error(-f, jacobian @ delta[free])
    report["relative_linear_residual"] = linear
    require(np.isfinite(delta).all() and linear <= 1e-8, "Linear solve gate failed")

    ratios = delta[scalar_ids]/q
    largest = float(np.max(abs(ratios)))
    divisor = max(1., largest)
    direction = delta/divisor
    approx = jacobian @ direction[free]
    block_positions = [np.searchsorted(free, ids) for ids in iteration.ids]
    audit = []
    derivatives = []
    start = perf_counter()
    for h in (1e-4, 1e-5, 1e-6):
        plus = target(state + h*direction)[free]/scale[free]
        minus = target(state - h*direction)[free]/scale[free]
        actual = (plus-minus)/(2*h)
        derivatives.append(actual)
        audit.append(dict(amplitude=h, discrepancy=relative_error(actual, approx),
                          block_discrepancy=[relative_error(actual[p], approx[p])
                                             for p in block_positions],
                          merit_derivative=float(f @ actual)))
    report["audit_seconds"] = perf_counter() - start
    report["normalization_divisor"] = divisor
    report["directional_audit"] = audit
    report["fd_relative_change_to_finest"] = [relative_error(derivatives[-1], d)
                                               for d in derivatives[:-1]]
    descent = all(item["merit_derivative"] < 0 for item in audit)
    slope = audit[-1]["merit_derivative"]*divisor
    report["merit_derivative_newton_scale"] = slope
    report["initial_merit"] = initial_merit

    def node(index):
        dof = int(scalar_ids[index])
        j = int(dof >= model.blocks[3].start)
        local = dof - model.blocks[j+2].start
        return dict(scalar=("k", "omega")[j], node=local, dof=dof,
                    coordinates=model.points[local].tolist(), relative_direction=float(ratios[index]))

    additive_limit = int(np.argmin(ratios))
    log_limit = int(np.argmax(abs(ratios)))
    for mapping, index, first in (
        ("additive", additive_limit, min(1., -.9/ratios[additive_limit])
         if ratios[additive_limit] < 0 else 1.),
        ("logarithmic", log_limit, min(1., 2/largest) if largest else 1.),
    ):
        trials = []
        start = perf_counter()
        for backtrack in range(8):
            step = first*.5**backtrack
            z = state + step*delta
            if mapping == "logarithmic":
                z[scalar_ids] = q*np.exp(step*ratios)
            require(np.array_equal(z[fixed], initial[fixed]) and np.all(z[scalar_ids] > 0)
                    and np.isfinite(z).all(), "Invalid trial state")
            value = target(z)
            blocks = iteration.norms(value)
            trial_merit = merit(value[free]/scale[free])
            armijo = bool(trial_merit <= initial_merit + 1e-4*step*slope) if descent else None
            record = dict(step=float(step), blocks=blocks, merit=trial_merit, armijo=armijo,
                          merit_reduction=1-trial_merit/initial_merit,
                          max_residual_reduction=1-max(blocks)/max(initial_blocks),
                          minima=[float(z[ids].min()) for ids in iteration.ids[2:]])
            trials.append(record)
            if armijo:
                break
        elapsed = perf_counter()-start
        final = trials[-1] if trials[-1]["armijo"] else min(trials, key=lambda x: x["merit"])
        report[mapping] = dict(initial_step=float(first), limiting_node=node(index),
                               trial_seconds=elapsed, trials=trials, final=final,
                               promising=bool(final["merit_reduction"] >= .01 and
                                              final["max_residual_reduction"] >= 0))


def run(checkpoint, output):
    report = dict(sol_job_id=os.getenv("SLURM_JOB_ID"), beta=BETA)

    def callback(model, initial, free, scale, max_steps):
        with np.load(checkpoint, allow_pickle=False) as cp:
            for name, expected in (("initial", initial), ("free", free), ("scale", scale)):
                np.testing.assert_array_equal(cp[name], expected)
            require(int(cp["next_step"]) == 200, "Expected sweep 200")
            state, dt = cp["z"], float(cp["dt"])
        probe(model, state, initial, free, scale, dt, report)
        return state, [], "Read-only diagnostic"

    try:
        report["checkpoint_sha256"] = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
        require(report["checkpoint_sha256"] == CHECKPOINT_HASH, "Checkpoint hash mismatch")
        segregated.plate.pseudo_transient = callback
        segregated.plate.run()
    except (OSError, ValueError, AssertionError, RuntimeError) as exc:
        report["failure"] = str(exc)
    output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(report, allow_nan=False))
    return "failure" not in report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    raise SystemExit(0 if run(args.checkpoint, args.output) else 2)
