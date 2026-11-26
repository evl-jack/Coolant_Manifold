"""Bounded graph-strength continuation from the verified beta=1 flat plate."""

import json
import os
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu

from . import project43_sst_bounded_newton as driver
from . import project43_sst_exact_graph as graph

exact = driver.exact


def run(checkpoint, output):
    output.mkdir()
    report = dict(job_id=os.getenv("SLURM_JOB_ID"), classification="FAILED", attempts=[])

    def persist():
        (output/"continuation.json").write_text(json.dumps(report, indent=2, allow_nan=False))

    def callback(model, initial, free, scale, max_steps):
        with np.load(checkpoint, allow_pickle=False) as archive:
            for key, expected in (("initial", initial), ("free", free), ("scale", scale)):
                np.testing.assert_array_equal(archive[key], expected)
            exact.require(int(archive["next_step"]) == 200, "Expected archived sweep 200")
            current, dt = archive["z"], float(archive["dt"])
        solver = driver.Run(model, initial, free, scale)
        beta, decrement = 1., .02

        def configure(value):
            exact.BETA = graph.BETA = value

        def snapshot(state):
            return exact.snapshot(solver.it, solver.edges, state, dt)

        def save():
            np.savez_compressed(output/"last_converged.npz", z=current, beta=beta,
                                initial=initial, free=free, scale=scale, dt=dt)
            report["smallest_converged_beta"] = beta
            persist()

        configure(beta)
        report["start"] = solver.metrics(current, snapshot(current))
        exact.require(solver.admissible(current) and max(report["start"]["stabilized"]) <= 1e-9,
                      "Archived beta=1 residual or admissibility gate failed")
        save()

        def tangent():
            start = perf_counter()
            base = snapshot(current)
            state_matrix = model.jacobian(current)+solver.edges.jacobian(current, base[1])
            jac = (diags(1/scale[free]) @ state_matrix[free][:, free]).tocsc()
            assembly = perf_counter()-start
            start = perf_counter()
            lu = splu(jac)
            factorization = perf_counter()-start
            direction = np.zeros_like(current)
            rhs = -base[0][1][free]/(beta*scale[free])
            start = perf_counter()
            direction[free] = lu.solve(rhs)
            solve = perf_counter()-start
            error = exact.relative_error(rhs, jac @ direction[free])
            exact.require(np.isfinite(direction).all() and error <= 1e-8,
                          f"Continuation tangent linear residual {error}")
            ids = solver.scalars
            ratio = direction[ids]/current[ids]
            limit = .9/float(np.max(ratio)) if np.max(ratio) > 0 else None
            row = int(ids[np.argmax(ratio)])
            window = model.blocks[3] if row >= model.blocks[3].start else model.blocks[2]
            node = row-window.start
            pivots = abs(lu.U.diagonal())
            audit = dict(beta=beta, linear_error=error, assembly_seconds=assembly,
                         factorization_seconds=factorization, solve_seconds=solve,
                         pivot_ratio=float(pivots.min()/pivots.max()), positivity_limit=limit,
                         limiter=dict(node=node, xy=model.points[node].tolist(),
                                      scalar="omega" if window == model.blocks[3] else "k",
                                      value=float(current[row]), tangent=float(direction[row])),
                         max_relative_scalar_tangent=float(np.max(abs(ratio))))
            for name, index, number in (("k1286", 2, 1286), ("omega287", 3, 287)):
                row = model.blocks[index].start+number
                audit[name] = dict(value=float(current[row]), tangent=float(direction[row]),
                                   original=float(base[0][0][row]), graph=float(base[0][1][row]/beta),
                                   stabilized=float(base[0][2][row]), xy=model.points[number].tolist())
            return direction, limit, audit

        for number in range(1, 21):
            configure(beta)
            direction, limit, audit = tangent()
            step = min(decrement, beta-.8, limit if limit is not None else decrement)
            record = dict(attempt=number, tangent=audit, decrement=step, corrections=[])
            report["attempts"].append(record)
            exact.require(step >= 1e-4 or (step == beta-.8 and step > 0),
                          "Positivity or corrector requires a decrement below 1e-4")
            target = max(.8, beta-step)
            record["target_beta"] = target
            trial = current-step*direction
            configure(target)
            try:
                for index in range(9):
                    exact.require(solver.admissible(trial), "Inadmissible predictor or corrector")
                    metrics = solver.metrics(trial, snapshot(trial))
                    record["last_metrics"] = metrics
                    if max(metrics["stabilized"]) <= 1e-9:
                        record["converged"] = True
                        current, beta = trial, target
                        save()
                        if beta <= .8:
                            report["classification"] = "VERIFIED"
                            persist()
                            return current, [], "Continuation complete"
                        decrement = min(.04, 1.25*step)
                        break
                    exact.require(index < 8, "Eight-correction limit reached")
                    info = {}
                    record["corrections"].append(info)
                    base, delta, slope = solver.solve(trial, dt, info)
                    info["line_search"] = search = {}
                    update = solver.search(trial, dt, base, delta, slope, "logarithmic", search)
                    exact.require(update is not None, search.get("failure", "Corrector line search failed"))
                    trial, info["accepted"] = update
            except (ValueError, RuntimeError) as error:
                record["failure"] = str(error)
                decrement = .5*step
                if decrement < 1e-4:
                    report["stop"] = "Corrector failure requires a decrement below 1e-4"
                    persist()
                    return current, [], report["stop"]
            persist()
            print(json.dumps(record), flush=True)
        report["stop"] = "Continuation attempt budget exhausted"
        persist()
        return current, [], report["stop"]

    try:
        report["checkpoint_sha256"] = driver.fingerprint(checkpoint)
        exact.require(report["checkpoint_sha256"] == exact.CHECKPOINT_HASH, "Checkpoint hash mismatch")
        exact.require(np.finfo(np.longdouble).eps < np.finfo(float).eps, "Extended precision unavailable")
        exact.segregated.plate.pseudo_transient = callback
        exact.segregated.plate.run()
    except (ValueError, RuntimeError, OSError, AssertionError) as error:
        report["stop"] = str(error)
    persist()
    return report["classification"] == "VERIFIED"


if __name__ == "__main__":
    raise SystemExit(0 if run(*map(Path, sys.argv[1:])) else 2)
