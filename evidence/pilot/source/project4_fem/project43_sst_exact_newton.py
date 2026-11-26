"""Two-checkpoint derivative verification and discarded Newton trials."""

import hashlib
import json
import os
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu

from .project43_sst_bounded import BETA, CHECKPOINT_HASH, Iteration, require, segregated
from .project43_sst_coupled_probe import relative_error
from .project43_sst_exact_graph import GraphEdges

AMPLITUDES = (1e-4, 1e-5, 1e-6, 1e-7)
NAMES = ("velocity", "pressure", "k", "omega")


def snapshot(iteration, edges, z, dt):
    raw = iteration.model.residual(z)
    graph = np.zeros_like(raw)
    operators = [iteration.operator(z, j, dt) for j in range(2)]
    for block, (_, matrix, _, _) in zip(iteration.model.blocks[2:], operators):
        graph[block] = BETA*(matrix @ z[block])
    flags = iteration.fields(z)[-1]
    branches = [edges.branches(item[0]) for item in operators]
    return (raw, graph, raw+graph), operators, flags, branches


def audit(iteration, edges, state, dt, base, matrices, rng, report):
    free, scale = iteration.free, iteration.scale
    states, _, flags, branches = base
    names = ("original", "graph", "complete")
    report["graph_ties"] = [int(np.count_nonzero(tied)) for _, tied in branches]
    all_pass = True
    for name, ids in zip((*NAMES, "mixed"), (*iteration.ids, free)):
        direction = np.zeros_like(state)
        # Small relative scalar directions keep all four requested amplitudes usable.
        magnitude = np.ones_like(state)
        for scalar in iteration.ids[2:]:
            magnitude[scalar] = .01*state[scalar]
        direction[ids] = rng.uniform(-1, 1, len(ids))*magnitude[ids]
        actions = [matrix @ direction[free] for matrix in matrices]
        records = []
        streak = np.zeros((3, 4), dtype=int)
        passed = np.zeros_like(streak, dtype=bool)
        for h in AMPLITUDES:
            sides = [snapshot(iteration, edges, state + sign*h*direction, dt)
                     for sign in (1, -1)]
            changes = []
            for _, _, other_flags, other_branches in sides:
                changes.append(dict(
                    sst=sum(int(np.count_nonzero(a != b)) for a, b in zip(flags, other_flags)),
                    graph=[int(np.count_nonzero(a != b))
                           for (a, _), (b, _) in zip(branches, other_branches)]))
            row = dict(amplitude=h, changes=changes, errors={})
            farther = {}
            for j, (kind, action) in enumerate(zip(names, actions)):
                diagnostics = []
                for block, block_ids in enumerate(iteration.ids):
                    p = np.searchsorted(free, block_ids)
                    w = scale[block_ids]
                    center = states[j][block_ids]/w
                    plus, minus = [side[0][j][block_ids]/w for side in sides]
                    changed = [
                        bool((j != 1 and block != 1 and m["sst"]) or
                             (j != 0 and block >= 2 and
                              (m["graph"][block-2] or m["sst"])))
                        for m in changes]
                    difference = (plus-minus)/(2*h)
                    method = "centered"
                    if any(changed):
                        # A second-order one-sided derivative stays on the base branch.
                        side = 0 if not changed[0] else 1
                        sign = 1 if side == 0 else -1
                        if side not in farther:
                            farther[side] = snapshot(iteration, edges, state + sign*2*h*direction, dt)
                        far = farther[side]
                        near = plus if side == 0 else minus
                        difference = sign*(-3*center+4*near-far[0][j][block_ids]/w)/(2*h)
                        method = "forward" if side == 0 else "backward"
                    err = float(np.linalg.norm(difference-action[p]))
                    norm = float(np.linalg.norm(action[p]))
                    cancellation = float(64*np.finfo(float).eps*
                                         np.linalg.norm(abs(plus)+abs(minus)+abs(center))/h)
                    zero = norm <= max(1e-12, cancellation*1e4)
                    relative = None if zero else relative_error(difference, action[p])
                    ok = err <= max(1e-10, cancellation) if zero else relative <= 1e-4
                    streak[j, block] = streak[j, block]+1 if ok else 0
                    passed[j, block] |= streak[j, block] >= 2
                    diagnostics.append(dict(relative=relative, absolute=err, action_norm=norm,
                                            cancellation=cancellation, method=method))
                row["errors"][kind] = diagnostics
            records.append(row)
        report[name] = dict(amplitudes=records, passed=passed.tolist())
        all_pass &= bool(passed.all())
    return all_pass


def node_diagnostics(iteration, z, dt, delta, target):
    m = iteration.model
    result = {}
    for j in range(2):
        matrix, graph, ids, local = iteration.operator(z, j, dt)
        system = (matrix+graph)[local][:, local]
        base_diagonal = system.diagonal()
        q = z[ids]
        sink = np.zeros_like(q)
        if j:
            fields = iteration.fields(z)
            k, omega = fields[5]
            grad_k, grad_omega = fields[6]
            distance, weights = m.geometry[2:4]
            cross = 2*.856*np.sum(grad_k*grad_omega, axis=-1)/omega
            near = np.sqrt(k)/(.09*omega*distance)
            viscous = 500*m.nu/(distance**2*omega)
            cap = 4*.856*k/(np.maximum(cross, 1e-10)*distance**2)
            blend = np.tanh(np.minimum(np.maximum(near, viscous), cap)**4)
            assembled = np.zeros(m.nb)
            np.add.at(assembled, m.triangles, (weights*np.maximum(-(1-blend)*cross, 0)) @ m.p)
            sink = assembled[local]/q
        positivity = np.maximum(0, target[ids]-system @ q-sink*q)/q
        _, candidate, audit_row = iteration.scalar(z, j, dt)
        selected = {"minimum_k": int(np.argmin(q))} if not j else {
            "maximum_scaled_omega_residual": int(np.argmax(abs(target[ids]/iteration.scale[ids]))),
            "upstream_omega_287": int(np.flatnonzero(local == 287)[0])}
        for label, i in selected.items():
            result[label] = dict(node=int(local[i]), xy=m.points[local[i]].tolist(),
                                 scalar=float(q[i]), scaled_target=float(target[ids[i]]/iteration.scale[ids[i]]),
                                 coupled_delta=float(delta[ids[i]]),
                                 segregated_delta=float((candidate[i]-q[i])/audit_row["relaxation"]),
                                 segregated_relaxation=audit_row["relaxation"],
                                 original_iteration_diagonal=float(base_diagonal[i]),
                                 sink_ratio=float(sink[i]/base_diagonal[i]),
                                 positivity_ratio=float(positivity[i]/base_diagonal[i]))
        result[NAMES[j+2]+"_maximum_diagonal_ratios"] = {
            "sink": float(np.max(sink/base_diagonal)),
            "positivity": float(np.max(positivity/base_diagonal))}
    return result


def trial_steps(iteration, edges, z, initial, dt, delta, target, report):
    m, free, scale = iteration.model, iteration.free, iteration.scale
    fixed = np.setdiff1d(np.arange(len(z)), free)
    ids = np.r_[iteration.ids[2], iteration.ids[3]]
    ratios = delta[ids]/z[ids]
    divisor = max(1., float(max(abs(ratios))))
    direction = delta/divisor
    f = target[free]/scale[free]
    merit = .5*float(f @ f)
    slopes = []
    for h in AMPLITUDES:
        def at(a):
            return snapshot(iteration, edges, z+a*direction, dt)[0][2][free]/scale[free]
        slope = float(f @ ((-3*f+4*at(h)-at(2*h))/(2*h)))*divisor
        slopes.append(slope)
    report["merit_slopes"] = slopes
    report["merit"] = merit
    report["linearized_merit_slope"] = -2*merit
    report["merit_slope_verified"] = bool(
        all(s < 0 for s in slopes[-2:]) and abs(slopes[-2]-slopes[-1]) <= 1e-4*abs(slopes[-1]))
    require(report["merit_slope_verified"], "No consistent negative forward merit derivative")
    maximum = max(iteration.norms(target))
    for mapping, index in (("additive", np.argmin(ratios)), ("logarithmic", np.argmax(abs(ratios)))):
        ratio = float(ratios[index])
        first = min(1., -.9/ratio) if mapping == "additive" and ratio < 0 else 1.
        if mapping == "logarithmic" and ratio != 0:
            first = min(1., 2/abs(ratio))
        row = int(ids[index])
        scalar = int(row >= m.blocks[3].start)
        local = row-m.blocks[2+scalar].start
        records = []
        start = perf_counter()
        for n in range(16):
            step = first*.5**n
            trial = z+step*delta
            if mapping == "logarithmic":
                trial[ids] = z[ids]*np.exp(step*ratios)
            good = bool(np.isfinite(trial).all() and np.all(trial[ids] > 0) and
                        np.array_equal(trial[fixed], initial[fixed]))
            record = dict(step=step, valid=good, promising=False)
            if good:
                nut = iteration.fields(trial)[7]
                if np.isfinite(nut).all() and np.all(nut >= 0):
                    residual = iteration.residual(trial, dt)[1]
                    blocks = iteration.norms(residual)
                    r = residual[free]/scale[free]
                    phi = .5*float(r @ r)
                    armijo = phi <= merit+1e-4*step*slopes[-1]
                    record.update(blocks=blocks, merit=phi, armijo=armijo,
                                  merit_reduction=1-phi/merit,
                                  max_reduction=1-max(blocks)/maximum,
                                  minima=[float(trial[i].min()) for i in iteration.ids[2:]],
                                  promising=armijo and phi <= .99*merit and max(blocks) <= maximum)
                else:
                    record["valid"] = False
            records.append(record)
            if record["promising"]:
                break
        report[mapping] = dict(trials=records, seconds=perf_counter()-start,
                               limiter=dict(scalar=NAMES[2+scalar], node=local,
                                            xy=m.points[local].tolist(), ratio=ratio))
    return any(trial["promising"] for mapping in ("additive", "logarithmic")
               for trial in report[mapping]["trials"])


def run(paths, evidence, output):
    report = dict(job_id=os.getenv("SLURM_JOB_ID"), checkpoints={})
    archived = json.loads(evidence.read_text())
    rng = np.random.default_rng(43)

    def callback(model, initial, free, scale, max_steps):
        iteration = Iteration(model, free, scale)
        edges = GraphEdges(iteration)
        prepared = []
        for step, path, expected in zip((200, 217), paths,
                                        (archived["fixed_state"], archived["history"][-1])):
            row = report["checkpoints"][str(step)] = {}
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            row["sha256"] = digest
            if step == 200:
                require(digest == CHECKPOINT_HASH, "Sweep-200 checkpoint hash mismatch")
            with np.load(path, allow_pickle=False) as cp:
                for name, reference in (("initial", initial), ("free", free), ("scale", scale)):
                    np.testing.assert_array_equal(cp[name], reference)
                require(int(cp["next_step"]) == step, "Wrong checkpoint sweep")
                z, dt = cp["z"], float(cp["dt"])
            start = perf_counter()
            base = snapshot(iteration, edges, z, dt)
            raw, _, target = base[0]
            for name, residual in (("original", raw), ("target", target)):
                row[name] = iteration.norms(residual)
                np.testing.assert_allclose(row[name], expected[name], rtol=.005, atol=2e-14)
            original = model.jacobian(z)
            graph = edges.jacobian(z, base[1])
            matrix = [diags(1/scale[free]) @ item[free][:, free]
                      for item in (original, graph, original+graph)]
            row["assembly_seconds"] = perf_counter()-start
            left, right = rng.standard_normal((2, len(free)))
            row["transpose_errors"] = [
                abs(float(left @ (a @ right)-(a.T @ left) @ right))/
                max(1., abs(float(left @ (a @ right)))) for a in matrix]
            require(max(row["transpose_errors"]) <= 1e-10, "Transpose identity failed")
            row["derivatives"] = {}
            row["verified"] = audit(iteration, edges, z, dt, base, matrix, rng, row["derivatives"])
            prepared.append((z, dt, row, matrix[-1], target))
        require(all(row["verified"] for _, _, row, _, _ in prepared), "Jacobian verification gate failed")
        for z, dt, row, matrix, target in prepared:
            start = perf_counter()
            lu = splu(matrix.tocsc())
            row["factorization_seconds"] = perf_counter()-start
            start = perf_counter()
            delta = np.zeros_like(z)
            rhs = -target[free]/scale[free]
            delta[free] = lu.solve(rhs)
            row["solve_seconds"] = perf_counter()-start
            row["linear_error"] = relative_error(rhs, matrix @ delta[free])
            require(np.isfinite(delta).all() and row["linear_error"] <= 1e-8, "Linear gate failed")
            row["nodes"] = node_diagnostics(iteration, z, dt, delta, target)
            row["feasible"] = trial_steps(iteration, edges, z, initial, dt, delta, target, row)
        require(all(row["feasible"] for _, _, row, _, _ in prepared), "Newton feasibility gate failed")
        return z, [], "Diagnostic only"

    try:
        segregated.plate.pseudo_transient = callback
        segregated.plate.run()
    except (OSError, ValueError, RuntimeError, AssertionError) as exc:
        report["failure"] = str(exc)
    output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key != "checkpoints"}), flush=True)
    return "failure" not in report


if __name__ == "__main__":
    first, second, archive, output = map(Path, sys.argv[1:])
    raise SystemExit(0 if run((first, second), archive, output) else 2)
