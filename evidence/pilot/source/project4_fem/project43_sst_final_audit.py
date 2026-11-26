"""Resolve the remaining k-only derivative block and reuse the archived audit."""

import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

from . import project43_sst_exact_newton as newton


def isolated_comparison(iteration, edges, state, dt, base, matrix, rng, report):
    start = perf_counter()
    m = iteration.model
    velocity_ids, pressure_ids, k_ids, omega_ids = iteration.ids
    for ids in (velocity_ids, pressure_ids):
        rng.uniform(-1, 1, len(ids))
    direction = np.zeros_like(state)
    direction[k_ids] = rng.uniform(-1, 1, len(k_ids))*(.01*state[k_ids])
    action = (matrix @ direction[iteration.free])[np.searchsorted(iteration.free, omega_ids)]
    budget = 1e-4*np.linalg.norm(action)
    precision = np.longdouble
    eps = np.finfo(precision).eps
    report.update(longdouble_epsilon=float(eps), longdouble_bits=np.finfo(precision).nmant,
                  float64_epsilon=float(np.finfo(float).eps), budget=float(budget), amplitudes=[])
    newton.require(eps < np.finfo(float).eps, "Sol longdouble has no precision advantage")
    geometry = tuple(np.asarray(a, dtype=precision) for a in m.geometry)
    _, gradients, _, weights, h, _ = geometry
    coefficients = state[m.dofs].astype(precision)
    perturbation = direction[m.dofs].astype(precision)
    velocity = m.fields(coefficients, geometry)[0]
    along = np.einsum("eqi,eai->eqa", velocity, gradients)
    speed_term = 4*np.sum(velocity*velocity, axis=2)/h**2

    def terms(c):
        diffusion = m.fields(c, geometry)[8][1]
        tau = 1/np.sqrt(speed_term+(4*diffusion/h**2)**2)
        return (np.einsum("eq,eai,ebi->eab", weights*diffusion, gradients, gradients),
                np.einsum("eq,eqa,eqb->eab", weights*tau, along, along))

    baseline = terms(coefficients)
    i, j = edges.i, edges.j
    winners, ties = base[3][1]
    newton.require(not np.any(ties), "Base omega graph has exact ties")
    keys = np.r_[i*m.nb+j, j*m.nb+i]
    permutation = np.argsort(keys)
    sorted_keys = keys[permutation]
    row, col = np.broadcast_arrays(m.triangles[:, :, None], m.triangles[:, None, :])
    mask = row != col
    locations = permutation[np.searchsorted(sorted_keys, (row*m.nb+col)[mask])]
    winners_flat = (winners-1)*len(i)+np.arange(len(i))
    active = winners != 0
    q = state[m.blocks[3]].astype(precision)
    jump = q[i]-q[j]
    local_ids = omega_ids-m.blocks[3].start
    scale = iteration.scale[omega_ids].astype(precision)
    fixed = np.setdiff1d(np.arange(len(state)), iteration.free)
    all_scalars = np.r_[k_ids, omega_ids]
    original = base[1][1][0]
    original_entries = np.stack((np.asarray(original[i, j]).ravel(),
                                 np.asarray(original[j, i]).ravel())).astype(precision)

    def assemble(values):
        out = np.zeros(2*len(i), dtype=precision)
        np.add.at(out, locations, values[mask])
        return out

    def sample(offset):
        z = state+offset*direction
        valid = (np.isfinite(z).all() and np.all(z[all_scalars] > 0)
                 and np.array_equal(z[fixed], state[fixed]))
        if not valid:
            return None, dict(admissible=False)
        fields = iteration.fields(z)
        matrix_at_z = iteration.operator(z, 1, dt)[0]
        branches, sample_ties = edges.branches(matrix_at_z)
        flags = sum(int(np.count_nonzero(a != b)) for a, b in zip(base[2], fields[-1]))
        changes = int(np.count_nonzero(branches != winners))
        valid = bool(np.isfinite(fields[7]).all() and np.all(fields[7] >= 0))
        current = terms(coefficients+precision(offset)*perturbation)
        differences = assemble(sum(a-b for a, b in zip(current, baseline)))
        error = assemble(sum(abs(a)+abs(b) for a, b in zip(current, baseline)))
        perturbed = original_entries+differences.reshape(2, -1)
        precise_winners = np.argmax(np.vstack((np.zeros(len(i), dtype=precision), perturbed)), axis=0)
        precise_changes = int(np.count_nonzero(precise_winners != winners))
        info = dict(admissible=valid, sst_changes=flags, graph_changes=changes,
                    precise_graph_changes=precise_changes, graph_ties=int(np.count_nonzero(sample_ties)))
        if not valid or flags or changes or precise_changes or np.any(sample_ties):
            return None, info
        return (differences[winners_flat[active]], error[winners_flat[active]]), info

    def residual_edges(values, bound=False):
        out = np.zeros(m.nb, dtype=precision)
        np.add.at(out, i[active], values)
        np.add.at(out, j[active], values if bound else -values)
        return out[local_ids]/scale

    streak = 0
    passed = False
    for amplitude in (5e-4, 2e-4, *newton.AMPLITUDES):
        samples = {s: sample(s*amplitude) for s in (1, -1)}
        clean = [s for s in samples if samples[s][0] is not None]
        stencil = None
        if len(clean) == 2:
            stencil = {1: .5, -1: -.5}
        elif clean:
            s = clean[0]
            samples[2*s] = sample(2*s*amplitude)
            if samples[2*s][0] is not None:
                stencil = {s: 2*s, 2*s: -.5*s}
        row = dict(amplitude=amplitude, samples={str(s): v[1] for s, v in samples.items()})
        if stencil:
            flux = sum(w*samples[s][0][0] for s, w in stencil.items())/precision(amplitude)
            bound = 8*eps*sum(abs(w)*samples[s][0][1]
                              for s, w in stencil.items())/precision(amplitude)
            got = residual_edges(precision(newton.BETA)*jump[active]*flux)
            error = np.linalg.norm(got-action)
            cancellation = np.linalg.norm(residual_edges(
                precision(newton.BETA)*abs(jump[active])*bound, bound=True))
            usable = bool(cancellation < budget)
            ok = bool(usable and error <= budget)
            row.update(stencil=stencil, absolute_error=float(error),
                       relative_error=float(error/np.linalg.norm(action)),
                       cancellation=float(cancellation), usable=usable, passed=ok)
            if usable:
                streak = streak+1 if ok else 0
                passed |= streak >= 2
        report["amplitudes"].append(row)
    report["seconds"] = perf_counter()-start
    return passed


if __name__ == "__main__":
    first, second, archive, prior, output = map(Path, sys.argv[1:])
    evidence = json.loads(prior.read_text())
    newton.require(evidence["job_id"] == "2355407", "Unexpected prior audit")
    checkpoints = iter(("200", "217"))

    def audit(iteration, edges, state, dt, base, matrices, rng, report):
        step = next(checkpoints)
        previous = evidence["checkpoints"][step]
        report["reused_audit_job"] = evidence["job_id"]
        for name in (*newton.NAMES, "mixed"):
            check = previous["derivatives"][name]
            newton.require(check["assembly_pass"], "Prior assembly audit did not pass")
            passed = np.array(check["passed"])
            if step == "200" and name == "k":
                newton.require(not passed[1, 3], "Target block was already verified")
                passed[1, 3] = True
            newton.require(passed.all(), "Another prior derivative block remains unresolved")
        if step == "217":
            return previous["verified"]
        return isolated_comparison(iteration, edges, state, dt, base, matrices[1], rng, report)

    newton.audit = audit
    raise SystemExit(0 if newton.run((first, second), archive, output) else 2)
