"""Gated feasibility and bounded Newton convergence at archived SST states."""

import hashlib
import json
import os
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu

from . import project43_sst_exact_newton as exact
from .project43_sst_exact_graph import GraphEdges, transport_volume


def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Run:
    def __init__(self, model, initial, free, scale):
        self.it = exact.Iteration(model, free, scale)
        self.edges = GraphEdges(self.it)
        self.m, self.free, self.scale = model, free, scale
        self.fixed = np.setdiff1d(np.arange(len(initial)), free)
        self.prescribed = initial[self.fixed]
        self.scalars = np.r_[self.it.ids[2], self.it.ids[3]]
        self.precision = np.longdouble
        self.geometry = tuple(np.asarray(a, dtype=self.precision) for a in model.geometry)
        i, j = self.edges.i, self.edges.j
        self.keys = np.r_[i*model.nb+j, j*model.nb+i]
        self.order = np.argsort(self.keys)
        self.maps = []
        for nodes in (model.triangles, model.edges[:, :2]):
            row, col = np.broadcast_arrays(nodes[:, :, None], nodes[:, None, :])
            mask = row != col
            index = self.order[np.searchsorted(self.keys[self.order], (row*model.nb+col)[mask])]
            self.maps.append((mask, index))

    def admissible(self, z):
        if not (np.isfinite(z).all() and np.all(z[self.scalars] > 0)
                and np.array_equal(z[self.fixed], self.prescribed)):
            return False
        nut = self.it.fields(z)[7]
        return bool(np.isfinite(nut).all() and np.all(nut >= 0))

    def precise_terms(self, z, geometry=None):
        m = self.m
        geometry = self.geometry if geometry is None else geometry
        c = z[m.dofs]
        boundary = z[m.edge_local_dofs]
        # The affine geometry arrays in transport_volume are promoted by z.
        volume = transport_volume(m, c)
        u, v = [boundary[:, offset:offset+3] @ m.edge_velocity.T for offset in (0, 3)]
        normals = geometry[-1]
        flux = u*normals[:, 0, None]+v*normals[:, 1, None]
        edge_values = np.einsum("eq,q,qa,qb->eab", flux, m.edge_w,
                                m.edge_scalar, m.edge_scalar)
        weights = []
        candidates = []
        for term in volume:
            entries = np.zeros(len(self.keys), dtype=self.precision)
            for values, (mask, index) in zip((term, edge_values), self.maps):
                np.add.at(entries, index, values[mask])
            choices = np.vstack((np.zeros(len(self.edges.i), dtype=self.precision),
                                 entries.reshape(2, -1)))
            candidates.append(choices)
            weights.append(choices.max(axis=0))
        return (m.element(c, geometry), m.edge_flux(boundary, normals),
                weights, candidates, m.fields(c, geometry)[-1])

    def difference(self, z, base, state):
        m = self.m
        result = np.zeros(len(state), dtype=self.precision)
        np.add.at(result, m.dofs, z[0]-base[0])
        np.add.at(result, m.edge_local_dofs[:, 6:], z[1]-base[1])
        i, j = self.edges.i, self.edges.j
        for block, a, b in zip(m.blocks[2:], z[2], base[2]):
            # Caller supplies the scalar jump increment at each sample.
            jump = state[block][i]-state[block][j]
            value = exact.BETA*((a-b)*jump)
            np.add.at(result[block], i, value)
            np.add.at(result[block], j, -value)
        return result

    def slope(self, z, delta, target, matrix, report):
        start = perf_counter()
        f = (target[self.free]/self.scale[self.free]).astype(self.precision)
        analytical = float(f @ (matrix @ delta[self.free]))
        report["predicted_slope"] = analytical
        exact.require(np.isfinite(analytical) and analytical < 0, "Newton action is not descent")
        divisor = max(1., float(np.max(abs(delta[self.scalars]/z[self.scalars]))))
        direction = delta.astype(self.precision)/divisor
        state = z.astype(self.precision)
        base = self.precise_terms(state)
        report["slope_samples"] = []
        cache = {}

        def sample(h):
            if h not in cache:
                state_h = state+self.precision(h)*direction
                if not self.admissible(state_h):
                    cache[h] = None, {"admissible": False}
                    return cache[h]
                terms = self.precise_terms(state_h)
                changes = sum(int(np.count_nonzero(a != b)) for a, b in zip(terms[4], base[4]))
                changes_graph = 0
                ties = 0
                tied_winners = []
                for b, a in zip(base[3], terms[3]):
                    active = b == b.max(axis=0)
                    tied = active.sum(axis=0) > 1
                    ties += int(np.count_nonzero(tied))
                    selected = np.argmax(a, axis=0)
                    tied_winners.extend(selected[tied].tolist())
                    changes_graph += int(np.count_nonzero(
                        ~active[selected, np.arange(len(selected))]))
                r = self.difference(terms, base, state)
                roundoff = np.zeros(len(state), dtype=self.precision)
                np.add.at(roundoff, self.m.dofs, abs(terms[0])+abs(base[0]))
                np.add.at(roundoff, self.m.edge_local_dofs[:, 6:],
                          abs(terms[1])+abs(base[1]))
                i, j = self.edges.i, self.edges.j
                for block, weight, old_weight in zip(self.m.blocks[2:], terms[2], base[2]):
                    flux = exact.BETA*weight*self.precision(h)*(direction[block][i]-direction[block][j])
                    np.add.at(r[block], i, flux)
                    np.add.at(r[block], j, -flux)
                    magnitude = exact.BETA*(abs(weight)+abs(old_weight))*abs(state[block][i]-state[block][j])
                    np.add.at(roundoff[block], i, magnitude)
                    np.add.at(roundoff[block], j, magnitude)
                df = r[self.free]/self.scale[self.free]
                merit_delta = f @ df + .5*(df @ df)
                df_bound = 8*np.finfo(self.precision).eps*roundoff[self.free]/self.scale[self.free]
                merit_bound = abs(f+df) @ df_bound + .5*(df_bound @ df_bound)
                cache[h] = merit_delta, dict(admissible=True, sst_changes=changes,
                                             graph_changes=changes_graph, base_graph_ties=ties,
                                             tied_winners=tied_winners,
                                             merit_roundoff_bound=float(merit_bound))
            return cache[h]

        previous = None
        for h in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10):
            near, near_info = sample(h)
            far, far_info = sample(2*h)
            item = dict(h=h, near=near_info, far=far_info)
            report["slope_samples"].append(item)
            if any(not info["admissible"] or info["sst_changes"] or info["graph_changes"]
                   for info in (near_info, far_info)):
                item["rejected"] = "Inadmissible or switched branch"
                continue
            if near_info["tied_winners"] != far_info["tied_winners"]:
                item["rejected"] = "Directional graph tie winner changed across the stencil"
                continue
            measured = float((2*near-.5*far)*divisor/h)
            bound = (2*near_info["merit_roundoff_bound"] +
                     .5*far_info["merit_roundoff_bound"])*divisor/h
            error = abs(measured-analytical)
            item.update(slope=measured, relative_to_prediction=error/abs(analytical),
                        absolute_slope_roundoff_bound=bound)
            if bound >= 1e-4*abs(analytical):
                item["rejected"] = "Estimated cancellation exceeds the error budget"
                continue
            if measured >= 0 or error+bound > 1e-4*abs(analytical):
                item["rejected"] = "Slope failed descent or agreement with the assembly"
                previous = None
                continue
            if (previous is not None and abs(measured-previous["slope"])+bound+
                    previous["absolute_slope_roundoff_bound"] <= 1e-4*abs(analytical)):
                report["accepted_slope_amplitudes"] = [previous["h"], h]
                report["slope_seconds"] = perf_counter()-start
                return measured
            previous = item
        report["slope_seconds"] = perf_counter()-start
        raise ValueError("No verified branch-consistent forward merit slope")

    def solve(self, z, dt, report):
        t = perf_counter()
        base = exact.snapshot(self.it, self.edges, z, dt)
        matrix = self.m.jacobian(z)+self.edges.jacobian(z, base[1])
        rows = self.free
        scaled = (diags(1/self.scale[rows]) @ matrix[rows][:, rows]).tocsc()
        report["assembly_seconds"] = perf_counter()-t
        t = perf_counter()
        lu = splu(scaled)
        report["factorization_seconds"] = perf_counter()-t
        rhs = -base[0][2][rows]/self.scale[rows]
        t = perf_counter()
        delta = np.zeros_like(z)
        delta[rows] = lu.solve(rhs)
        report["solve_seconds"] = perf_counter()-t
        error = exact.relative_error(rhs, scaled @ delta[rows])
        report["linear_error"] = error
        exact.require(np.isfinite(delta).all() and error <= 1e-8, "Linear solve gate failed")
        report["scalar_corrections"] = {}
        for name, block, ids in zip(exact.NAMES[2:], self.m.blocks[2:], self.it.ids[2:]):
            ratios = delta[ids]/z[ids]
            dof = int(ids[np.argmax(abs(ratios))])
            node = dof-block.start
            report["scalar_corrections"][name] = dict(
                node=node, xy=self.m.points[node].tolist(), scalar=float(z[dof]),
                delta=float(delta[dof]), relative_delta=float(delta[dof]/z[dof]))
        slope = self.slope(z, delta, base[0][2], scaled, report)
        return base, delta, slope

    def metrics(self, z, base):
        raw, _, target = base[0]
        scaled = target[self.free]/self.scale[self.free]
        nut = self.it.fields(z)[7]
        return dict(original=self.it.norms(raw), stabilized=self.it.norms(target),
                    merit=.5*float(scaled @ scaled), skin_friction=float(self.m.cf_at(z)),
                    minima=[float(np.min(z[i])) for i in self.it.ids[2:]],
                    eddy_viscosity=[float(nut.min()), float(nut.max())])

    def search(self, z, dt, base, delta, slope, mapping, report, feasibility=False):
        start = perf_counter()
        initial = self.metrics(z, base)
        ids = self.scalars
        ratios = delta[ids]/z[ids]
        index = int(np.argmax(abs(ratios))) if mapping == "logarithmic" else int(np.argmin(ratios))
        ratio = float(ratios[index])
        if mapping == "logarithmic":
            cap = min(1., 2/abs(ratio)) if ratio else 1.
        else:
            cap = min(1., -.9/ratio) if ratio < 0 else 1.
        dof = int(ids[index])
        scalar = int(dof >= self.m.blocks[3].start)
        node = dof-self.m.blocks[2+scalar].start
        report["limiter"] = dict(scalar=exact.NAMES[2+scalar], node=node,
                                 xy=self.m.points[node].tolist(), relative_delta=ratio)
        report["trials"] = []
        for count in range(16):
            step = cap*.5**count
            if step < 1e-6:
                report["seconds"] = perf_counter()-start
                report["failure"] = "Step length below 1e-6"
                return None
            trial = z+step*delta
            if mapping == "logarithmic":
                trial[ids] = z[ids]*np.exp(step*ratios)
            exact.require(np.isfinite(trial).all() and
                          np.array_equal(trial[self.fixed], self.prescribed),
                          "Nonfinite trial or prescribed values changed")
            entry = dict(step=step, admissible=self.admissible(trial))
            report["trials"].append(entry)
            if not entry["admissible"]:
                continue
            snapshot = exact.snapshot(self.it, self.edges, trial, dt)
            entry.update(self.metrics(trial, snapshot))
            entry["sst_limiter_changes"] = [int(np.count_nonzero(a != b))
                                            for a, b in zip(snapshot[2], base[2])]
            entry["graph_branch_changes"] = [int(np.count_nonzero(a[0] != b[0]))
                                             for a, b in zip(snapshot[3], base[3])]
            entry["armijo"] = entry["merit"] <= initial["merit"]+1e-4*step*slope
            entry["merit_reduction"] = 1-entry["merit"]/initial["merit"]
            entry["max_residual_reduction"] = 1-max(entry["stabilized"])/max(initial["stabilized"])
            if entry["armijo"] and (not feasibility or
                    (entry["merit_reduction"] >= .01 and entry["max_residual_reduction"] >= 0)):
                report["seconds"] = perf_counter()-start
                return trial, entry
        report["seconds"] = perf_counter()-start
        report["failure"] = "No admissible Armijo step met the gate"
        return None


def run(first, prior, output):
    evidence = json.loads(prior.read_text())
    report = dict(job_id=os.getenv("SLURM_JOB_ID"), prior_job=evidence["job_id"],
                  stage_a={"reused_feasibility_job": evidence["job_id"]}, stage_b={}, accepted=[])
    output.mkdir(exist_ok=True)

    def callback(model, initial, free, scale, max_steps):
        solver = Run(model, initial, free, scale)
        sweep = 200
        archived = evidence["checkpoints"]["200"]
        exact.require(archived["verified"] and fingerprint(first) == archived["sha256"],
                      "Verified sweep-200 input hash mismatch")
        exact.require(any(trial.get("promising") for trial in
                          archived["logarithmic"]["trials"]),
                      "Sweep-200 logarithmic feasibility evidence is missing")
        with np.load(first, allow_pickle=False) as cp:
            for key, value in (("free", free), ("scale", scale), ("initial", initial)):
                np.testing.assert_array_equal(cp[key], value)
            exact.require(int(cp["next_step"]) == sweep, "Wrong source sweep")
            z, dt = cp["z"], float(cp["dt"])
        base = exact.snapshot(solver.it, solver.edges, z, dt)
        metrics = solver.metrics(z, base)
        np.testing.assert_allclose(metrics["stabilized"], archived["target"], rtol=.005, atol=2e-14)
        exact.require(solver.admissible(z), "Archived state is inadmissible")
        status = report["stage_b"]
        status.update(start_sweep=sweep, initial=metrics)
        report["accepted"].append(dict(iteration=0, **metrics))
        status["last_verified_iteration"] = 0
        current = z
        try:
            for accepted in range(31):
                last = report["accepted"][-1]
                if max(last["stabilized"]) <= 1e-9:
                    confirmations = [solver.metrics(current,
                        exact.snapshot(solver.it, solver.edges, current.copy(), dt)) for _ in range(2)]
                    exact.require(solver.admissible(current) and all(
                        max(item["stabilized"]) <= 1e-9 for item in confirmations),
                        "Independent tolerance confirmation failed")
                    status.update(passed=True, independent_confirmations=confirmations)
                    return current, [], "Bounded diagnostic complete"
                exact.require(accepted < 30, "Iteration budget exhausted")
                if accepted >= 5:
                    exact.require(last["merit"] <= .99*report["accepted"][-6]["merit"],
                                  "Merit stagnated over five accepted iterations")
                attempt = status["current_attempt"] = {"iteration": accepted+1}
                base, delta, slope = solver.solve(current, dt, attempt)
                info = attempt["line_search"] = {}
                result = solver.search(current, dt, base, delta, slope, "logarithmic", info)
                exact.require(result is not None, info.get("failure", "Line search failed"))
                current, entry = result
                report["accepted"].append(dict(iteration=accepted+1, **entry, diagnostics=attempt))
                status["last_verified_iteration"] = accepted+1
                print(json.dumps(dict(iteration=accepted+1, **entry)), flush=True)
        finally:
            np.savez_compressed(output/"checkpoint.npz", z=current, initial=initial, free=free,
                                scale=scale, dt=dt, source_sweep=int(sweep),
                                accepted=status["last_verified_iteration"])

    try:
        exact.require(evidence["job_id"] == "2355980", "Unexpected verification evidence")
        exact.require(np.finfo(np.longdouble).eps < np.finfo(float).eps,
                      "Extended precision unavailable on Sol")
        exact.segregated.plate.pseudo_transient = callback
        exact.segregated.plate.run()
    except (ValueError, RuntimeError, OSError, AssertionError) as error:
        report["failure"] = str(error)
    (output/"verification.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    print(json.dumps(dict(job_id=report["job_id"], failure=report.get("failure"),
                          accepted=max(0, len(report["accepted"])-1))), flush=True)
    return "failure" not in report


if __name__ == "__main__":
    raise SystemExit(0 if run(*map(Path, sys.argv[1:])) else 2)
