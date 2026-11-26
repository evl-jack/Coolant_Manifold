"""Bounded defect correction for the archived flat-plate SST problem."""

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from . import segregated

CHECKPOINT_HASH = "dcfe6d3520ddd2ca637cd2fa5c1ee1baef1cfba5873deed327d09ff9b9106a30"
BETA = .8
TOL = 1e-9


def require(condition, message):
    if not condition:
        raise ValueError(message)


class Iteration:
    def __init__(self, model, free, scale):
        self.model, self.free, self.scale = model, free, scale
        shape = (len(model.triangles), 3, 3)
        self.rows = np.broadcast_to(model.triangles[:, :, None], shape).ravel()
        self.cols = np.broadcast_to(model.triangles[:, None, :], shape).ravel()
        self.ids = [free[(free >= b.start) & (free < b.stop)] for b in model.blocks]
        self.flow = free[free < model.fem.size]
        self.mass = model.fem.drag_matrix(1.)[self.flow][:, self.flow]
        self.timing = dict(assembly_seconds=0., solve_seconds=0.)

    def operator(self, z, j, dt):
        return segregated.scalar_operator(self.model, z, j, self.model.blocks[j+2],
                                          self.free, dt, self.rows, self.cols)

    def residual(self, z, dt):
        raw = self.model.residual(z)
        target = raw.copy()
        for j in range(2):
            _, graph, ids, local = self.operator(z, j, dt)
            target[ids] += BETA * (graph @ z[self.model.blocks[j+2]])[local]
        return raw, target

    def norms(self, residual):
        require(np.isfinite(residual[self.free]).all(), "Nonfinite steady residual")
        return [float(np.max(abs(residual[i] / self.scale[i]))) for i in self.ids]

    def fields(self, z):
        return self.model.fields(z[self.model.dofs], self.model.geometry)

    def scalar(self, z, j, dt):
        start = perf_counter()
        m = self.model
        matrix, graph, ids, local = self.operator(z, j, dt)
        q = z[ids]
        require(np.isfinite(q).all() and np.all(q > 0), "Invalid free scalar")
        target = m.residual(z)[ids] + BETA * (graph @ z[m.blocks[j+2]])[local]
        system = (matrix + graph)[local][:, local].tocsr()
        if j == 1:
            fields = self.fields(z)
            k, omega = fields[5]
            grad_k, grad_omega = fields[6]
            distance, weights = m.geometry[2:4]
            cross = 2 * .856 * np.sum(grad_k * grad_omega, axis=-1) / omega
            near = np.sqrt(k) / (.09 * omega * distance)
            viscous = 500 * m.nu / (distance**2 * omega)
            cap = 4 * .856 * k / (np.maximum(cross, 1e-10) * distance**2)
            blend = np.tanh(np.minimum(np.maximum(near, viscous), cap)**4)
            element = (weights * np.maximum(-(1-blend)*cross, 0)) @ m.p
            sink = np.zeros(m.nb)
            np.add.at(sink, m.triangles, element)
            system += diags(sink[local] / q)
        absolute_rhs = system @ q - target
        correction = np.maximum(0., -absolute_rhs) / q
        system += diags(correction)
        diagonal = system.diagonal()
        off = system.copy()
        off.setdiag(0)
        radius = np.asarray(abs(off).sum(axis=1)).ravel()
        roundoff = 32 * np.finfo(float).eps * (abs(diagonal) + radius)
        entries = off.tocoo()
        positive = int(np.count_nonzero(entries.data > roundoff[entries.row]))
        margin = diagonal - radius
        require(np.isfinite(system.data).all() and np.isfinite(target).all(),
                "Nonfinite scalar coefficients")
        require(np.all(diagonal > 0) and positive == 0 and np.all(margin > roundoff),
                f"Scalar {j} matrix gate failed: off-diagonals={positive}, min margin={margin.min()}")
        self.timing["assembly_seconds"] += perf_counter() - start
        start = perf_counter()
        delta = spsolve(system.tocsc(), -target)
        self.timing["solve_seconds"] += perf_counter() - start
        require(np.isfinite(delta).all(), "Nonfinite scalar update")
        error = float(np.max(abs(system @ delta + target) / self.scale[ids]))
        require(error < 1e-10, f"Scalar {j} linear residual {error}")
        falling = delta < 0
        factor = min(.7, min(-.9*q[falling]/delta[falling], default=1.))
        candidate = q + factor*delta
        require(np.isfinite(candidate).all() and np.all(candidate > 0), "Invalid scalar candidate")
        return ids, candidate, dict(relaxation=factor, linear_residual=error,
                                    minimum=float(candidate.min()),
                                    full_step_minimum=float((q+delta).min()),
                                    positivity_nodes=int(np.count_nonzero(correction)),
                                    positive_offdiagonals=positive,
                                    diagonal_margin_min=float(margin.min()))

    def advance_flow(self, z, dt):
        start = perf_counter()
        flow = self.flow
        rhs = -self.model.residual(z)[flow]
        system = self.model.jacobian(z)[flow][:, flow] + self.mass/dt
        require(np.isfinite(system.data).all() and np.isfinite(rhs).all(),
                "Nonfinite flow coefficients")
        self.timing["assembly_seconds"] += perf_counter() - start
        start = perf_counter()
        delta = spsolve(system.tocsc(), rhs)
        self.timing["solve_seconds"] += perf_counter() - start
        require(np.isfinite(delta).all(), "Nonfinite flow correction")
        z[flow] += .5*delta


def run(checkpoint, output):
    from .test_project43_sst_bounded import gate_a

    output.mkdir(parents=True, exist_ok=True)
    report = dict(gate_a="failed", gate_b="not_run", beta_target=BETA, beta_iteration=1.)
    history = []

    def verify(model, initial, free, scale, max_steps):
        with np.load(checkpoint, allow_pickle=False) as cp:
            z, dt = cp["z"].copy(), float(cp["dt"])
            for key, expected in (("initial", initial), ("free", free), ("scale", scale)):
                np.testing.assert_array_equal(cp[key], expected)
            require(int(cp["next_step"]) == 200, "Expected sweep-200 checkpoint")
        iteration = Iteration(model, free, scale)
        fixed = np.setdiff1d(np.arange(len(z)), free)
        completed = 200
        try:
            report["fixed_state"] = {}
            gate_a(iteration, z, dt, report["fixed_state"])
            report["gate_a"], report["gate_b"] = "passed", "failed"
            small = [0, 0]
            held = 0
            for sweep in range(201, 321):
                trial = z.copy()
                iteration.advance_flow(trial, dt)
                diagnostics = []
                for j in range(2):
                    ids, candidate, audit = iteration.scalar(trial, j, dt)
                    trial[ids] = candidate
                    diagnostics.append(audit)
                    small[j] = small[j]+1 if audit["relaxation"] < 1e-4 else 0
                nut = iteration.fields(trial)[7]
                require(np.isfinite(trial).all() and np.isfinite(nut).all() and np.all(nut >= 0),
                        "Invalid state or sampled eddy viscosity")
                require(np.array_equal(trial[fixed], initial[fixed]), "Prescribed trace changed")
                raw, target = iteration.residual(trial, dt)
                original, stabilized = iteration.norms(raw), iteration.norms(target)
                cf = float(model.cf_at(trial))
                require(np.isfinite(cf), "Nonfinite skin friction")
                z, completed = trial, sweep
                dt = min(1., 1.05*dt)
                record = dict(sweep=sweep, original=original, target=stabilized, cf=cf,
                              nut_min=float(nut.min()), nut_max=float(nut.max()), scalars=diagnostics)
                history.append(record)
                print(json.dumps(record, allow_nan=False), flush=True)
                maximum = max(stabilized)
                held = held+1 if maximum <= TOL else 0
                if held >= 3:
                    cf_window = [item["cf"] for item in history[-3:]]
                    if np.ptp(cf_window) < 1e-4*abs(cf_window[0]):
                        report["gate_b"] = "passed"
                        break
                require(max(small) < 5, "Scalar relaxation below 1e-4 for five sweeps")
                if len(history) > 15 and maximum > TOL:
                    require(maximum <= .99*max(history[-16]["target"]),
                            "Target residual decreased by less than 1% over 15 sweeps")
            if report["gate_b"] != "passed":
                report["failure"] = "Sweep budget exhausted"
        except (ValueError, AssertionError, RuntimeError) as exc:
            report["failure"] = str(exc)
        finally:
            report["timing"] = iteration.timing
            np.savez_compressed(output/"checkpoint.npz", z=z, dt=dt, next_step=completed,
                                initial=initial, free=free, scale=scale)
        return z, [], report.get("failure")

    try:
        digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
        report["input_sha256"] = digest
        require(digest == CHECKPOINT_HASH, "Checkpoint SHA-256 mismatch")
        segregated.plate.pseudo_transient = verify
        segregated.plate.run()
    except (OSError, ValueError, AssertionError, RuntimeError) as exc:
        report["failure"] = str(exc)
    report["history"] = history
    (output/"Project43_SST_bounded_verification.json").write_text(
        json.dumps(report, indent=2, allow_nan=False))
    return report["gate_b"] == "passed"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    raise SystemExit(0 if run(args.checkpoint, args.output) else 2)
