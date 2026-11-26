"""Branch-consistent finite differences and an independent sparse-action check."""

import sys
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.sparse import coo_matrix

from . import project43_sst_exact_newton as newton
from .project43_sst_exact_graph import transport_volume

KINDS = ("original", "graph", "complete")
FLAGS = ("eddy_viscosity", "production", "cross_diffusion_floor",
         "F1_near_viscous", "F1_cross_limit", "F2_near_viscous")
EPS = np.finfo(float).eps


def directional_complex_step(edges, state, direction, base):
    """Perturb all DOFs at once and apply graph winners from real global assembly."""
    m = edges.model
    h = 1e-28
    z = state + 1j*h*direction
    raw = m.residual(z).imag/h
    graph = np.zeros_like(raw)
    volume = transport_volume(m, z[m.dofs])
    boundary = m.edge_local_dofs
    u, v = [z[boundary[:, offset:offset+3]] @ m.edge_velocity.T for offset in (0, 3)]
    normals = m.geometry[-1]
    flux = u*normals[:, 0, None] + v*normals[:, 1, None]
    boundary_values = np.einsum("eq,q,qa,qb->eab", flux, m.edge_w,
                                m.edge_scalar, m.edge_scalar)

    def assemble(nodes, data):
        i, j = np.broadcast_arrays(nodes[:, :, None], nodes[:, None, :])
        return coo_matrix((data.ravel(), (i.ravel(), j.ravel())),
                          shape=(m.nb, m.nb)).tocsr()

    edge_matrix = assemble(m.edges[:, :2], boundary_values)
    ties_report = []
    for j, block in enumerate(m.blocks[2:]):
        i, k = edges.i, edges.j
        matrix = assemble(m.triangles, volume[j]) + edge_matrix
        real = base[1][j][0]
        coefficients = np.stack((np.zeros(len(i), dtype=complex),
                                 np.asarray(matrix[i, k]).ravel(),
                                 np.asarray(matrix[k, i]).ravel()))
        winners, ties = base[3][j]
        selected = coefficients[winners, np.arange(len(i))]
        # The forward derivative of max at a tie takes the largest tangent
        # among the branches attaining the real base value.
        real_entries = np.stack((np.zeros(len(i)), np.asarray(real[i, k]).ravel(),
                                  np.asarray(real[k, i]).ravel()))
        tangents = np.where(real_entries == real_entries.max(axis=0),
                            coefficients.imag/h, -np.inf)
        difference = (tangents.max(axis=0)-selected.imag/h)*(state[block][i]-state[block][k])
        ties_report.append(dict(count=int(np.count_nonzero(ties)),
                                max_forward_minus_selected=float(
                                    newton.BETA*np.max(abs(difference[ties]), initial=0.))))
        contribution = (newton.BETA*selected*(z[block][i]-z[block][k])).imag/h
        np.add.at(graph[block], i, contribution)
        np.add.at(graph[block], k, -contribution)
    return (raw, graph, raw+graph), ties_report


def audit(iteration, edges, state, dt, base, matrices, rng, report):
    start = perf_counter()
    model = iteration.model
    free, scale = iteration.free, iteration.scale
    fixed = np.setdiff1d(np.arange(len(state)), free)
    scalar_ids = np.r_[iteration.ids[2], iteration.ids[3]]
    base_values, _, base_flags, base_branches = base
    report["graph_ties"] = [int(np.count_nonzero(t)) for _, t in base_branches]
    all_pass = True

    def local_terms(z, operators):
        volume = model.element(z[model.dofs], model.geometry)
        boundary = model.edge_flux(z[model.edge_local_dofs], model.geometry[-1])
        graph = []
        for block, operator in zip(model.blocks[2:], operators):
            i, k = edges.i, edges.j
            matrix = operator[0]
            weights = np.maximum(0., np.maximum(np.asarray(matrix[i, k]).ravel(),
                                                 np.asarray(matrix[k, i]).ravel()))
            graph.append(newton.BETA*weights*(z[block][i]-z[block][k]))
        return volume, boundary, *graph

    term_rows = (model.dofs, model.edge_local_dofs[:, 6:], edges.i+model.blocks[2].start,
                 edges.i+model.blocks[3].start)
    opposite_rows = (None, None, edges.j+model.blocks[2].start, edges.j+model.blocks[3].start)
    base_terms = local_terms(state, base[1])

    def local_action(stencil, samples, h):
        original = np.zeros(model.size)
        graph = np.zeros_like(original)
        original_bound = np.zeros_like(original)
        graph_bound = np.zeros_like(original)
        for index, (rows, opposite) in enumerate(zip(term_rows, opposite_rows)):
            value, bound = (original, original_bound) if opposite is None else (graph, graph_bound)
            populated = [(weight, samples[factor][2][index])
                         for factor, weight in stencil.items() if factor]
            center = base_terms[index]
            derivative = sum(w*(term-center) for w, term in populated)/h
            roundoff = 8*EPS*sum(abs(w)*(abs(term)+abs(center))
                                 for w, term in populated)/h
            np.add.at(value, rows, derivative)
            np.add.at(bound, rows, roundoff)
            if opposite is not None:
                np.add.at(value, opposite, -derivative)
                np.add.at(bound, opposite, roundoff)
        return ((original, graph, original+graph),
                (original_bound, graph_bound, original_bound+graph_bound))

    def evaluate(z):
        valid = (np.isfinite(z).all() and np.all(z[scalar_ids] > 0)
                 and np.array_equal(z[fixed], state[fixed]))
        if not valid:
            return None, dict(admissible=False), None
        values, operators, flags, branches = newton.snapshot(iteration, edges, z, dt)
        nut = iteration.fields(z)[7]
        valid = bool(np.isfinite(nut).all() and np.all(nut >= 0)
                     and all(np.isfinite(v[free]).all() for v in values))
        switched = [a != b for a, b in zip(base_flags, flags)]
        graph = [a != b for (a, _), (b, _) in zip(base_branches, branches)]
        return values, dict(admissible=valid,
                            sst={label: int(np.count_nonzero(mask))
                                 for label, mask in zip(FLAGS, switched)},
                            graph=[int(np.count_nonzero(mask)) for mask in graph],
                            first_switched_elements=np.flatnonzero(
                                np.any(np.logical_or.reduce(switched), axis=1))[:8].tolist(),
                            first_switched_edges=[
                                np.c_[edges.i[mask], edges.j[mask]][:8].tolist() for mask in graph]
                            ), local_terms(z, operators) if valid else None

    def usable(info, kind, block):
        if not info["admissible"]:
            return False
        uses_sst = (kind != 1 and block != 1) or (kind != 0 and block >= 2)
        return (not (uses_sst and any(info["sst"].values()))
                and not (kind != 0 and block >= 2 and info["graph"][block-2]))

    for direction_name, ids in zip((*newton.NAMES, "mixed"), (*iteration.ids, free)):
        direction = np.zeros_like(state)
        magnitude = np.ones_like(state)
        magnitude[scalar_ids] = .01*state[scalar_ids]
        direction[ids] = rng.uniform(-1, 1, len(ids))*magnitude[ids]
        actions = [a @ direction[free] for a in matrices]
        complex_actions, ties = directional_complex_step(edges, state, direction, base)
        assembly = {}
        assembly_pass = True
        for kind, action, computed in zip(KINDS, actions, complex_actions):
            d = computed[free]/scale[free]
            error = np.abs(d-action)
            i = int(np.argmax(error))
            assembly[kind] = dict(absolute=float(np.linalg.norm(error)),
                                  reference_norm=float(np.linalg.norm(d)),
                                  max_error_row=int(free[i]))
            assembly_pass &= np.linalg.norm(error) <= 1e-10*max(1., np.linalg.norm(d))
        attempts = []
        streak = np.zeros((3, 4), dtype=int)
        passed = np.zeros((3, 4), dtype=bool)
        amplitudes = newton.AMPLITUDES + ((1e-8, 1e-9)
                                        if direction_name in ("velocity", "mixed") else ())
        for h in amplitudes:
            if h < 1e-7 and passed.all():
                break
            samples = {factor: evaluate(state+factor*h*direction) for factor in (1, -1)}
            record = dict(amplitude=h, errors={})
            local_actions = {}
            for j, (kind, action) in enumerate(zip(KINDS, actions)):
                blocks = []
                for block, block_ids in enumerate(iteration.ids):
                    p = np.searchsorted(free, block_ids)
                    self_norm = float(np.linalg.norm(action[p]))
                    clean = [s for s in (1, -1) if usable(samples[s][1], j, block)]
                    stencil = {1: .5, -1: -.5} if len(clean) == 2 else None
                    if len(clean) == 1:
                        s = clean[0]
                        if 2*s not in samples:
                            samples[2*s] = evaluate(state+2*s*h*direction)
                        if usable(samples[2*s][1], j, block):
                            stencil = {0: -1.5*s, s: 2*s, 2*s: -.5*s}
                    result = dict(action_norm=self_norm, usable=False)
                    if stencil is not None:
                        values = [(factor, weight*(base_values if factor == 0 else
                                   samples[factor][0])[j][block_ids]/scale[block_ids])
                                  for factor, weight in stencil.items()]
                        diff = sum(v for _, v in values)/h
                        absolute = abs(diff-action[p])
                        error = float(np.linalg.norm(absolute))
                        cancellation = float(8*EPS*np.linalg.norm(sum(abs(v) for _, v in values))/h)
                        zero = self_norm <= 1e-12
                        budget = 1e-10 if zero else 1e-4*self_norm
                        if cancellation >= budget or error > budget:
                            key = tuple(stencil)
                            if key not in local_actions:
                                local_actions[key] = local_action(stencil, samples, h)
                            local_value, local_bound = local_actions[key]
                            absolute = abs(local_value[j][block_ids]/scale[block_ids]-action[p])
                            error = float(np.linalg.norm(absolute))
                            cancellation = float(np.linalg.norm(local_bound[j][block_ids]/scale[block_ids]))
                            result["local_differences"] = True
                        dominated = cancellation >= budget
                        relative = None if zero else float(error/self_norm)
                        i = int(np.argmax(absolute))
                        result.update(stencil=list(stencil), absolute=error, relative=relative,
                                      cancellation=cancellation, usable=not dominated,
                                      cancellation_dominated=dominated,
                                      max_error_row=int(block_ids[i]),
                                      max_row_absolute=float(absolute[i]))
                        if not dominated:
                            ok = error <= budget
                            streak[j, block] = streak[j, block]+1 if ok else 0
                            passed[j, block] |= streak[j, block] >= 2
                    blocks.append(result)
                record["errors"][kind] = blocks
            record["sample_checks"] = {str(key): value[1] for key, value in samples.items()}
            attempts.append(record)
        unresolved = []
        for j, block in zip(*np.where(~passed)):
            candidates = [a["errors"][KINDS[j]][block] for a in attempts]
            rows = [item["max_error_row"] for item in candidates if "max_error_row" in item]
            row = rows[-1] if rows else int(iteration.ids[block][
                np.argmax(abs(actions[j][np.searchsorted(free, iteration.ids[block])]))])
            elements = np.flatnonzero(np.any(model.dofs == row, axis=1))
            unresolved.append(dict(matrix=KINDS[j], block=newton.NAMES[block], row=row,
                                   elements=elements.tolist(),
                                   base_limiter_counts={n: int(np.count_nonzero(a[elements]))
                                                        for n, a in zip(FLAGS, base_flags)}))
        report[direction_name] = dict(amplitudes=attempts, passed=passed.tolist(),
                                      complex_step=assembly, selected_vs_forward_ties=ties,
                                      assembly_pass=bool(assembly_pass), unresolved=unresolved)
        all_pass &= bool(passed.all() and assembly_pass)
    report["audit_seconds"] = perf_counter()-start
    return all_pass


if __name__ == "__main__":
    newton.audit = audit
    first, second, archive, output = map(Path, sys.argv[1:])
    raise SystemExit(0 if newton.run((first, second), archive, output) else 2)
