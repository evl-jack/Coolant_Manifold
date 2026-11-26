"""Fixed binary T-junction with Brinkman momentum and fluid-side SST transport.

Run from the package parent with python -m project4_fem.project43_sst.
Scalar wall traces are applied on the known, mesh-aligned interface. Velocity
can cross that interface under a finite Brinkman penalty. This reference
therefore requires a frozen interface.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import block_diag, diags
from scipy.sparse.linalg import spsolve

from .fillet import pressure_weights
from .next_rans_adjoint.sst import FittedRANS, maximum, minimum, blend
from .project43_baseline import Geometry
from .project43_brinkman import create_problem
from .project43_rans_limit import physical_wall_segments


class FrozenSST(FittedRANS):
    def __init__(self, fem, geometry, re, gamma, drag, backflow, lump_destruction=False,
                 wall_floor_fraction=0.):
        if not np.isfinite(wall_floor_fraction) or wall_floor_fraction < 0.:
            raise ValueError("Wall floor fraction must be finite and nonnegative")
        self.wall_floor_fraction = wall_floor_fraction
        super().__init__(fem.points, fem.triangles, re, physical_wall_segments(fem, geometry),
                         outlet_backflow=backflow)
        self.solid = gamma[:, 0] == 0
        self.drag = drag
        self.lump_destruction = lump_destruction

    def prepare(self, points):
        """Regularize SST blending length only, using the local minimum altitude."""
        geometry = super().prepare(points)
        raw = geometry[2]
        vertices = points[self.triangles]
        edges = np.roll(vertices, -1, axis=1) - vertices
        lengths = np.sqrt(np.sum(edges**2, axis=2))
        longest = lengths[np.arange(len(lengths)), np.argmax(lengths.real, axis=1)]
        twice_area = 2 * geometry[3].sum(axis=1)
        h_min = twice_area / longest
        effective = maximum(raw, self.wall_floor_fraction * h_min[:, None])
        # Only the reference geometry supplies persistent reporting quantities.
        if points is self.points:
            self.geometric_distance = raw
            self.local_altitude = h_min
        return (*geometry[:2], effective, *geometry[3:])

    def wall_distance_report(self, height_mm):
        fluid = ~self.solid
        raw = self.geometric_distance[fluid]
        effective = self.geometry[2][fluid]
        clipped = effective > raw
        weights = self.geometry[3][fluid]
        floor = self.wall_floor_fraction * self.local_altitude[fluid]
        return dict(floor_fraction=self.wall_floor_fraction, h_local="minimum triangle altitude",
                    reference_length_mm=height_mm,
                    raw_fluid_min_mm=float(raw.min()*height_mm),
                    effective_fluid_min_mm=float(effective.min()*height_mm),
                    fluid_floor_range_mm=(float(floor.min()*height_mm), float(floor.max()*height_mm)),
                    clipped_fluid_quadrature_points=int(clipped.sum()),
                    total_fluid_quadrature_points=int(clipped.size),
                    clipped_fluid_area_fraction=float(weights[clipped].sum()/weights.sum()))

    def reaction_coefficients(self, coefficients, geometry):
        """SST blending quantities used to split production and destruction."""
        _, gp, distance, _, _, _ = geometry
        kn, wn = coefficients[:, 15:18], coefficients[:, 18:21]
        k = kn @ self.p.T
        omega = wn @ self.p.T
        grad_k = np.einsum("ea,eaj->ej", kn, gp)
        grad_w = np.einsum("ea,eaj->ej", wn, gp)
        cross = (2*.856*np.sum(grad_k*grad_w, axis=1)[:, None])/omega
        near = np.sqrt(k)/(.09*omega*distance)
        viscous = 500*self.nu/(omega*distance**2)
        limit = 4*.856*k/(maximum(cross, 1e-10)*distance**2)
        f1 = np.tanh(minimum(maximum(near, viscous), limit)**4)
        return .09*omega, blend(f1, .075, .0828)*omega, (1-f1)*cross

    def element(self, coefficients, geometry):
        # Remove scalar physics in solid cells while retaining viscous flow.
        local = coefficients.copy()
        local[self.solid, 15:18] = 0
        local[self.solid, 18:21] = 1
        residual = super().element(local, geometry)
        if self.lump_destruction:
            rates = self.reaction_coefficients(local, geometry)
            for block, rate in zip((slice(15,18), slice(18,21)), rates[:2]):
                nodes = local[:, block]
                values = nodes @ self.p.T
                # Row-sum lump the destruction coefficient matrix only:
                # replace integral(phi_i c z) by z_i integral(phi_i c).
                correction = nodes[:, None, :] - values[:, :, None]
                residual[:, block] += np.einsum("eq,qa,eqa->ea", geometry[3]*rate, self.p, correction)
        residual[self.solid, 15:] = 0
        u = np.stack((local[:, :6] @ self.v.T, local[:, 6:12] @ self.v.T), axis=-1)
        drag = np.einsum("eq,qa,eqi->eia", geometry[3]*self.drag, self.v, u)
        residual[:, :12] += drag.reshape(-1, 12)
        return residual


def scalar_diagnostics(model, state, fixed, monitored_vertex):
    """Reassemble scalar rows by term, including upwind exterior traces."""
    coefficients = state[model.dofs]
    geo = model.geometry
    u, _, _, _, _, scalars, grad, _, diffusivity, sources, _ = model.fields(coefficients, geo)
    rates = model.reaction_coefficients(coefficients, geo)
    _, gp, _, weights, diameter, normals = geo
    weights = weights*(~model.solid[:, None])
    along = np.einsum("eqj,eaj->eqa", u, gp)
    ux, uy = (state[model.edge_local_dofs[:, i:i+3]] @ model.edge_velocity.T for i in (0,3))
    edge_flow = ux*normals[:, 0, None] + uy*normals[:, 1, None]
    r = model.residual(state)
    summary = {}
    for i, name in enumerate(("k","omega")):
        nodes = coefficients[:,15+3*i:18+3*i]
        diff = diffusivity[i]
        gradient = grad[i]
        cross = rates[2] if i else np.zeros_like(diff)
        destruction = rates[i]*scalars[i]
        prod = sources[i] + destruction - cross
        terms = dict(diffusion=np.einsum("eq,eaj,ej->ea", weights*diff,gp,gradient[:,0]),
                     convection_volume=-np.einsum("eq,eqa,eq->ea",weights,along,scalars[i]),
                     production=-(weights*prod) @ model.p,
                     destruction=(nodes*((weights*rates[i]) @ model.p) if model.lump_destruction
                                  else (weights*destruction) @ model.p),
                     cross_diffusion=-(weights*cross) @ model.p)
        tau = 1/np.sqrt(4*np.sum(u*u,axis=2)/diameter**2+(4*diff/diameter**2)**2)
        terms["streamline"] = np.einsum("eq,eqa,eq->ea", weights*tau,along,np.sum(u*gradient,axis=2))
        assembled = {}
        for key, local in terms.items():
            row = np.zeros(model.nb)
            np.add.at(row, model.triangles, local)
            assembled[key] = row
        ext = model.edge_bcs[:,i,None]
        z = state[model.edge_local_dofs[:,6+2*i:8+2*i]] @ model.edge_scalar.T
        back = (edge_flow < 0) & np.isfinite(ext)
        for key, trace_flux in (("boundary_transport", np.where(back,0.,edge_flow*z)),
                                ("prescribed_backflow",np.where(back,edge_flow*ext,0.))):
            row = np.zeros(model.nb)
            np.add.at(row, model.edges[:,:2], (trace_flux*model.edge_w) @ model.edge_scalar)
            assembled[key] = row
        actual = r[model.blocks[2+i]]
        err = float(np.max(np.abs(sum(assembled.values())-actual)))
        if err > 1e-10*max(1., np.max(abs(actual))):
            raise RuntimeError(f"{name} decomposed and inherited residuals disagree")
        # Reaction at constrained nodes includes the immersed-wall exchange.
        reaction = actual[fixed].sum()
        boundary = assembled["boundary_transport"].sum()+assembled["prescribed_backflow"].sum()
        volume = np.sum(weights*sources[i])
        defect = abs(boundary-volume-reaction)/max(abs(boundary)+abs(volume)+abs(reaction),1e-14)
        summary[name] = dict(at_monitored_vertex={k: float(v[monitored_vertex]) for k,v in assembled.items()},
                             row_sum=float(actual[monitored_vertex]), max_reassembly_error=err,
                             external_advective_flux=float(boundary), integrated_net_production=float(volume),
                             dirichlet_reaction=float(reaction), relative_balance_defect=float(defect))
        free_vertices = np.setdiff1d(np.arange(model.nb), fixed)
        worst = int(free_vertices[np.argmax(np.abs(actual[free_vertices]))])
        summary[name]["worst_free_vertex"] = dict(
            vertex=worst, global_dof=int(model.blocks[2+i].start+worst), residual=float(actual[worst]),
            signed_terms={key: float(values[worst]) for key, values in assembled.items()})
    return summary


def check_tangent(model, state, free):
    # Probe each active residual block on a perturbed admissible state.
    rng = np.random.default_rng(4304)
    amplitude = np.maximum(np.abs(state), 1e-3)
    probe = state.copy()
    probe[free] += 1e-3*amplitude[free]*rng.uniform(-1, 1, len(free))
    direction = np.zeros(model.size)
    direction[free] = amplitude[free]*rng.uniform(-1, 1, len(free))
    tangent = model.jacobian(probe) @ direction
    step = 1e-6
    centered = (model.residual(probe+step*direction)-model.residual(probe-step*direction))/(2*step)
    errors = []
    for block in model.blocks:
        ids = free[(free >= block.start) & (free < block.stop)]
        denom = np.max(np.abs(tangent[ids]))
        errors.append(float(np.max(np.abs(tangent[ids]-centered[ids]))/max(denom, 1e-14)))
    if not np.all(np.isfinite(errors)) or max(errors) > 2e-5:
        raise RuntimeError(f"Directional tangent check failed {errors}")
    return errors


def refine_newton(model, initial, free, scale, max_steps=10, return_best=False):
    """Permit nonmonotone local Newton steps in a bounded residual neighborhood."""
    state = initial.copy()
    best_state, best_norm = state, np.inf
    ids = free[free >= model.fem.size]
    history = []
    for step in range(max_steps+1):
        residual = model.residual(state)
        norm = float(np.max(np.abs(residual[free]/scale[free])))
        history.append(norm)
        if norm < best_norm:
            best_state, best_norm = state, norm
        if norm <= 1e-11:
            return state, history, None
        if step == max_steps or not np.isfinite(norm) or norm > 1e-8:
            break
        jac = model.jacobian(state)
        update = np.zeros(model.size)
        update[free] = spsolve(jac[free][:,free].tocsc(), -residual[free])
        falling = ids[update[ids] < 0]
        factor = min(1., min(-.9*state[falling]/update[falling], default=1.))
        trial = state+factor*update
        if np.min(trial[ids]) <= 0 or not np.all(np.isfinite(trial)):
            break
        trial_norm = float(np.max(abs(model.residual(trial)[free]/scale[free])))
        if not np.isfinite(trial_norm) or trial_norm > 1e-8:
            break
        state = trial
    return best_state if return_best else state, history, "Local Newton refinement did not meet the residual tolerance"


def correct_omega(model, state, free, scale):
    """Try one omega-block correction against the full coupled steady residual."""
    residual = model.residual(state)
    scaled = np.abs(residual[free]/scale[free])
    norm = float(scaled.max())
    worst = int(free[np.argmax(scaled)])
    block = model.blocks[3]
    trace = dict(entry_norm=norm, entry_worst_dof=worst, trials=[], accepted=False)
    if worst < block.start or norm <= 1e-11:
        return state, trace
    ids = free[free >= block.start]
    jac = model.jacobian(state)
    delta = spsolve(jac[ids][:, ids].tocsc(), -residual[ids])
    if not np.isfinite(delta).all():
        trace["failure"] = "Nonfinite omega-block correction"
        return state, trace
    predicted = jac[:, ids] @ delta
    trace["linear_defect_max_raw"] = float(np.max(abs(predicted[ids]+residual[ids])))
    for halving in range(12):
        factor = 2.**(-halving)
        trial = state.copy()
        trial[ids] += factor*delta
        if not np.isfinite(trial).all() or np.min(trial[free[free >= model.fem.size]]) <= 0:
            trace["trials"].append(dict(factor=factor, admissible=False))
            continue
        trial_residual = model.residual(trial)
        trial_norm = float(np.max(abs(trial_residual[free]/scale[free])))
        if not np.isfinite(trial_norm):
            trace["trials"].append(dict(factor=factor, admissible=True, scaled_steady_trial=None))
            continue
        trace["trials"].append(dict(factor=factor, admissible=True, scaled_steady_trial=trial_norm,
                                    predicted_row_change=float(factor*predicted[worst]),
                                    observed_row_change=float(trial_residual[worst]-residual[worst])))
        if trial_norm < norm:
            trace["accepted"] = True
            return trial, trace
    return state, trace


def pseudo_transient(model, initial, free, scale, max_steps=24):
    """One backward-Euler linearization per artificial time step, in H/U units."""
    best_state = initial.copy()
    fluid_weights = model.geometry[3] * (~model.solid[:, None])
    scalar_mass = np.zeros(model.nb)
    np.add.at(scalar_mass, model.triangles, fluid_weights @ model.p)
    # This velocity mass has identically zero pressure rows and columns, so
    # every correction still satisfies the linearized continuity constraint.
    mass = block_diag((model.fem.drag_matrix(1.), diags(scalar_mass),
                       diags(scalar_mass)), format="csr")[free][:, free]
    positive = free >= model.fem.size
    dt = .01
    history = []
    for _ in range(max_steps):
        steady = model.residual(best_state)[free]
        norm = float(np.max(np.abs(steady/scale[free])))
        if norm <= 1e-11:
            break
        jac = model.jacobian(best_state)[free][:, free]
        update = spsolve((jac + mass/dt).tocsc(), -steady)
        falling = positive & (update < 0)
        factor = min(1., min(-.9*best_state[free[falling]]/update[falling], default=1.))
        trial = best_state.copy()
        trial[free] += factor*update
        valid = np.isfinite(trial).all() and np.min(trial[free[positive]]) > 0
        trial_norm = (float(np.max(abs(model.residual(trial)[free]/scale[free])))
                      if valid else np.inf)
        accepted = np.isfinite(trial_norm) and trial_norm < norm
        history.append(dict(dt=dt, factor=float(factor), scaled_steady_before=norm,
                            scaled_steady_trial=trial_norm if np.isfinite(trial_norm) else None,
                            accepted=bool(accepted)))
        if factor < 1.:
            candidates = np.flatnonzero(falling)
            i = candidates[np.argmin(-best_state[free[candidates]]/update[candidates])]
            block, node = divmod(int(free[i]) - model.fem.size, model.nb)
            history[-1]["positivity_limiter"] = dict(
                field=("k", "omega")[block], node=node,
                xy=model.fem.points[node].tolist(), value=float(best_state[free[i]]),
                proposed_change=float(update[i]))
        if accepted:
            best_state = trial
            dt *= 2.
        else:
            dt *= .25
        if dt < 1e-10:
            break
    final_norm = float(np.max(abs(model.residual(best_state)[free]/scale[free])))
    failure = ("Minimum pseudo-time step reached" if dt < 1e-10 else "Pseudo-transient step budget reached")
    return best_state, history, None if final_norm <= 1e-11 else failure


def run(spacing=.5, lump_destruction=False, refine=False, use_pseudo_time=False, wall_floor_fraction=0.):
    g = Geometry()
    re, rho, mu, alpha = 50., 998.2, .001002, 1e6
    h_m = g.height_mm*1e-3
    bulk = re*mu/(rho*h_m)
    fem, gamma, vel_fixed, vel_values, ports, flux, _ = create_problem(g, spacing)
    resistance = alpha*h_m*h_m/mu*(1-gamma)**2
    intensity, length_mm = .05, .56
    kin = 1.5*intensity**2
    win = np.sqrt(kin)/(.09**.25*(length_mm/g.height_mm))
    backflow = {tuple(e[:2]): (kin, win) for port in ports[1:] for e in port}
    model = FrozenSST(fem, g, re, gamma, resistance/re, backflow, lump_destruction, wall_floor_fraction)
    dist = model.distance(fem.points, fem.points)
    solid_nodes = np.unique(fem.triangles[model.solid])
    physical_wall = np.flatnonzero(dist < 1e-12)
    scalar_wall = np.union1d(solid_nodes, physical_wall)
    inlet = np.setdiff1d(np.unique(ports[0][:, :2]), scalar_wall)
    scalar_fixed = np.union1d(scalar_wall, inlet)
    h = spacing/g.height_mm
    omega_wall = 60*model.nu/(.075*h*h)
    k = kin*np.clip(dist/h, 0, 1)
    omega = win + 6*model.nu/(.075*(dist+h/np.sqrt(10.))**2)
    k[scalar_wall], omega[scalar_wall] = 0, omega_wall
    k[inlet], omega[inlet] = kin, win

    fem.inertia = re
    liquid = fem.solve(vel_fixed, vel_values, fem.load(), resistance=resistance)
    initial = np.r_[liquid.values[:2*fem.n_u], liquid.values[2*fem.n_u:]/re, k, omega]
    fixed = np.r_[vel_fixed, fem.size+scalar_fixed, fem.size+model.nb+scalar_fixed]
    free = np.setdiff1d(np.arange(model.size), fixed)
    forcing = np.zeros(model.size)
    raw = model.residual(initial)
    scale = np.ones(model.size)
    for block in model.blocks:
        ids = free[(free >= block.start) & (free < block.stop)]
        scale[block] = max(1., np.max(np.abs(raw[ids])))
    tangent_errors = check_tangent(model, initial, free)
    state, _, history, failure = model.solve(initial, fixed, initial[fixed], forcing, scale, steps=24)
    refined_history = []
    if refine and lump_destruction and failure and history[-1] < 1e-8:
        state, refined_history, failure = refine_newton(model, state, free, scale, return_best=use_pseudo_time)
    pseudo_history = []
    omega_trace = None
    if use_pseudo_time and failure:
        state, omega_trace = correct_omega(model, state, free, scale)
        state, pseudo_history, failure = pseudo_transient(model, state, free, scale)
    residual = model.residual(state)
    norms = [float(np.max(np.abs(residual[ids]))) for block in model.blocks
             for ids in [free[(free >= block.start) & (free < block.stop)]]]
    report = dict(reynolds_height=re, spacing_mm=spacing, alpha_max_kg_m3_s=alpha,
                  closure_wall_distance=model.wall_distance_report(g.height_mm),
                  geometry="20 x 4 mm main plus centered 1 x 3 mm branch",
                  inlet_intensity=intensity, inlet_length_mm=length_mm, k_in_over_bulk2=kin,
                  omega_in_H_over_bulk=win, wall_omega_H_over_bulk=omega_wall,
                  scalar_domain="Fluid triangles only. Molecular viscosity in the solid",
                  destruction_discretization="row-sum lumped coefficient matrix" if lump_destruction else "consistent Galerkin",
                  scalar_boundary="k=0 and omega=60 nu/(beta1 h^2) on all physical walls including the frozen interface; wall data win at inlet corners",
                  outlets="Zero modified traction; zero diffusive scalar flux; interior outgoing scalar trace and inlet scalar data on reverse flow",
                  active_residual_blocks=["momentum","continuity","k","omega"],
                  free_block_max_residual=norms, block_scales=[float(scale[b.start]) for b in model.blocks],
                  tangent_centered_relative_errors=tangent_errors, newton_scaled_residuals=history,
                  additional_newton_scaled_residuals=refined_history,
                  converged=failure is None, failure=failure, total_triangles=len(fem.triangles),
                  fluid_triangles=int(np.sum(~model.solid)), free_dofs=len(free),
                  minimum_free_scalar={name: float(np.min(state[free[(free >= b.start) & (free < b.stop)]]))
                                       for name, b in zip(("k", "omega"), model.blocks[2:])})
    vertex = np.flatnonzero(np.all(np.isclose(fem.points*g.height_mm, [g.length_mm, spacing]), axis=1)).item()
    if use_pseudo_time:
        report.update(omega_correction=omega_trace, pseudo_transient_steps=pseudo_history,
                      last_pseudo_trial_scaled_steady_residual=(pseudo_history[-1]["scaled_steady_trial"]
                                                                 if pseudo_history else None),
                      best_returned_scaled_steady_residual=float(np.max(abs(residual[free]/scale[free]))),
                      final_scaled_steady_residual=float(np.max(abs(residual[free]/scale[free]))))
    report["monitored_vertex_mm"] = [g.length_mm, spacing]
    report["scalar_diagnostics"] = scalar_diagnostics(model, state, scalar_fixed, vertex)
    for diag in report["scalar_diagnostics"].values():
        row = diag["worst_free_vertex"]
        row["coordinates_mm"] = (fem.points[row["vertex"]]*g.height_mm).tolist()
    worst_dof = int(free[np.argmax(abs(residual[free]/scale[free]))])
    block_index = next(i for i, b in enumerate(model.blocks) if b.start <= worst_dof < b.stop)
    local_index = worst_dof-model.blocks[block_index].start
    coords = fem.velocity_points[local_index % fem.n_u] if block_index == 0 else fem.points[local_index]
    report["worst_free_row"] = dict(block=report["active_residual_blocks"][block_index],
                                    global_dof=worst_dof, coordinates_mm=(coords*g.height_mm).tolist(),
                                    signed_raw=float(residual[worst_dof]),
                                    signed_scaled=float(residual[worst_dof]/scale[worst_dof]))
    if failure:
        jac = model.jacobian(state)
        delta = spsolve(jac[free][:, free].tocsc(), -residual[free])
        constrained = (free >= fem.size) & (delta < 0)
        if np.any(constrained):
            steps = -state[free[constrained]]/delta[constrained]
            idx = np.argmin(steps)
            dof = free[constrained][idx]
            report["positivity_limiter"] = dict(
                scalar="k" if dof < fem.size+model.nb else "omega",
                coordinates_mm=(fem.points[(dof-fem.size) % model.nb]*g.height_mm).tolist(),
                value=float(state[dof]), newton_update=float(delta[constrained][idx]),
                largest_positive_step=float(steps[idx]))
    if failure is None:
        for name in ("k","omega"):
            if (report["minimum_free_scalar"][name] <= 0
                    or report["scalar_diagnostics"][name]["relative_balance_defect"] > 1e-8):
                raise RuntimeError(f"{name} failed positivity or scalar balance")
        static = state[:fem.size].copy()
        static[2*fem.n_u:] -= (2./3.)*state[model.blocks[2]]
        qin, qm, qs = flux @ static
        if min(qs, qm) <= 0 or abs(qin-qs-qm)/qin > 1e-10:
            raise RuntimeError("Port balance failed")
        p_in, p_out = pressure_weights(fem, ports)
        loss = (p_in @ static - static @ (p_out @ static)/(qs+qm))*rho*bulk**2
        nut = model.fields(state[model.dofs], model.geometry)[7][~model.solid]
        report.update(side_share=float(qs/qin), static_pressure_loss_Pa=float(loss),
                      relative_mass_imbalance=float(abs(qin-qs-qm)/qin), max_fluid_nut_over_nu=float(nut.max()/model.nu))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spacing-mm", type=float, choices=(.5, .25), default=.5)
    parser.add_argument("--wall-floor-fraction", type=float, default=0.,
                        help="Lower bound for SST closure distance in local minimum triangle altitudes (0 = original)")
    parser.add_argument("--lump-destruction", action="store_true", help="Row-sum lump both SST destruction terms")
    parser.add_argument("--refine-newton", action="store_true", help="Attempt local Newton refinement on the lumped variant")
    parser.add_argument("--pseudo-transient", action="store_true",
                        help="After Newton refinement, use bounded pseudo-time steps and test the steady residual")
    parser.add_argument("--output", type=Path, help="Fresh JSON report path")
    args = parser.parse_args()
    if args.pseudo_transient and not (args.lump_destruction and args.refine_newton):
        parser.error("--pseudo-transient requires --lump-destruction and --refine-newton")
    report = run(args.spacing_mm, args.lump_destruction, args.refine_newton, args.pseudo_transient,
                 args.wall_floor_fraction)
    result = json.dumps(report, allow_nan=False, indent=2)
    if args.output:
        with args.output.open("x", encoding="utf-8") as out:
            out.write(result+"\n")
    print(result)
    raise SystemExit(0 if report["converged"] else 2)
