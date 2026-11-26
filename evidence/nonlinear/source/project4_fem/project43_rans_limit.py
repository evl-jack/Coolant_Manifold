"""Check the laminar momentum block of the SST FEM on the frozen Brinkman T.

Run from the repository root with python -m project4_fem.project43_rans_limit.
Velocity uses inlet bulk speed and length uses inlet height. The SST element
uses inertial pressure units, while the frozen FEM uses viscous pressure units.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

from .fillet import pressure_weights
from .next_rans_adjoint.sst import FittedRANS
from .project43_baseline import Geometry
from .project43_brinkman import create_problem


def physical_wall_segments(fem, g):
    """Bottom wall plus the four sides of the original immersed interface."""
    left = g.branch_center_mm - g.branch_width_mm / 2
    right = g.branch_center_mm + g.branch_width_mm / 2
    top = g.height_mm + g.branch_length_mm
    corners = [(0, 0), (g.length_mm, 0), (0, g.height_mm), (left, g.height_mm),
               (right, g.height_mm), (g.length_mm, g.height_mm), (left, top), (right, top)]
    nodes = []
    for corner in corners:
        hits, = np.where(np.all(np.isclose(fem.points, np.asarray(corner) / g.height_mm,
                                           rtol=0, atol=1e-12), axis=1))
        if hits.size != 1:
            raise ValueError(f"Expected one vertex at {corner}")
        nodes.append(hits[0])
    return np.array(nodes)[[[0, 1], [2, 3], [4, 5], [3, 6], [4, 7]]]


class LaminarSSTMomentum:
    """Use the SST element's velocity-pressure rows with prescribed zero k.

    A finite auxiliary omega keeps closure denominators defined. The two
    scalar test equations are inactive in laminar mode. Assemble drag directly
    by quadrature in inertial units. Return rows in the frozen FEM convention.
    """

    def __init__(self, fem, geometry, reynolds, drag_inertial):
        self.sst = FittedRANS(fem.points, fem.triangles, reynolds,
                              physical_wall_segments(fem, geometry))
        self.fem = fem
        self.reynolds = reynolds
        self.drag_inertial = drag_inertial

    def residual(self, w):
        fem, sst = self.fem, self.sst
        local_state = np.zeros((len(fem.triangles), 21))
        local_state[:, :15] = w[fem.dofs]
        local_state[:, 12:15] /= self.reynolds
        local_state[:, 18:21] = 1.
        residual = np.zeros(fem.size)
        # The parent element also evaluates inactive scalar rows. Batching
        # bounds temporary field allocations on the baseline mesh.
        for start in range(0, len(local_state), 256):
            idx = slice(start, start + 256)
            geometry = tuple(field[idx] for field in sst.geometry[:-1]) + (sst.geometry[-1],)
            coeff = local_state[idx]
            local = sst.element(coeff, geometry)[:, :15]
            u_nodes = np.stack((coeff[:, :6], coeff[:, 6:12]), axis=-1)
            u = np.einsum("qa,eai->eqi", sst.v, u_nodes)
            drag = np.einsum("eq,qa,eqi->eia", geometry[3] * self.drag_inertial[idx], sst.v, u)
            local[:, :12] += drag.reshape(-1, 12)
            local[:, :12] *= self.reynolds
            np.add.at(residual, fem.dofs[idx], local)
        return residual


def run(reference):
    source = json.loads((reference / "report.json").read_text(encoding="utf-8"))
    g = Geometry(**source["geometry_mm"])
    re = source["reynolds_height"]
    fem, gamma, fixed, values, ports, flux, _ = create_problem(g, source["mesh_max_axis_spacing_mm"])
    viscous_pressure = source["pressure_scale_Pa"]
    speed = source["bulk_inlet_speed_m_s"]
    height = g.height_mm * 1e-3
    with np.load(reference / "field.npz", allow_pickle=False) as field:
        np.testing.assert_array_equal(field["triangles"], fem.triangles)
        np.testing.assert_allclose(field["vertices_mm"], fem.points*g.height_mm, atol=1e-12, rtol=0)
        np.testing.assert_array_equal(field["velocity_dofs"], fem.velocity_dofs)
        np.testing.assert_array_equal(field["fluid_fraction_per_triangle"], gamma[:, 0])
        saved = np.r_[(field["velocity_m_s"] / speed).T.ravel(), field["pressure_Pa"] / viscous_pressure]
    np.testing.assert_allclose(saved[fixed], values, atol=1e-12, rtol=0)
    alpha = source["alpha_max_kg_m3_s"] * (1 - gamma)**2
    rans = LaminarSSTMomentum(fem, g, re, alpha*height / (source["fluid_density_kg_m3"]*speed))
    fem.inertia = re
    operator = fem.viscous_pressure + fem.drag_matrix(alpha*height**2/source["dynamic_viscosity_Pa_s"])
    free = np.setdiff1d(np.arange(fem.size), fixed)
    load = fem.load()
    # Excite velocity in both fluid and solid and perturb pressure. Agreement
    # away from the solution tests the operators as well as their common root.
    state = saved.copy()
    rng = np.random.default_rng(43)
    state[free] += 1.e-3 * rng.uniform(-1, 1, size=len(free))
    reference_residual, _ = fem.residual(state, operator, load, tangent=False)
    probe = rans.residual(state)
    probe_error = float(np.max(np.abs(probe-reference_residual)))
    if probe_error > 1e-10:
        raise RuntimeError(f"Disagreement in off-solution residuals: {probe_error}")

    residuals = []
    for _ in range(5):
        residual = rans.residual(state)
        norm = float(np.max(np.abs(residual[free])))
        residuals.append(norm)
        if norm < 1e-11:
            break
        _, jac = fem.residual(state, operator, load)
        state[free] -= splu(jac[free][:, free].tocsc()).solve(residual[free])
    else:
        raise RuntimeError(f"Closure-off Newton replay failed: {residuals}")
    n_vel = 2*fem.n_u
    velocity_error = float(np.max(np.abs(state[:n_vel]-saved[:n_vel])))
    pressure_error = float(np.max(np.abs(state[n_vel:]-saved[n_vel:])))
    if max(velocity_error, pressure_error) > 1e-8:
        raise RuntimeError("Closure-off solve disagrees with saved dimensionless field")
    inlet, main, side = flux @ state
    imbalance = abs(inlet-main-side)/inlet
    if min(main, side) <= 0 or imbalance > 1.e-10:
        raise RuntimeError("Closure-off solve failed port flux checks")
    p_in, p_out = pressure_weights(fem, ports)
    loss = viscous_pressure * (p_in @ state - state @ (p_out @ state)/(main+side))
    return dict(case="Frozen T-junction SST momentum, prescribed zero eddy viscosity",
                reference=str(reference), reynolds_height=re, spacing_mm=source["mesh_max_axis_spacing_mm"],
                alpha_max_kg_m3_s=source["alpha_max_kg_m3_s"], turbulence_transport_active=False,
                dimensionless_k=0., auxiliary_dimensionless_omega=1., eddy_viscosity=0.,
                off_solution_max_row_discrepancy=probe_error, closure_off_newton_residuals=residuals,
                max_dimensionless_velocity_difference=velocity_error,
                max_viscous_scaled_pressure_difference=pressure_error,
                side_share=float(side/inlet), static_pressure_loss_Pa=float(loss),
                relative_mass_imbalance=float(imbalance))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path,
                        default=Path("output/project43_frozen_t_re50_h025_alpha1e6"))
    parser.add_argument("--output", type=Path, help="New JSON report file")
    args = parser.parse_args()
    result = json.dumps(run(args.reference), indent=2, allow_nan=False)
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(result + "\n")
    print(result)
