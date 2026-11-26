"""Sol-only, independent fixed-mask benchmark for a displacement correction.

Scale the nominal channel half-width H to one, and the velocity by G H²/μ
for an imposed axial pressure gradient G. On the half-domain 0 <= y <= B,
solve -u'' + I[y >= a] u/epsilon² = 1 with u'(0)=0 and u(B)=0.
The forcing covers the whole immersed domain, including the penalized pad.
Here epsilon=sqrt(mu/alpha)/H, and a fixes the step in the computational
penalty mask. The target wall stays at y=1. The independent fitted profile is
(1-y²)/2 with full flux 2/3, in these units.

With outer padding D=B-a, continuity across a gives u'(a)=-a and
u(a)=epsilon²+a*epsilon*tanh(D/epsilon)-epsilon²*sech(D/epsilon).
The extrapolated core parabola vanishes at sqrt(a²+2*u(a)). Thus, for
D >> epsilon and a=1, the mask's apparent no-slip wall is outward by
epsilon+O(epsilon²). Test a=1-epsilon to extend the penalty into the nominal
fluid and cancel this leading tangential displacement. The fitted geometry
and the nominal integration boundary stay at y=1 in both cases.

This binary, parallel-flow diagnostic introduces no changes to the manifold
PDE, adjoint, optimizers or production acceptance criteria. Axial transport
outside the nominal channel is reported; this fixture has no wall-normal
flow and cannot validate manifold connectivity or wall-normal leakage.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def interface_speed(a, outer, epsilon):
    ratio = (outer - a) / epsilon
    decay = np.exp(-ratio)
    sech = 2 * decay / (1 + decay**2)
    return epsilon**2 * (1 - sech) + a * epsilon * np.tanh(ratio)


def _sinh_ratio(numerator, denominator):
    """Evaluate sinh(x)/sinh(z) without overflow for 0 <= x <= z."""
    return (np.exp(numerator - denominator) * -np.expm1(-2 * numerator)
            / -np.expm1(-2 * denominator))


def exact_velocity(y, a, outer, epsilon):
    y = np.asarray(y)
    slip = interface_speed(a, outer, epsilon)
    pad = (outer - a) / epsilon
    t = np.clip((y - a) / epsilon, 0, pad)
    inside = slip + (a**2 - y**2) / 2
    outside = (epsilon**2 + (slip - epsilon**2) * _sinh_ratio(pad - t, pad)
               - epsilon**2 * _sinh_ratio(t, pad))
    return np.where(y <= a, inside, outside)


def exact_flux_to(y, a, outer, epsilon):
    """Full symmetric axial flux across -y <= coordinate <= y."""
    core = min(y, a)
    slip = interface_speed(a, outer, epsilon)
    flux = core * (slip + a**2 / 2) - core**3 / 6
    if y > a:
        z = (outer - a) / epsilon
        t = (y - a) / epsilon
        base = -np.expm1(-2 * z)
        # Stable integrals of the two hyperbolic interpolation functions.
        near = -np.expm1(-t) * -np.expm1(t - 2 * z) / base
        far = np.exp(t - z) * np.expm1(-t)**2 / base
        flux += epsilon**2 * (y - a) + epsilon * ((slip - epsilon**2) * near - epsilon**2 * far)
    return 2 * flux


def _mesh(cuts, maximum_spacing):
    pieces = [np.linspace(a, b, int(np.ceil((b - a) / maximum_spacing)) + 1)[:-1]
              for a, b in zip(cuts[:-1], cuts[1:])]
    return np.concatenate([*pieces, cuts[-1:]])


def _p2_solution(edges, a, epsilon):
    """Independent scalar P2 assembly with a fitted coefficient jump."""
    widths = np.diff(edges)
    elements = len(widths)
    stencil = 2 * np.arange(elements)[:, None] + np.arange(3)
    stiffness = np.array([[7, -8, 1], [-8, 16, -8], [1, -8, 7]]) / 3
    mass = np.array([[4, 2, -1], [2, 16, 2], [-1, 2, 4]]) / 30
    drag = (edges[:-1] >= a) / epsilon**2
    matrices = (stiffness[None, :, :] / widths[:, None, None]
                + (widths * drag)[:, None, None] * mass[None, :, :])
    rows = np.repeat(stencil, 3, axis=1)
    columns = np.tile(stencil, (1, 3))
    size = 2 * elements + 1
    operator = coo_matrix((matrices.ravel(), (rows.ravel(), columns.ravel())), shape=(size, size)).tocsr()
    rhs = np.zeros(size)
    np.add.at(rhs, stencil, widths[:, None] * np.array([1, 4, 1]) / 6)
    # Natural symmetry at the first DOF, essential zero velocity at the last.
    velocity = np.append(spsolve(operator[:-1, :-1], rhs[:-1]), 0.)
    if not np.isfinite(velocity).all():
        raise RuntimeError("The scalar P2 channel solve produced nonfinite values")
    return velocity, stencil


def _discrete_record(edges, a, outer, epsilon, fitted_flux):
    velocity, stencil = _p2_solution(edges, a, epsilon)
    local = velocity[stencil]
    widths = np.diff(edges)
    nominal = edges[1:] <= 1.
    integrals = widths * (local @ np.array([1, 4, 1])) / 3
    flux = float(integrals[nominal].sum())
    solid_flux = float(integrals[~nominal].sum())
    vertices = velocity[::2]
    wall_index = int(np.searchsorted(edges, 1.))
    interface_index = int(np.searchsorted(edges, a))
    exact_flux = exact_flux_to(1., a, outer, epsilon)
    exact_solid_flux = exact_flux_to(outer, a, outer, epsilon) - exact_flux
    q, w = np.polynomial.legendre.leggauss(8)
    xi = (1 + q) / 2
    basis = np.array([(1 - xi) * (1 - 2 * xi), 4 * xi * (1 - xi), xi * (2 * xi - 1)])
    y = edges[:-1, None] + widths[:, None] * xi[None, :]
    numerical = local @ basis
    exact = exact_velocity(y, a, outer, epsilon)
    fitted = (1 - y[nominal] ** 2) / 2
    return {
        "elements_half_channel": int(len(widths)),
        "max_element_width_over_layer": float(widths.max() / epsilon),
        "nominal_wall_velocity": float(vertices[wall_index]),
        "exact_nominal_wall_velocity": float(exact_velocity(1., a, outer, epsilon)),
        "apparent_wall_offset_from_target": float(np.sqrt(a**2 + 2 * vertices[interface_index]) - 1),
        "exact_apparent_wall_offset_from_target": float(np.sqrt(a**2 + 2 * interface_speed(a, outer, epsilon)) - 1),
        "nominal_axial_flux": flux,
        "exact_penalized_nominal_axial_flux": float(exact_flux),
        "nominal_flux_discretization_error_relative_to_fitted": float((flux - exact_flux) / fitted_flux),
        "nominal_flux_modeling_error_relative_to_fitted": float((exact_flux - fitted_flux) / fitted_flux),
        "nominal_flux_total_error_relative_to_fitted": float((flux - fitted_flux) / fitted_flux),
        "nominal_solid_axial_flux_relative_to_fitted": solid_flux / fitted_flux,
        "exact_nominal_solid_axial_flux_relative_to_fitted": float(exact_solid_flux / fitted_flux),
        "velocity_l2_error_vs_exact_penalized": float(np.sqrt(np.sum(widths[:, None] * w * (numerical - exact)**2))),
        "nominal_velocity_l2_error_vs_exact_fitted": float(np.sqrt(np.sum(
            widths[nominal, None] * w * (numerical[nominal] - fitted)**2))),
    }


def run(half_width_mm=2., viscosity=1.002e-3, outer=1.5, penalties=(1e6, 1e7),
        spacing_over_layer=(1., .5, .25)):
    """Return a JSON-ready paired comparison at held nominal geometry."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run numerical execution as a Sol Slurm job")
    parameters = [half_width_mm, viscosity, outer, *penalties, *spacing_over_layer]
    if not np.isfinite(parameters).all() or min(parameters) <= 0 or outer <= 1:
        raise ValueError("Require finite positive channel and penalty parameters and outer > 1")
    if not penalties or not spacing_over_layer:
        raise ValueError("Provide at least one penalty and one mesh resolution")
    fitted_flux = 2 / 3
    records = []
    for alpha in penalties:
        epsilon = float(np.sqrt(viscosity / alpha) / (half_width_mm * 1e-3))
        corrected = 1 - epsilon
        if corrected <= 0:
            raise ValueError("The Brinkman length must be less than the nominal half-width")
        for spacing in spacing_over_layer:
            edges = _mesh([0., corrected, 1., outer], epsilon * spacing)
            records.append({
                "alpha_kg_per_m3_s": alpha,
                "brinkman_length_over_half_width": epsilon,
                "brinkman_length_mm": epsilon * half_width_mm,
                "target_max_element_width_over_layer": spacing,
                "uncorrected": {"mask_start": 1., **_discrete_record(edges, 1., outer, epsilon, fitted_flux)},
                "corrected": {"mask_start": corrected, **_discrete_record(edges, corrected, outer, epsilon, fitted_flux)},
            })
    return {
        "scope": "Fixed binary planar channel displacement benchmark only",
        "sol_job_id": os.environ["SLURM_JOB_ID"],
        "nominal_half_width_mm": half_width_mm,
        "dynamic_viscosity_pa_s": viscosity,
        "outer_boundary_over_half_width": outer,
        "fitted_nominal_full_axial_flux": fitted_flux,
        "coordinates_scaled_by": "nominal half-width H",
        "velocity_scaled_by": "G*H^2/mu",
        "flux_per_unit_depth_scaled_by": "G*H^3/mu",
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", type=Path, help="New JSON result path")
    arguments = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical benchmark as a Sol Slurm job")
    result = run()
    with arguments.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
