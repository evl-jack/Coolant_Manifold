"""Analytical fixtures for the independent KKT multiplier audit on Sol."""

from dataclasses import asdict, dataclass
import json

import numpy as np

from .kkt_audit import audit_kkt


PRIMAL_TOLERANCE = 1e-6
DUAL_TOLERANCE = 1e-8
STATIONARITY_TOLERANCE = 1e-4
COMPLEMENTARITY_TOLERANCE = 1e-5
WITNESS_TOLERANCE = 1e-8


@dataclass(frozen=True)
class Fixture:
    name: str
    design: tuple
    objective_gradient: tuple
    inequalities: tuple
    constraint_jacobian: tuple
    expected_verified: bool
    rationale: str
    exact_multipliers: tuple | None = None


def _fixtures():
    t = STATIONARITY_TOLERANCE
    return (
        Fixture("lower_box_active", (0.,), (2.,), (), (), True,
                "Positive gradient at the lower bound satisfies stationarity.", ()),
        Fixture("upper_box_active", (1.,), (-2.,), (), (), True,
                "Negative gradient at the upper bound satisfies stationarity.", ()),
        Fixture("lower_box_wrong_sign", (0.,), (-2.,), (), (), False,
                "Negative gradient at the lower bound admits an inward descent direction."),
        Fixture("upper_box_wrong_sign", (1.,), (2.,), (), (), False,
                "Positive gradient at the upper bound admits an inward descent direction."),
        Fixture("lower_box_near_inside_tolerance", (0.75 * t,), (2.,), (), (), True,
                "Projection hits zero and the genuine residual is 0.75 times the tolerance.", ()),
        Fixture("upper_box_near_inside_tolerance", (1. - 0.75 * t,), (-2.,), (), (), True,
                "Projection hits one and the genuine residual is 0.75 times the tolerance.", ()),
        Fixture("lower_box_near_outside_tolerance", (1.25 * t,), (2.,), (0.,), ((1.,),), False,
                "Projection hits zero and every nonnegative multiplier gives residual 1.25 times tolerance."),
        Fixture("upper_box_near_outside_tolerance", (1. - 1.25 * t,), (-2.,), (0.,), ((-1.,),), False,
                "Projection hits one and every nonnegative multiplier gives residual 1.25 times tolerance."),
        Fixture("nearly_lower_box_needs_opposite_sign_bound", (0.75 * t,), (-2.,), (0.,), ((1.,),), True,
                "A positive multiplier is required to prevent an inward step; lambda=2 makes the gradient zero."),
        Fixture("nearly_upper_box_needs_opposite_sign_bound", (1. - 0.75 * t,), (2.,), (0.,), ((-1.,),), True,
                "A positive multiplier is required to prevent an inward step; lambda=2 makes the gradient zero."),
        Fixture("interior_active_coupled", (0.3, 0., 1., 0.6), (-3., 1., -4., 5.), (0., 0., -0.2),
                ((1., 2., 1., -1.), (1., -1., -2., -3.), (0., 0., 0., 0.)), True,
                "Two interior equations uniquely give lambda=(2,1) on the active constraints and preserve both bound signs."),
        Fixture("interior_wrong_dual_sign", (0.5,), (1.,), (0.,), ((1.,),), False,
                "The exact interior gradient balance needs lambda=-1."),
        Fixture("complementarity_cap_allows_witness", (0.5,), (-2.,), (-4e-6,), ((1.,),), True,
                "The stationarity-margin optimizer can use lambda=2, with complementarity residual 8e-6.", (2.,)),
        Fixture("complementarity_cap_blocks_witness", (0.5,), (-2.,), (-6e-6,), ((1.,),), False,
                "Complementarity caps lambda at 5/3; every admissible witness leaves projected stationarity at least 1/3."),
        Fixture("constraint_violated", (0.5,), (0.,), (2e-6,), ((1.,),), False,
                "The inequality violation exceeds primal tolerance even though lambda=0 certifies stationarity."),
        Fixture("rounded_positive_constraint_has_complementarity_cap", (0.5,), (-20.,), (0.75e-6,), ((1.,),), False,
                "Primal tolerance admits the positive constraint, but complementarity caps lambda at 40/3."),
        Fixture("box_violated", (-2e-6,), (1.,), (), (), False,
                "The actual design is outside the box by more than primal tolerance."),
        Fixture("gradient_witness_must_satisfy_all_rows", (0.5, 0.5), (-1., -1.01), (0.,), ((1., 1.),), False,
                "The sole multiplier cannot simultaneously satisfy two interior stationarity intervals separated by 0.01."),
    )


def _witness_residuals(design, objective_gradient, inequalities, jacobian, multipliers):
    # Construct the box projection directly from the trial point, independently
    # of the LP's inequalities, its objective, and its internal margin.
    gradient = objective_gradient + multipliers @ jacobian
    residual = []
    for coordinate, partial in zip(design, gradient):
        trial = coordinate - partial
        if trial < 0.:
            residual.append(float(coordinate))
        elif trial > 1.:
            residual.append(float(coordinate - 1.))
        else:
            residual.append(float(partial))
    inequality = max([0., *inequalities])
    box = max([0., *(-design), *(design - 1.)])
    return {"gradient_lagrangian": gradient.tolist(), "projected_residual": residual,
            "primal": float(max(inequality, box)), "dual": float(max([0., *(-multipliers)])),
            "stationarity": max(abs(value) for value in residual),
            "complementarity": float(max([0., *np.abs(multipliers * inequalities)]))}


def _check_fixture(fixture):
    design = np.asarray(fixture.design)
    gradient = np.asarray(fixture.objective_gradient)
    g = np.asarray(fixture.inequalities)
    jacobian = np.asarray(fixture.constraint_jacobian).reshape(g.size, design.size)
    report = audit_kkt(design, gradient, g, jacobian)
    witness = _witness_residuals(design, gradient, g, jacobian, report.multipliers)
    witness_verified = (witness["primal"] <= PRIMAL_TOLERANCE and witness["dual"] <= DUAL_TOLERANCE
                        and witness["stationarity"] <= STATIONARITY_TOLERANCE
                        and witness["complementarity"] <= COMPLEMENTARITY_TOLERANCE)
    residual_error = max(abs(report.primal_violation - witness["primal"]),
                         abs(report.dual_violation - witness["dual"]),
                         abs(report.stationarity_inf - witness["stationarity"]),
                         abs(report.complementarity_inf - witness["complementarity"]),
                         float(np.max(np.abs(report.projected_residual - witness["projected_residual"]))))
    exact_multiplier_error = None
    if fixture.exact_multipliers is not None:
        exact_multiplier_error = float(np.max(np.abs(report.multipliers - fixture.exact_multipliers), initial=0.))
    active_multiplier_error = None
    if fixture.name == "interior_active_coupled":
        active_multiplier_error = float(np.max(np.abs(report.multipliers[:2] - [2., 1.])))
    passed = (report.verified == witness_verified == fixture.expected_verified and residual_error <= WITNESS_TOLERANCE
              and (exact_multiplier_error is None or exact_multiplier_error <= WITNESS_TOLERANCE)
              and (active_multiplier_error is None or active_multiplier_error <= WITNESS_TOLERANCE))
    return {"fixture": asdict(fixture), "audit_verified": bool(report.verified),
            "lp_status": report.lp_status, "multiplier_source": report.multiplier_source,
            "multipliers": report.multipliers.tolist(), "independent_witness": witness,
            "independent_witness_verified": bool(witness_verified), "residual_error": residual_error,
            "exact_multiplier_error": exact_multiplier_error, "active_multiplier_error": active_multiplier_error,
            "passed": bool(passed)}


def run_kkt_fixtures():
    records = [_check_fixture(fixture) for fixture in _fixtures()]
    return {"passed": all(record["passed"] for record in records), "fixtures": records,
            "limits": {"primal": PRIMAL_TOLERANCE, "dual": DUAL_TOLERANCE, "stationarity": STATIONARITY_TOLERANCE,
                       "complementarity": COMPLEMENTARITY_TOLERANCE, "witness_agreement": WITNESS_TOLERANCE}}


if __name__ == "__main__":
    results = run_kkt_fixtures()
    print(json.dumps(results, indent=2, allow_nan=False))
    raise SystemExit(0 if results["passed"] else 1)
