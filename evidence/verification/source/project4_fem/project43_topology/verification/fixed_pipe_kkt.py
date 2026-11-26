"""Re-evaluate archived fixed-pipe designs and audit reduced KKT conditions."""

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib
import numpy as np

from ..evaluator import Stage
from ..fixed_pipe_pilot_model import SIDE_CONSTRAINT_SCALE, FixedPipeEvaluator, build_fixed_pipe
from ..gate2 import environment, write_report
from ..kkt_audit import audit_kkt

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402


TOLERANCES = dict(primal_tolerance=1e-6, dual_tolerance=1e-8,
                  stationarity_tolerance=1e-4, complementarity_tolerance=1e-5)


def archived_result(path):
    record = json.loads(path.read_text(encoding="utf-8"))
    if not record["completed"]:
        raise ValueError(f"An archived result is incomplete at {path}")
    return record["results"]


def load_checkpoint(path):
    with np.load(path, allow_pickle=False) as fields:
        return SimpleNamespace(values=fields["state"]), fields["design"]


def audit_point(design, value, evaluator):
    report = audit_kkt(design, value.gradient, value.constraints, value.jacobian, **TOLERANCES)
    g = float(value.constraints[0])
    multiplier = float(report.multipliers[0])
    slack = 100.0 * (value.side_fraction - evaluator.side_reference)
    return {
        **asdict(report),
        "design": design,
        "normalized_power": value.objective,
        "power_w_per_m": value.power_w_per_m,
        "power_reference_w_per_m": evaluator.power_reference,
        "side_fraction": value.side_fraction,
        "reference_side_fraction": evaluator.side_reference,
        "side_slack_percentage_points": slack,
        "scaled_side_constraint": g,
        "side_constraint_activity": (
            "inactive" if g < -TOLERANCES["primal_tolerance"] else
            "violated" if g > TOLERANCES["primal_tolerance"] else "active within tolerance"
        ),
        "marginal_power_w_per_m_per_pp": multiplier * evaluator.power_reference * 0.01 / SIDE_CONSTRAINT_SCALE,
        "mass_imbalance": value.mass_imbalance,
        "nonlinear_residual": value.nonlinear_residual,
        "lower_count": int(np.count_nonzero(design == 0)),
        "upper_count": int(np.count_nonzero(design == 1)),
        "interior_count": int(np.count_nonzero((design > 0) & (design < 1))),
        "objective_gradient": value.gradient,
        "constraint_gradient": value.jacobian[0],
        "lagrangian_gradient": value.gradient + multiplier * value.jacobian[0],
    }


def audit_stage(domain, pipeline, stage_path, replay_history):
    stage_summary = archived_result(stage_path / "summary.json")
    stage = Stage(**stage_summary["settings"])
    original_file = stage_path / "original" / "fields.npz"
    selected_file = stage_path / "selected" / "fields.npz"
    original_state, original = load_checkpoint(original_file)
    selected_state, selected = load_checkpoint(selected_file)
    np.testing.assert_array_equal(original, pipeline.initial)
    baseline_evaluator = FixedPipeEvaluator(domain, pipeline, stage, 1.0, 0.0, initial_state=original_state)
    baseline = baseline_evaluator.evaluate(original)
    evaluator = FixedPipeEvaluator(domain, pipeline, stage, baseline.power_w_per_m, baseline.side_fraction,
                                   initial_state=selected_state)
    endpoint = evaluator.evaluate(selected)
    for label, replayed in (("original", baseline), ("selected", endpoint)):
        expected = stage_summary[label]
        np.testing.assert_allclose(
            [replayed.power_w_per_m, replayed.side_fraction],
            [expected["power_w_per_m"], expected["side_fraction"]],
            rtol=1e-8, atol=1e-13,
        )
    report = audit_point(selected, endpoint, evaluator)
    report.update(settings=asdict(stage), native_status=stage_summary["mma_native_status"],
                  native_message=stage_summary["mma_native_message"], selected_from=stage_summary["selected_from"],
                  original_fields_sha256=sha256(original_file.read_bytes()).hexdigest(),
                  selected_fields_sha256=sha256(selected_file.read_bytes()).hexdigest())

    if replay_history:
        audit_evaluator = FixedPipeEvaluator(domain, pipeline, stage, baseline.power_w_per_m,
                                             baseline.side_fraction, initial_state=original_state)
        optimizer = archived_result(stage_path / "mma.json")
        history = []
        for iteration, record in enumerate(optimizer["evaluation_history"], start=1):
            design = np.asarray(record["design"])
            value = audit_evaluator.evaluate(design)
            np.testing.assert_allclose(
                [value.objective, value.constraints[0]],
                [record["objective"], record["inequalities"][0]], rtol=1e-7, atol=1e-8,
            )
            point = audit_point(design, value, audit_evaluator)
            fields = ("normalized_power", "side_slack_percentage_points", "scaled_side_constraint",
                      "stationarity_inf", "complementarity_inf", "primal_violation", "dual_violation",
                      "multipliers", "verified")
            history.append({"unique_evaluation": iteration, **{field: point[field] for field in fields}})
            print(f"beta={stage.beta:g} evaluation={iteration} rKKT={point['stationarity_inf']:.6e}", flush=True)
        report["mma_history"] = history
    return report


def save_figure(fig, output, name):
    fig.savefig(output / f"{name}.png", dpi=250)
    fig.savefig(output / f"{name}.pdf")
    plt.close(fig)


def plot_density_audit(output, pipeline, endpoint):
    x, y = pipeline.grid.centers_mm[pipeline.free].T
    x_edges = np.r_[np.unique(x) - 0.5 * pipeline.grid.spacing_mm, x.max() + 0.5 * pipeline.grid.spacing_mm]
    y_edges = np.r_[np.unique(y) - 0.5 * pipeline.grid.spacing_mm, y.max() + 0.5 * pipeline.grid.spacing_mm]
    shape = (len(y_edges) - 1, len(x_edges) - 1)
    design = endpoint["design"]
    bounds = np.select([design == 0, design == 1], [0, 2], default=1).reshape(shape)
    residual = np.abs(endpoint["projected_residual"]).reshape(shape)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.0), layout="constrained")
    categories = ListedColormap(["#344c70", "#d8ba64", "#d4e3e3"])
    image = axes[0].pcolormesh(x_edges, y_edges, bounds, cmap=categories,
                              norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5], 3))
    fig.colorbar(image, ax=axes[0], ticks=[0, 1, 2], format=plt.FuncFormatter(
        lambda value, _: ("Lower", "Interior", "Upper")[int(value)]))
    peak = np.unravel_index(np.argmax(residual), residual.shape)
    image = axes[1].pcolormesh(x_edges, y_edges, residual * 1e8, cmap="magma")
    axes[1].plot(x.reshape(shape)[peak], y.reshape(shape)[peak], "o", ms=6, mfc="none", mec="white")
    fig.colorbar(image, ax=axes[1], label=r"$|m_i-\Pi_{[0,1]}(m_i-w_i)|$  ($10^{-8}$)")
    for axis, title in zip(axes, ("Active density bounds", "Projected stationarity by cell")):
        axis.set(title=title, xlabel="x (mm)", ylabel="y (mm)", aspect="equal")
    save_figure(fig, output, "fixed_pipe_kkt_spatial")


def plot_mma_history(output, reports):
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 5.1), layout="constrained")
    fields = ("normalized_power", "side_slack_percentage_points", "stationarity_inf", "complementarity_inf")
    titles = ("Normalized power", "Side-flow slack (percentage points)",
              "Projected stationarity", "Complementarity")
    for axis, field, title in zip(axes.flat, fields, titles):
        for report in reports:
            history = report["mma_history"]
            axis.plot([point["unique_evaluation"] for point in history], [point[field] for point in history],
                      label=f"β = {report['settings']['beta']:g}")
        axis.set(xlabel="Unique MMA evaluation within stage", title=title)
        axis.grid(alpha=0.22)
    axes[0, 1].axhline(0, color="black", linestyle=":", lw=0.8)
    for axis, term in zip(axes[1], ("stationarity", "complementarity")):
        axis.set_yscale("symlog", linthresh=1e-10)
        axis.axhline(TOLERANCES[f"{term}_tolerance"], color="black", linestyle=":", lw=0.8, label="Tolerance")
    axes[1, 0].legend(frameon=False, fontsize=8)
    save_figure(fig, output, "fixed_pipe_kkt_history")


def plot_archived_responses(output, archive_run):
    """Plot archived responses without treating them as iteration-level KKT evidence."""
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.0), layout="constrained")
    for stage in sorted(archive_run.glob("stage_*")):
        beta = archived_result(stage / "summary.json")["settings"]["beta"]
        history = archived_result(stage / "mma.json")["evaluation_history"]
        steps = np.arange(1, len(history) + 1)
        axes[0].plot(steps, [point["objective"] for point in history], label=f"β = {beta:g}")
        axes[1].plot(steps, [-100 * SIDE_CONSTRAINT_SCALE * point["inequalities"][0] for point in history])
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].axhline(0, color="black", linestyle=":", lw=0.8)
    for axis, title in zip(axes, ("Normalized power", "Side-flow slack (percentage points)")):
        axis.set(title=title, xlabel="Unique MMA evaluation within stage")
        axis.grid(alpha=0.22)
    save_figure(fig, output, "fixed_pipe_mma_archived_responses")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive_run", type=Path, help="Directory containing stage folders and mesh.npz")
    parser.add_argument("output", type=Path)
    parser.add_argument("--replay-history", action="store_true", help="Re-evaluate every archived MMA trial")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)

    def run():
        domain, pipeline = build_fixed_pipe()
        with np.load(args.archive_run / "mesh.npz", allow_pickle=False) as mesh:
            np.testing.assert_array_equal(domain.fem.points, mesh["points"])
            np.testing.assert_array_equal(domain.fem.triangles, mesh["triangles"])
            np.testing.assert_array_equal(domain.cell_index, mesh["cell_index"])
        stages = sorted(args.archive_run.glob("stage_*")) if args.replay_history else [args.archive_run / "stage_03_beta_8"]
        reports = [audit_stage(domain, pipeline, stage, args.replay_history) for stage in stages]
        final = reports[-1]
        if final["settings"] != {"q": 0.10, "beta": 8.0, "alpha_max": 1e5}:
            raise ValueError("Expected the final archived pilot stage at q=0.10, beta=8, alpha_max=1e5")
        plot_density_audit(args.output, pipeline, final)
        if args.replay_history:
            plot_mma_history(args.output, reports)
        else:
            plot_archived_responses(args.output, args.archive_run)
        return {"environment": environment(), "mesh_dofs": domain.fem.size, "tolerances": TOLERANCES,
                "side_constraint_scale": SIDE_CONSTRAINT_SCALE, "stages": reports}

    if not write_report(args.output, "fixed_pipe_kkt", run):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
