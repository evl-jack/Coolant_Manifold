"""Regenerate the fixed-pipe derivative tables and Taylor figure from JSON."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


RESPONSE_LABELS = ("Normalized pumping power", "Side-flow fraction")
PLOT_CASES = ("gray", "beta_8")
CASE_LABELS = {"original": "Original binary design", "gray": "Interior gray design",
               "beta_0": "Selected β = 0", "beta_2": "Selected β = 2", "beta_8": "Selected β = 8"}


def write_derivative_table(report, destination):
    columns = ("case", "beta", "direction", "response", "difference", "step", "adjoint", "finite_difference",
               "absolute_error", "relative_error", "scaled_error", "fixed_state_resistance_tangent_error",
               "E0", "E1", "max_forward_residual", "max_mass_imbalance", "directional_test_passed")
    with destination.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(columns)
        for name, case in report["cases"].items():
            for direction, check in case["directions"].items():
                for index, response in enumerate(report["responses"]):
                    for row in check["records"]:
                        writer.writerow((name, case["beta"], direction, response, check["difference"], row["step"],
                                         check["adjoint"][index], row["fd"][index], row["absolute_error"][index],
                                         row["relative_error"][index], row["screen_scaled_error"][index],
                                         row["fixed_state_resistance_tangent_error"], row["e0"][index], row["e1"][index],
                                         check["maximum_final_residual"], check["maximum_mass_imbalance"], check["passed"]))


def plot_taylor(report, destination):
    cases = [name for name in PLOT_CASES if name in report["cases"]]
    if not cases:
        cases = list(report["cases"])[:2]
    figure, axes = plt.subplots(2, len(cases), squeeze=False, sharex=True,
                               figsize=(4.0 * len(cases), 5.4), layout="constrained")
    for col, name in enumerate(cases):
        directions = report["cases"][name]["directions"]
        direction_name = "random" if "random" in directions else next(iter(directions))
        case = directions[direction_name]
        records = case["records"]
        steps = np.array([row["step"] for row in records])
        axes[0, col].set_title(f"{CASE_LABELS[name]}\n{case['difference'].capitalize()} differences")
        for row, response in enumerate(RESPONSE_LABELS):
            axis = axes[row, col]
            e0 = np.array([entry["e0"][row] for entry in records])
            e1 = np.array([entry["e1"][row] for entry in records])
            axis.loglog(steps, np.where(e0 > 0, e0, np.nan), "o-", color="#165a8b", label=r"$E_0$")
            axis.loglog(steps, np.where(e1 > 0, e1, np.nan), "s-", color="#b24f20", label=r"$E_1$")
            reference = steps[:4]
            if e0[0] > 0:
                axis.loglog(reference, 1.6 * e0[0] * reference / reference[0], "--", color="#165a8b",
                            alpha=.6, label="Slope 1")
            if e1[0] > 0:
                axis.loglog(reference, .35 * e1[0] * (reference / reference[0]) ** 2, "--", color="#b24f20",
                            alpha=.6, label="Slope 2")
            axis.axhline(case["nominal_taylor_noise_floor"], color="#666666", linewidth=.8, linestyle=":",
                         label="Taylor screening floor")
            axis.grid(visible=True, which="major", color="#d4dce1", linewidth=.5)
            if col == 0:
                axis.set_ylabel(f"{response}\nAbsolute Taylor remainder")
            if row == 1:
                axis.set_xlabel("Raw-density step size")
    axes[0, -1].legend(loc="best", fontsize=7, frameon=False)
    figure.savefig(destination.with_suffix(".pdf"))
    figure.savefig(destination.with_suffix(".png"), dpi=230)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--additional-reports", type=Path, nargs="*", default=())
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    for path in args.additional_reports:
        additional = json.loads(path.read_text(encoding="utf-8"))
        for key in ("responses", "seed", "steps", "fem", "snapshot_root"):
            if additional[key] != report[key]:
                raise ValueError(f"Mismatched {key} in {path}")
        if report["cases"].keys() & additional["cases"].keys():
            raise ValueError(f"Duplicate design case in {path}")
        report["cases"].update(additional["cases"])
        if "filter_and_transfer" in additional:
            report["filter_and_transfer"] = additional["filter_and_transfer"]
        report["passed"] = report["passed"] and additional["passed"]
    if args.additional_reports:
        (args.output / "adjoint_merged.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n",
                                                        encoding="utf-8")
    write_derivative_table(report, args.output / "adjoint_directional_derivatives.csv")
    plot_taylor(report, args.output / "adjoint_taylor_remainders")


if __name__ == "__main__":
    main()
