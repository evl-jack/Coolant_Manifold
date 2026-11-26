"""Write the report summary from the archived Sol evidence."""

import json
from pathlib import Path
import sys

root = Path(sys.argv[1])


def read(name):
    return json.loads((root / name).read_text(encoding="utf-8"))


mesh = read("mesh-evidence/mesh_report.json")
adjoint = read("adjoint-evidence/adjoint_merged.json")
extra = read("adjoint-coordinate-supplement.json")
timings = [
    read(f"timing-2408534/{case}/benchmark.json")
    for case in ("pilot-128", "scale-32", "scale-128", "scale-512")
]
kkt = read("kkt-results/fixed_pipe_kkt.json")["results"]
pilot = timings[0]
measured = pilot["including_common_baseline"]
endpoint = kkt["stages"][-1]
primary_screen = mesh["paired_comparisons"]["fitted-raw"]["fitted_performance_verified"]
checks = [d for c in adjoint["cases"].values() for d in c["directions"].values()]
retry = [d for c in extra["cases"].values() for d in c["directions"].values()]
orders = [w["three_point_order"] for d in retry for w in d["taylor_windows"].values() if w]

lines = [
    "# Fixed-pipe Brinkman pilot verification",
    "",
    "The paired mesh screen " + ("verified" if primary_screen else "did not verify") + " a pumping-power saving for the raw reconstructed wall while maintaining the original fitted side-flow fraction.",
    "",
    "## Discrete adjoint",
    "",
    r"For a converged discrete state satisfying $R(z,m)=0$, solve $R_z^T\lambda_F=F_z$ and evaluate $\nabla_mF=F_m-R_m^T\lambda_F$. The fixed-port responses have zero direct design term. The transpose chain passes through the cell-to-quadrature map, rational Brinkman interpolation, projection and anchored Helmholtz filter. The nonlinear tangent differentiates both velocity factors in convection. See [the derivation](source/project4_fem/project43_topology/verification/ADJOINT_METHOD.md).",
    "",
    f"Finite differences passed for {sum(d['derivative_passed'] for d in checks)} of {len(checks)} tested directions. The initial Taylor order screen passed for {sum(d['taylor_passed'] for d in checks)} directions. The larger-step follow-up passed for {sum(d['passed'] for d in retry)} of {len(retry)} previously flagged directions and retained the original Taylor floor. Measured follow-up orders ran from {min(orders):.4f} to {max(orders):.4f}.",
    "",
    f"The full 128-variable centered gradient used {measured['fd_forward_solve_count']} forward solves. Its normalized-power and side-fraction infinity errors were {pilot['gradient_error']['normalized_power']['absolute_error_linf']:.3e} and {pilot['gradient_error']['side_fraction']['absolute_error_linf']:.3e}.",
    "",
    "| Cost case | Controls | FEM DOFs | Adjoint wall time (s) | FD wall time (s) | Aggregate FD CPU time (s) |",
    "|---|---:|---:|---:|---:|---:|",
]
for d in timings:
    totals = d["including_common_baseline"]
    lines.append(f"| {d['case']} | {d['design_variables']} | {d['fem_dofs']} | {d['complete_adjoint']['wall_seconds']:.3f} | {totals['fd_elapsed_wall_seconds']:.3f} | {totals['fd_aggregate_compute_cpu_seconds']:.3f} |")
lines += [
    "",
    f"The serial pilot speedup was {measured['measured_serial_fd_over_adjoint_speedup']:.2f}× on the same node. The three scaling cases used eight FD workers on one fixed analysis mesh. Each perturbation supplied both responses.",
    "",
    "## KKT audit",
    "",
    f"The archived β = 8 design achieved normalized power {endpoint['normalized_power']:.9f} and scaled side-flow constraint {endpoint['scaled_side_constraint']:.8f}. The side requirement was {endpoint['side_constraint_activity']}, with {endpoint['side_slack_percentage_points']:.6f} percentage points of slack and multiplier {endpoint['multipliers'][0]:.2g}.",
    f"Projected stationarity was {endpoint['stationarity_inf']:.3e}. The primal, dual and complementarity residuals were {endpoint['primal_violation']:.3g}, {endpoint['dual_violation']:.3g} and {endpoint['complementarity_inf']:.3g}. The design had {endpoint['lower_count']} lower-bound, {endpoint['upper_count']} upper-bound and {endpoint['interior_count']} interior variables.",
    "The multiplier gives zero local marginal pumping cost for a sufficiently small increase in the side-flow requirement. The KKT result applies to the immersed β = 8 problem. MMA continuation was unnecessary.",
    f"The history replayed {sum(len(stage['mma_history']) for stage in kkt['stages'])} archived unique MMA evaluations.",
    "",
    "## Mesh and hydraulic comparison",
    "",
    "| Geometry | Target spacing (mm) | Triangles | FEM DOFs | Power (µW/m) | Resistance (Pa·s/m²) | Side flow (%) |",
    "|---|---:|---:|---:|---:|---:|---:|",
]
for name, rows in mesh["series"].items():
    for d in rows:
        if d["state_passed"]:
            lines.append(f"| {name} | {d['mesh_target_mm']:.5g} | {d['triangles']} | {d['dofs']} | {1e6*d['mechanical_power_W_per_m']:.6f} | {d['hydraulic_resistance_pa_s_per_m2']:.3f} | {100*d['side_fraction']:.6f} |")
lines += [
    "",
    "| Comparison at finest paired level | Power saving (µW/m) | Saving (%) | Side gain (pp) | Power uncertainty screen (µW/m) | Raw energy residual sum (µW/m) | Mesh screen |",
    "|---|---:|---:|---:|---:|---:|---|",
]
for name, d in mesh["paired_comparisons"].items():
    a = d["by_mesh"][-1]
    lines.append(f"| {name} | {a['power_saving_W_per_m']*1e6:.6f} | {a['power_saving_percent']:.4f} | {a['side_increase_pp']:.6f} | {d['paired_power_uncertainty_screen_W_per_m']*1e6:.6f} | {a['combined_raw_energy_residual_W_per_m']*1e6:.6f} | {'Pass' if d['screen_passed'] else 'Fail'} |")
lines += [
    "",
    "The power uncertainty screen uses successive paired and individual changes. It is a descriptive discretization estimate. The finite-resistance immersed model has a minimum nominal penetration length of 0.10010 mm and the 0.0625 mm target provides about 1.60 nominal intervals across that length.",
    "",
    "The smoothed wall's coarse-to-medium side-flow change was 0.02505 percentage points against a 0.020 target. Its result remains provisional under the two-comparison criterion. The primary raw fitted geometry met the three-level screen, so the secondary smoothed geometry was left at the planned three levels.",
    "",
    "The fitted areas are 83.000, 85.507 and 85.458 mm² for the original, raw and smoothed geometries. The corresponding maximum wall turns over 0.25 mm are 90.000°, 54.678° and 55.558°. Area and outlet share both change, so the study cannot isolate the effect of turning angle. No compatible matched-area and matched-flow control was available.",
    "",
    "## Reproducing the evidence",
    "",
    "The complete scripts, frozen inputs, machine-readable results, PDF and PNG figures, Slurm logs and recorded job identifiers are included. See the method notes next to the verification scripts. The staging launcher reproduces the main Slurm studies when run in a fresh copy of the staged input and source tree. Run `bash rebuild_figures.sh` to regenerate the plots and tables with the pilot Python dependencies available.",
    "",
    "All 28 common pilot input files matched the supplied Fixed_Pipe_Pilot.zip byte for byte. The archive fingerprints and the diffs for the three modified source files are included. The changes expose independent mesh and density-grid spacing and allow response-only evaluation for finite differences. The original defaults reproduced the archived mesh and responses. The topology and physical smoothed contour were frozen for these studies.",
]
(root / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
