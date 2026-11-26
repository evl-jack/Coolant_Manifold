"""Re-evaluate frozen Stokes and Navier–Stokes designs on the common pilot mesh."""

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from shapely.geometry import shape

from ..evaluator import HEIGHT_M, REYNOLDS, VELOCITY_SCALE, VISCOSITY, Stage
from ..fixed_pipe_pilot_model import (
    SIDE_CONSTRAINT_SCALE, FixedPipeEvaluator, build_fixed_pipe, metrics_for,
    original_geometry, reconstruct_geometry,
)
from ..reconstruction import write_reconstruction
from ..restoration import FEASIBILITY_TOLERANCE
from ..stokes_fixed_pipe_model import StokesFixedPipeEvaluator, build_stokes_fixed_pipe


FIELD_NAMES = ("raw", "filtered", "projected", "resistance")
DESIGNS = ("original", "stokes_optimized", "navier_stokes_optimized")
MODELS = ("stokes", "navier_stokes")
CSV_FIELDS = (
    "design", "model", "power_w_per_m", "normalized_power",
    "hydraulic_resistance_pa_s_per_m2", "inlet_flow_m2_per_s", "side_flow_m2_per_s",
    "side_fraction", "side_flow_percent", "main_fraction", "side_requirement",
    "scaled_constraint", "constraint_satisfied", "physical_area_mm2",
    "max_turn_degrees", "mass_imbalance", "nonlinear_residual", "fields_file",
)


def read_result(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    if "completed" in data:
        if not data["completed"]:
            raise ValueError(f"Incomplete archived result at {path}")
        return data["results"]
    return data


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def read_design(path, pipeline, stage):
    with np.load(path, allow_pickle=False) as archive:
        design = archive["design"]
        density = pipeline.evaluate(design, q=stage.q, beta=stage.beta, alpha_max=stage.alpha_max)
        for field in FIELD_NAMES:
            np.testing.assert_allclose(getattr(density, field), archive[field], rtol=1e-12, atol=1e-12)
    return design, density


def assert_archived_mesh(archive, domain):
    with np.load(archive, allow_pickle=False) as mesh:
        np.testing.assert_array_equal(domain.fem.points, mesh["points"])
        np.testing.assert_array_equal(domain.fem.triangles, mesh["triangles"])
        np.testing.assert_array_equal(domain.fem.velocity_dofs, mesh["velocity_dofs"])
        np.testing.assert_array_equal(domain.cell_index, mesh["cell_index"])


def write_fields(path, domain, pipeline, design, density, state):
    fem = domain.fem
    np.savez_compressed(
        path,
        pressure_points_mm=fem.points * (1000 * HEIGHT_M),
        velocity_points_mm=fem.velocity_points * (1000 * HEIGHT_M),
        triangles=fem.triangles,
        velocity_dofs=fem.velocity_dofs,
        pressure_pa=state[2 * fem.n_u:] * (VISCOSITY * VELOCITY_SCALE / HEIGHT_M),
        velocity_m_per_s=state[:2 * fem.n_u].reshape(2, -1).T * VELOCITY_SCALE,
        state=state,
        design=design,
        density_centers_mm=pipeline.grid.centers_mm,
        density_x_mm=pipeline.grid.x_edges_mm,
        density_y_mm=pipeline.grid.y_edges_mm,
        editable_mask=pipeline.free,
        **{name: getattr(density, name) for name in FIELD_NAMES},
    )


def compare(archived_root, stokes_fields, output):
    final_stage = archived_root / "stage_03_beta_8"
    stage = Stage(**read_result(final_stage / "summary.json")["settings"])
    reference = read_result(final_stage / "original" / "metrics.json")
    power_reference = reference["power_w_per_m"]
    side_reference = reference["side_fraction"]

    ns_domain, pipeline = build_fixed_pipe()
    stokes_domain, stokes_pipeline = build_stokes_fixed_pipe()
    for domain in (ns_domain, stokes_domain):
        assert_archived_mesh(archived_root / "mesh.npz", domain)
    np.testing.assert_array_equal(pipeline.free, stokes_pipeline.free)

    sources = {
        "original": final_stage / "original" / "fields.npz",
        "stokes_optimized": stokes_fields,
        "navier_stokes_optimized": final_stage / "selected" / "fields.npz",
    }
    data = {name: read_design(path, pipeline, stage) for name, path in sources.items()}
    np.testing.assert_array_equal(data["original"][0], pipeline.initial)
    geometry = {
        name: original_geometry() if name == "original" else reconstruct_geometry(pipeline, density)
        for name, (_, density) in data.items()
    }
    for design, folder in (("original", "original"), ("navier_stokes_optimized", "selected")):
        archived = shape(read_result(final_stage / folder / "geometry.json")["geometry"])
        if geometry[design].geometry.symmetric_difference(archived).area > 1e-11:
            raise ValueError(f"Replay did not reproduce the archived {design} geometry")

    evaluators = {
        "stokes": StokesFixedPipeEvaluator(stokes_domain, stokes_pipeline, stage, power_reference, side_reference),
        "navier_stokes": FixedPipeEvaluator(ns_domain, pipeline, stage, power_reference, side_reference),
    }
    rows = []
    geometries = {}
    for name, (design, density) in data.items():
        contour = name + ".geometry.json"
        record = write_reconstruction(geometry[name], output / contour)
        geometries[name] = {
            "geometry_file": contour,
            "geometry_sha256": record["geometry_sha256"],
            "physical_area_mm2": record["physical_fluid_area_mm2"],
            "connected_ports": record["connected_ports"],
        }
        for model, evaluator in evaluators.items():
            # The first solve for each model starts cold. Subsequent solves use
            # only a previously converged state from that same model.
            value = evaluator.evaluate(design)
            primal = evaluator.state.values
            fields = f"{name}.{model}.npz"
            write_fields(output / fields, evaluator.domain, evaluator.pipeline, design, density, primal)
            metrics = metrics_for(value, geometry[name], evaluator.pipeline)
            scaled_constraint = float(value.constraints[0])
            rows.append({
                "design": name,
                "model": model,
                **metrics,
                "side_flow_percent": 100 * value.side_fraction,
                "side_flow_m2_per_s": float(metrics["inlet_flow_m2_per_s"] * value.side_fraction),
                "side_requirement": side_reference,
                "side_margin_percentage_points": 100 * (value.side_fraction - side_reference),
                "scaled_constraint": scaled_constraint,
                "constraint_satisfied": scaled_constraint <= FEASIBILITY_TOLERANCE,
                "fields_file": fields,
            })
            geometries[name]["max_turn_degrees"] = metrics["max_turn_degrees"]
            print(f"{name} / {model}: P={1e6 * value.power_w_per_m:.9f} µW/m, "
                  f"side={100 * value.side_fraction:.7f}%, g={scaled_constraint:+.4g}", flush=True)

    for name, folder in (("original", "original"), ("navier_stokes_optimized", "selected")):
        archived = read_result(final_stage / folder / "metrics.json")
        replay = next(row for row in rows if row["design"] == name and row["model"] == "navier_stokes")
        for key in ("power_w_per_m", "side_fraction"):
            np.testing.assert_allclose(replay[key], archived[key], rtol=1e-8, atol=1e-12)

    difference = data["stokes_optimized"][0] - data["navier_stokes_optimized"][0]
    measures = {
        "mean_absolute_raw_density_difference": float(np.abs(difference).mean()),
        "rms_raw_density_difference": float(np.sqrt(np.mean(difference ** 2))),
        "maximum_absolute_raw_density_difference": float(np.abs(difference).max()),
        "symmetric_difference_area_mm2": float(geometry["stokes_optimized"].geometry.symmetric_difference(
            geometry["navier_stokes_optimized"].geometry).area),
        "raw_density_difference_scope": "128 editable cells in the archived density-grid order",
        "geometries": geometries,
    }
    differences_by_evaluator = {}
    for model in MODELS:
        by_design = {row["design"]: row for row in rows if row["model"] == model}
        stokes = by_design["stokes_optimized"]
        ns = by_design["navier_stokes_optimized"]
        differences_by_evaluator[model] = {
            "power_stokes_minus_navier_stokes_w_per_m": stokes["power_w_per_m"] - ns["power_w_per_m"],
            "side_stokes_minus_navier_stokes_percentage_points": 100 * (stokes["side_fraction"] - ns["side_fraction"]),
            "resistance_stokes_minus_navier_stokes_pa_s_per_m2": (
                stokes["hydraulic_resistance_pa_s_per_m2"] - ns["hydraulic_resistance_pa_s_per_m2"]
            ),
        }
    result = {
        "stage": {"q": stage.q, "beta": stage.beta, "alpha_max": stage.alpha_max},
        "mesh": {"target_mm": ns_domain.max_spacing_mm, "dofs": ns_domain.fem.size,
                 "triangles": len(ns_domain.fem.triangles)},
        "reference_power_w_per_m": power_reference,
        "reference_side_fraction": side_reference,
        "side_constraint_scale": SIDE_CONSTRAINT_SCALE,
        "side_constraint_feasibility_tolerance_scaled": FEASIBILITY_TOLERANCE,
        "inlet_velocity_scale_m_per_s": VELOCITY_SCALE,
        "model_inertia_coefficients": {"stokes": 0, "navier_stokes": REYNOLDS},
        "port_power_kinetic_coefficient": REYNOLDS,
        "source_fields_sha256": {name: file_hash(path) for name, path in sources.items()},
        "cases": rows,
        "geometry_comparison": measures,
        "differences_by_evaluator": differences_by_evaluator,
    }
    with (output / "cross_evaluation.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    with (output / "cross_evaluation.csv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archived_pilot", type=Path)
    parser.add_argument("stokes_fields", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    compare(args.archived_pilot, args.stokes_fields, args.output)


if __name__ == "__main__":
    main()
