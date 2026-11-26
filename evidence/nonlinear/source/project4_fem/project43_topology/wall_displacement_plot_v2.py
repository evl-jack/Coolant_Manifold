"""Plot the actual Sol checkpoint and evaluated trials of the single v2 experiment."""

import argparse
import csv
import json
import os
from pathlib import Path

import matplotlib
import numpy as np

from .baffle_free_export import _graphics, _state_archive
from .reconstruction import write_reconstruction
from .wall_displacement_experiment_v2 import STAGE, domain_for, physical_geometry, setup
from .wall_displacement_export_v1 import _cell_image
from .wall_displacement_fem_v1 import WallDisplacementEvaluatorV1
from .wall_displacement_geometry_v1 import reconstruct_clipped
from .wall_displacement_optimize_v1 import sharp_reference
from .wall_displacement_smooth_v2 import GeometricConstraint


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))["results"]


def collect(lane):
    rows = []
    for archive in sorted(lane.glob("restoration_*.json")):
        rows.extend((archive.stem, item["physical_objective"], item["inequalities"], item["design"])
                    for item in read(archive)["trial_history"])
    if (lane / "mma.json").is_file():
        rows.extend(("mma", item["objective"], item["inequalities"], item["design"])
                    for item in read(lane / "mma.json")["evaluation_history"])
    return rows


def run(output, lane, seed, checkpoint, common):
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pipeline, initial, _ = setup(checkpoint)
    reference, _ = sharp_reference(common)
    rows = collect(lane)
    if (lane / "terminal/checkpoint.json").exists():
        design = np.asarray(read(lane / "terminal/checkpoint.json")["design"])
        selection = "mma_terminal"
    elif rows:
        chosen = min(rows, key=lambda row: max(row[2]))
        design = np.asarray(chosen[3])
        selection = f"least_maximum_violation_{chosen[0]}"
    else:
        design = initial
        selection = "verification_seed"
    current = output / selection
    current.mkdir()
    evaluator = WallDisplacementEvaluatorV1(domain_for(pipeline), pipeline, STAGE, reference)
    response = evaluator.evaluate(design)
    density = evaluator.density(design)
    geometry = reconstruct_clipped(pipeline.grid, density.projected, pipeline.free)
    write_reconstruction(geometry, current / "geometry.json")
    smoothness = GeometricConstraint(pipeline, STAGE)
    geometry_metrics = physical_geometry(smoothness, smoothness.sample(design))
    fields = _state_archive(current / "immersed_fields.npz", evaluator.domain.fem, evaluator.state)
    _cell_image(output / "editable_mask.png", pipeline.grid, pipeline.free, mask=True)
    _cell_image(current / "raw_density.png", pipeline.grid, pipeline.raw(design))
    _graphics(output, pipeline.grid, [dict(folder=current, geometry=geometry, density=density.projected, fields=fields)], {})
    for old, new in (("fitted_pressure", "immersed_pressure"), ("fitted_velocity_magnitude", "immersed_velocity_magnitude")):
        path = current / f"{old}.png"
        path.rename(current / f"{new}.png")
    preflight = read(seed)
    values = [("seed", preflight["seed_power_W_per_m"] / reference, preflight["seed_inequalities"], initial), *rows]
    with (output / "history.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["evaluation", "phase", "normalized_power", "side_flow_percent", "g1", "g2", "g3"])
        for i, (phase, p, g, _) in enumerate(values):
            writer.writerow([i, phase, p, 10. - .1 * g[0], *g])
    power = [row[1] for row in values]
    flow = [10. - .1 * row[2][0] for row in values]
    angle = [row[2][2] for row in values]
    for name, label, series, limit in (("objective_history", "Normalized pumping power", power, None),
                                       ("flow_history", "Side-outlet flow (%)", flow, 10.),
                                       ("geometric_constraint_history", "Scaled wall turning constraint", angle, 0.)):
        fig, ax = plt.subplots(figsize=(9, 3.6), dpi=200)
        ax.plot(range(len(series)), series, ".-", color="#13596b")
        if limit is not None:
            ax.axhline(limit, color="#93483b", linestyle="--", linewidth=1)
        ax.set(xlabel="Archived design evaluation", ylabel=label)
        fig.tight_layout()
        fig.savefig(output / f"{name}.png")
        plt.close(fig)
    prior = read(checkpoint)
    (output / "plots.json").write_text(json.dumps({"selection": selection, "history_count": len(rows),
        "side_fraction": float(response.side_fraction), "power_W_per_m": float(response.power_w_per_m),
        "physical_geometry": geometry_metrics,
        "previous_beta4": {"side_fraction": prior["side_fraction"], "power_W_per_m": prior["power_w_per_m"],
                           "max_turn_degrees": prior["wall_geometry_audit"]["max_turn_degrees"],
                           "fluid_area_mm2": prior["wall_geometry_audit"]["fluid_area_mm2"]},
        "flow_field": "Immersed Navier–Stokes at q=0.10 and beta=8"}, indent=2), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("output", "lane", "seed", "checkpoint", "common"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run the replay and plotting on Sol")
    args.output.mkdir(parents=True)
    run(args.output, args.lane, args.seed, args.checkpoint, args.common)


if __name__ == "__main__":
    main()
