"""Run the in-house 3D RANS solver or its square-duct verification case."""

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from .mesh import duct, from_polymesh
from .solver import Settings, solve


def poiseuille(points, width, mean):
    odd = np.arange(1, 80, 2)
    a, b = odd[:, None], odd[None, :]
    coefficient = 16 / (np.pi**4 * a*b * (a*a+b*b))
    average = np.sum(coefficient * (4/(np.pi**2 * a*b)))
    modes = np.sin(np.pi*points[:, 1:, None]*odd/width)
    profile = np.einsum("pi,ij,pj->p", modes[:, 0], coefficient, modes[:, 1])
    return mean*profile/average, average


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, help="Path to the validated Sol polyMesh directory")
    parser.add_argument("--laminar", action="store_true")
    parser.add_argument("--nx", type=int, default=12)
    parser.add_argument("--ny", type=int, default=12)
    parser.add_argument("--speed", type=float, help="Duct mean or manifold inlet speed in m/s")
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--order", type=int, choices=(1, 2), default=1)
    parser.add_argument("--cfl", type=float, default=3.)
    parser.add_argument("--tolerance", type=float, default=1e-5)
    parser.add_argument("--intensity", type=float, default=.05, help="Inlet turbulence intensity as a fraction")
    parser.add_argument("--length-scale-mm", type=float, default=.56)
    parser.add_argument("--restart", type=Path, help="Read fields.npz from a run on the same mesh")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (min(args.nx, args.ny) < 2 or min(args.cfl, args.tolerance, args.intensity, args.length_scale_mm) <= 0
            or (args.speed is not None and args.speed <= 0)):
        parser.error("Grid counts must be at least two and all scale parameters must be positive")
    settings = Settings(order=args.order, pseudo_cfl=args.cfl, tolerance=args.tolerance, turbulence=not args.laminar,
                        intensity=args.intensity, length_scale=args.length_scale_mm*.001)
    width = .0035
    mesh = from_polymesh(args.mesh) if args.mesh else duct(args.nx, args.ny, width=width)
    speed = args.speed if args.speed is not None else (.05 if args.laminar and not args.mesh else 2.)
    inlet = np.array([speed, 0., 0.])
    exact = args.mesh is None and args.laminar
    if exact:
        profile, average = poiseuille(mesh.faces[mesh.patches == "inlet"], width, speed)
        inlet = profile[:, None] * np.array([1., 0., 0.])
    print(f"Solving {len(mesh.volumes)} cells with convection order {args.order}", flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "history.csv").open("w", newline="") as file:
        history = csv.DictWriter(file, fieldnames=("iteration", "field_change", "continuity", "mass_imbalance",
            "pressure_drop_pa", "outlet_flow_drift", "pressure_drop_drift", "inlet_m3_s", "main_m3_s",
            "side_m3_s", "main_reverse_m3_s", "side_reverse_m3_s"))
        history.writeheader()

        def monitor(row):
            flows, reverse = row["flow_m3_s"], row["reverse_flow_m3_s"]
            history.writerow({key: row[key] for key in history.fieldnames[:7]} | dict(
                inlet_m3_s=-flows["inlet"], main_m3_s=flows.get("main_outlet", 0.),
                side_m3_s=flows.get("side_outlet", 0.), main_reverse_m3_s=reverse.get("main_outlet", 0.),
                side_reverse_m3_s=reverse.get("side_outlet", 0.)))
            if row["iteration"] == 1 or row["iteration"] % 25 == 0:
                print("iter {iteration:4d}: change={field_change:.2e}, continuity={continuity:.2e}, "
                      "dp={pressure_drop_pa:.6g} Pa".format(**row), flush=True)
                file.flush()

        if args.restart:
            with np.load(args.restart) as saved:
                summary, fields = solve(mesh, inlet, args.steps, settings, monitor, restart=saved)
        else:
            summary, fields = solve(mesh, inlet, args.steps, settings, monitor)
    if exact:
        reference, _ = poiseuille(mesh.centers, width, speed)
        x = mesh.centers[:, 0]
        interior = (x > .005) & (x < .015)
        pressure_gradient = -np.polyfit(x[interior], fields["pressure_pa"][interior], 1)[0]
        exact_gradient = settings.viscosity*speed/(width**2*average)
        summary["duct_verification"] = dict(
            relative_velocity_l2=float(np.linalg.norm(fields["velocity"][interior, 0]-reference[interior]) /
                                       np.linalg.norm(reference[interior])),
            pressure_gradient_pa_per_m=float(pressure_gradient), exact_gradient_pa_per_m=float(exact_gradient),
            relative_pressure_gradient_error=float(abs(pressure_gradient/exact_gradient-1)))
    flows = summary["flow_m3_s"]
    if "side_outlet" in flows and flows["main_outlet"] > 0:
        summary["side_to_main_ratio"] = flows["side_outlet"]/flows["main_outlet"]
    summary.update(cells=len(mesh.volumes), turbulence="SST-style k-omega" if settings.turbulence else "laminar",
                   configuration=asdict(settings), inlet_speed_m_s=speed,
                   source_mesh=str(args.mesh) if args.mesh else "square duct",
                   restart_file=str(args.restart) if args.restart else None)
    report = json.dumps(summary, indent=2)
    (args.output / "summary.json").write_text(report + "\n", encoding="utf-8")
    np.savez_compressed(args.output / "fields.npz", **fields)
    print(report)


if __name__ == "__main__":
    main()
