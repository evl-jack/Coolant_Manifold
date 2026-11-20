"""Run the explicit reference and density models; write a compact comparison."""
import csv
from pathlib import Path

try:
    from .solve import solve
except ImportError:
    from solve import solve


def run(iterations=None):
    rows = []
    for name, include_side in (("explicit", True), ("A", False), ("B", False), ("C", False)):
        result = solve(iterations, include_side_outlet=include_side, model=name)
        flows = {tag: sum(f for f, face in zip(result["boundary_flux"], result["boundary"])
                          if face["tag"] == tag) for tag in ("inlet", "main_outlet", "side_outlet")}
        rows.append({"model": name, "converged": result["converged"], "iterations": result["iterations"],
                     "mass_residual": result["mass_residual"], "global_mass_residual": result["global_mass_residual"],
                     "inlet_kg_s_m": -flows["inlet"], "main_outlet_kg_s_m": flows["main_outlet"],
                     "side_outlet_kg_s_m": flows["side_outlet"],
                     "side_split": flows["side_outlet"] / max(-flows["inlet"], 1e-30)})
    path = Path(__file__).resolve().parent / "results" / "model_comparison.csv"
    path.parent.mkdir(exist_ok=True)
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Comparison written to {path}")


if __name__ == "__main__":
    run()
