"""Build display-only pipe sections from the archived contour polygons."""

from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "assets/viewer"
INPUT = "evidence/verification/input"
STAGES = ("stage_01_beta_0", "stage_02_beta_2", "stage_03_beta_8")
CONTOURS = {
    "original": f"{INPUT}/{STAGES[2]}/original/geometry.json",
    "raw": f"{INPUT}/{STAGES[2]}/selected/geometry.json",
    "smoothed": f"{INPUT}/smoothed_geometry.json",
    "stokes": f"evidence/nonlinear/stokes_run/{STAGES[2]}/selected/geometry.json",
}


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sections(path):
    record = read(path)
    polygon = record.get("geometry", record)
    ring, = polygon["coordinates"]
    walls = [(a, b) for a, b in zip(ring, ring[1:]) if a[1] != b[1]]
    levels = sorted({4.0, 6.0, 7.0, *(y for _, y in ring if 4 < y < 7)})
    rows = []
    for lo, hi in zip(levels, levels[1:]):
        middle = (lo + hi) / 2
        edges = [(a, b) for a, b in walls if min(a[1], b[1]) < middle < max(a[1], b[1])]
        if len(edges) < 2 or (len(edges) != 2 and "stokes_run" not in path):
            raise ValueError(f"{path}: expected two branch walls at y={middle}")
        for y in (lo, hi):
            crossings = [a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
                         for a, b in edges]
            # Stokes contains a short re-entrant slice. Use its outer envelope for the loft.
            left, right = min(crossings), max(crossings)
            row = [round(s, 10) for s in (y, left, right)]
            if not rows or row != rows[-1]:
                rows.append(row)
    if "stokes_run" in path:
        for index in range(1, len(rows)):
            previous, current = rows[index - 1:index + 1]
            if previous[0] == current[0] and current[0] < 6:
                # Bridge re-entrant slice changes over 0.01 mm for display.
                previous[0] -= min(0.005, (previous[0] - rows[index - 2][0]) / 3)
                current[0] += min(0.005, (rows[index + 1][0] - current[0]) / 3)
    if rows[-1] != [7, 9.5, 10.5] or any(
        row[1:] != [9.5, 10.5] for row in rows if row[0] > 6
    ):
        raise ValueError(f"{path}: fixed outlet dimensions changed")
    return {"source": path, "sections_mm": rows}


def build():
    profiles = {name: sections(path) for name, path in CONTOURS.items()}
    history = []
    for position, stage in enumerate(STAGES):
        for variant in (("original", "selected") if position == 0 else ("selected",)):
            prefix = f"evidence/pilot/local/{stage}/images"
            history.append({
                "label": "Original" if variant == "original" else f"β = {(0, 2, 8)[position]}",
                "geometry": f"{prefix}/{'original' if variant == 'original' else 'optimized'}_physical_geometry.png",
                "density": f"{prefix}/{'initial' if variant == 'original' else 'optimized'}_density.png",
            })
    (OUTPUT / "pipe-profiles.json").write_text(
        json.dumps({"wall_mm": 0.22, "profiles": profiles, "history": history},
                   ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    build()
