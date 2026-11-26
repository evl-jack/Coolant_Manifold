"""Build the static report from archived measurements without running CFD."""

import argparse
import hashlib
import html
import json
import math
from pathlib import Path
import re
import shutil
import sys

import numpy as np
from shapely import constrained_delaunay_triangles
from shapely.geometry import shape

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EVIDENCE = ROOT / "output"


def read(path, wrapped=False):
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if wrapped:
        if not value["completed"]:
            raise ValueError(f"Incomplete record: {path}")
        return value["results"]
    return value


def dump(path, data):
    path.write_text(json.dumps(data, allow_nan=False, separators=(",", ":")), encoding="utf-8")


def contour_data(geometry):
    polygon = shape(geometry)
    if not polygon.is_valid or polygon.geom_type != "Polygon":
        raise ValueError("Expected a valid single fluid polygon")
    caps = constrained_delaunay_triangles(polygon)
    if not math.isclose(sum(cap.area for cap in caps.geoms), polygon.area, rel_tol=1e-10):
        raise ValueError("Cap triangulation changed fluid area")
    rings = geometry["coordinates"]
    levels = sorted({y for ring in rings for _, y in ring if 4 <= y <= 7} | {4., 7.})
    profile = []
    for y in levels:
        # Take the upper one-sided limit at the shoulder; exclude the main-pipe ledge.
        sample_y = y + 1e-8 if y == 4 else y - 1e-8 if y == 7 else y
        xs = [x1 + (sample_y - y1) * (x2 - x1) / (y2 - y1)
              for ring in rings for (x1, y1), (x2, y2) in zip(ring, ring[1:])
              if min(y1, y2) <= sample_y < max(y1, y2)]
        if len(xs) != 2:
            raise ValueError(f"Expected two transition walls at y={y}, found {len(xs)}")
        left, right = sorted(xs)
        profile.append([y, left, right, (right - left) / 2])
    return {"rings": rings, "profile": profile,
            "caps": [list(cap.exterior.coords)[:3] for cap in caps.geoms]}


def metric(value, units, source, experiment, status, settings):
    return dict(value=value, units=units, source=source, experiment=experiment,
                status=status, settings=settings)


def export_data(evidence, assets):
    pilot = evidence / "fixed_pipe_pilot_local"
    fitted = evidence / "fixed_pipe_fitted_smoothing"
    settings = read(pilot / "configuration.json", True)
    # The run metadata contains a workstation name, which the public site does not need.
    settings.pop("environment", None)
    checkpoint_rows = read(pilot / "summary.json", True)["stages"]
    checkpoints = []
    for row in checkpoint_rows:
        stage = row["stage"]
        beta = row["settings"]["beta"]
        source = pilot / f"stage_{stage:02}_beta_{beta:g}"
        target = f"checkpoint-{stage}.json"
        original = read(source / "original/metrics.json", True)
        selected = read(source / "selected/metrics.json", True)
        dump(assets / target, {
            key: contour_data(read(source / key / "geometry.json")["geometry"])
            for key in ("original", "selected")
        })
        checkpoints.append({"stage": stage, "beta": beta, "geometry": f"assets/{target}",
                            "status": row["mma_native_message"], "evaluations": row["mma_unique_evaluations"],
                            "reference": original, "selected": selected})
    final = pilot / "stage_03_beta_8"
    with np.load(pilot / "mesh.npz") as mesh:
        points, triangles, dofs = (mesh[name] for name in ("points", "triangles", "velocity_dofs"))
        nv = len(mesh["velocity_points"])
        assert np.array_equal(dofs[:, :3], triangles)
    velocity_scale = 50 * .001002 / (998.2 * .004)
    pressure_scale = .001002 * velocity_scale / .004
    fields = {}
    for key in ("original", "selected"):
        with np.load(final / key / "fields.npz") as archive:
            state = archive["state"]
            if len(state) != 2 * nv + len(points) or not np.isfinite(state).all():
                raise ValueError("Unexpected Taylor–Hood state layout")
            fields[key] = dict(ux=(state[:nv] * velocity_scale).tolist(),
                               uy=(state[nv:2*nv] * velocity_scale).tolist(),
                               pressure=(state[2*nv:] * pressure_scale).tolist())
    dump(assets / "fields.json", dict(points=(4 * points).tolist(), triangles=triangles.tolist(),
                                      velocityDofs=dofs.tolist(), fields=fields))
    smoothing = read(fitted / "results.json")
    measurements = {}
    mappings = {"area": ("physical_fluid_area_mm2", 1, "mm²"), "turn": ("max_turn_degrees", 1, "°"),
                "power": ("mechanical_power_W_per_m", 1e6, "µW/m"),
                "resistance": ("resistance", 1, "Pa·s/m²"), "flow": ("side_share", 100, "%"),
                "mesh": ("triangles", 1, "triangles")}
    for design, record in smoothing.items():
        measurements[design] = {
            key: metric(record[field] * factor, units, "fixed_pipe_fitted_smoothing/results.json",
                        "fixed-pipe-fitted-comparison", "preliminary",
                        dict(reynolds=50, target_spacing_mm=record["spacing_mm"]))
            for key, (field, factor, units) in mappings.items()
        }
    journal = [json.loads(line) for line in (pilot / "progress.jsonl").read_text(encoding="utf-8").splitlines() if line]
    kkt_path = evidence / "fixed_pipe_kkt_endpoint/fixed_pipe_kkt.json"
    kkt = read(kkt_path, True) if kkt_path.exists() else None
    kkt_rows = [{key: row[key] for key in ("verified", "multipliers", "stationarity_inf",
                 "complementarity_inf", "primal_violation", "dual_violation")}
                for row in kkt["stages"]] if kkt else None
    shutil.copyfile(fitted / "three_panel.png", assets / "three-panel.png")
    smoke = evidence / "fixed_pipe_independent_verification/adjoint_taylor_remainders.png"
    if smoke.exists():
        shutil.copyfile(smoke, assets / "adjoint-smoke.png")
    evidence_file = evidence / "fixed_pipe_site_release.json"
    release = read(evidence_file) if evidence_file.exists() else {}
    required = {"finalTaylor", "cost", "scaling", "independentKkt", "immersedMesh", "fittedMesh", "stokes"}
    for name, record in release.items():
        if name not in required or not all(record.get(k) is not None for k in
                                           ("source", "experiment", "status", "settings", "columns", "rows")):
            raise ValueError(f"Invalid numerical release section: {name}")
        if any(len(row) != len(record["columns"]) for row in record["rows"]):
            raise ValueError(f"Release table width mismatch: {name}")
        figure = record.get("figure")
        if figure:
            if Path(figure).name != figure or Path(figure).suffix.lower() not in (".png", ".webp", ".jpg"):
                raise ValueError("Release figure must be an approved raster basename")
            shutil.copyfile(evidence / figure, assets / figure)
    dump(assets / "results.json", dict(schema=1, settings=settings, fitted=measurements,
         checkpoints=checkpoints, history=journal, fields="assets/fields.json", kkt=kkt_rows,
         kktTolerance=kkt["tolerances"] if kkt else None, adjointSmoke=smoke.exists(),
         pending={name: release.get(name) for name in sorted(required)},
         sources={"pilot": "fixed_pipe_pilot_local", "history": "fixed_pipe_pilot_local/progress.jsonl",
                  "kkt": "fixed_pipe_kkt_endpoint/fixed_pipe_kkt.json",
                  "fieldPressureScalePa": pressure_scale, "fieldVelocityScaleMps": velocity_scale,
                  "fieldSHA256": hashlib.sha256((assets / "fields.json").read_bytes()).hexdigest()}))


def build(evidence, output):
    local_dependencies = HERE.parent / "project4/.build-deps"
    if local_dependencies.exists():
        sys.path.insert(0, str(local_dependencies))
    import markdown2
    from latex2mathml.converter import convert

    output.mkdir(parents=True, exist_ok=True)
    assets = output / "assets"
    assets.mkdir(exist_ok=True)
    export_data(evidence, assets)
    for name in ("site.css", "viewer.js", "site.js"):
        shutil.copyfile(HERE / name, output / name)
    (output / ".nojekyll").touch()
    template = (HERE / "layout.html").read_text(encoding="utf-8")
    for source, destination, title in (("REPORT.md", "index.html", "Fixed-pipe coolant manifold"),
                                        ("HANDOFF.md", "handoff.html", "Project handoff"),
                                        ("OPTIMIZATION.md", "optimization.html", "Optimization formulation")):
        equations = []

        def preserve(match):
            key = f"PROJECTFOURMATH{len(equations):04d}"
            # Bracket-delimited display equations avoid conflicts with Markdown.
            equations.append((key, convert(match[1].strip(), display="block")))
            return key

        markdown = (HERE / source).read_text(encoding="utf-8")
        markdown = re.sub(r"\\\[([\s\S]+?)\\\]", preserve, markdown)
        body = markdown2.markdown(markdown, extras=["tables", "fenced-code-blocks", "header-ids"])
        for key, equation in equations:
            body = body.replace(f"<p>{key}</p>", equation).replace(key, equation)
        (output / destination).write_text(template.replace("{{title}}", html.escape(title))
                                         .replace("{{content}}", body), encoding="utf-8")


def offline_preview(output, target):
    """Write a single-file preview for devices that cannot serve localhost."""
    page = (output / "index.html").read_text(encoding="utf-8")
    css = (output / "site.css").read_text(encoding="utf-8")
    viewer = (output / "viewer.js").read_text(encoding="utf-8").replace("export ", "")
    script = (output / "site.js").read_text(encoding="utf-8").split("\n", 1)[1]
    script = script.replace('  const response = await fetch(path);',
                            '  if (Object.hasOwn(embeddedData, path)) return JSON.parse(embeddedData[path]);\n'
                            '  const response = await fetch(path);')
    assets = output / "assets"
    embedded = {f"assets/{path.name}": path.read_text(encoding="utf-8")
                for path in assets.glob("*.json")}
    script = "const embeddedData = " + json.dumps(embedded).replace("<", "\\u003c") + ";\n" + viewer + script
    page = page.replace('<link rel="stylesheet" href="site.css">', "<style>" + css + "</style>")
    page = page.replace('<script type="module" src="site.js"></script>', '<script type="module">' + script + '</script>')
    # Keep the images in the corresponding packaged asset directory.
    if target.parent == output:
        target.write_text(page, encoding="utf-8")
    else:
        raise ValueError("Offline preview must sit next to the generated assets")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=EVIDENCE)
    parser.add_argument("--output", type=Path, default=HERE / "dist")
    parser.add_argument("--offline-preview", action="store_true")
    args = parser.parse_args()
    build(args.evidence, args.output)
    if args.offline_preview:
        offline_preview(args.output, args.output / "offline-preview.html")
