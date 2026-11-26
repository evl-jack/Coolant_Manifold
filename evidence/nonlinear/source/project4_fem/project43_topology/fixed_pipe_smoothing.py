"""Postprocess the archived beta=8 contour and replay three fitted walls."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, box, mapping, shape

from .fixed_pipe_pilot_model import original_geometry, DESIGN_BOUNDS_MM, MAIN_BOUNDS_MM, COLLAR_BOUNDS_MM
from .wall_displacement_geometry_v1 import _turning_report
from .reconstruction import Reconstruction, solve_fitted
from .evaluator import HEIGHT_M, VELOCITY_SCALE


def smooth(raw):
    vertices = np.array(raw.exterior.coords[:-1])
    anchors = [(8., 4.), (9.5, 6.), (10.5, 6.), (12., 4.)]
    ids = [np.flatnonzero(np.all(vertices == p, axis=1)).item() for p in anchors]
    for a, b in zip(ids[::2], ids[1::2]):
        wall = vertices[a:b+1]
        for _ in range(10):
            wall[1:-1] = .5 * wall[1:-1] + .25 * (wall[:-2] + wall[2:])
            # Constrain the end segments to the adjoining fixed wall tangents.
            for end, near in ((0, 1), (-1, -2)):
                coord = 1 if wall[end, 1] == 4. else 0
                wall[near, coord] = wall[end, coord]
    return Polygon(vertices), ids


def main():
    root = Path('output/fixed_pipe_pilot_local/stage_03_beta_8/selected')
    folder = Path('output/fixed_pipe_fitted_smoothing')
    folder.mkdir(exist_ok=True)
    original = original_geometry().geometry
    raw = shape(json.loads((root / 'geometry.json').read_text())['geometry'])
    smoothed, ends = smooth(raw)
    rectangle = box(*DESIGN_BOUNDS_MM)
    cases = {'Original': original, 'Raw optimized': raw, 'Smoothed optimized': smoothed}
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6), layout='constrained', sharex=True, sharey=True)
    for ax, (label, polygon) in zip(axes, cases.items()):
        if not (polygon.is_valid and polygon.covers(box(*MAIN_BOUNDS_MM)) and
                polygon.covers(box(*COLLAR_BOUNDS_MM)) and
                polygon.symmetric_difference(original).difference(rectangle).area == 0.):
            raise ValueError(f'Invalid fixed geometry in {label}')
        ax.fill(*polygon.exterior.xy, color='#acd8da', edgecolor='#17546a')
        ax.set(title=label, aspect='equal', xlim=(7.5,12.5), ylim=(3.5,7.2), xlabel='x (mm)')
    axes[0].set_ylabel('y (mm)')
    fig.savefig(folder / 'three_panel.png', dpi=165)
    plt.close(fig)
    (folder / 'smoothed_geometry.json').write_text(json.dumps(mapping(smoothed)))
    def evaluate(item):
        name, polygon = item
        rec = Reconstruction(polygon, (polygon,), (0,0,0), ())
        result, _, _ = solve_fitted(rec, .25)
        assert result['state_passed']
        power = result['mechanical_power_W_per_m']
        turn = _turning_report(original, rec, rectangle)
        result.update(resistance=power/(result['inlet_flux']*HEIGHT_M*VELOCITY_SCALE)**2,
                      max_turn_degrees=turn['max_turn_degrees'])
        if name != 'Original':
            v = np.array(polygon.exterior.coords[:-1])
            d1 = v[ends] - v[np.array(ends)-1]
            d2 = v[np.array(ends)+1] - v[ends]
            cosine = (d1*d2).sum(axis=1)/(np.linalg.norm(d1,axis=1)*np.linalg.norm(d2,axis=1))
            result['join_degrees'] = np.degrees(np.arccos(np.clip(cosine,-1,1))).tolist()
        return name, result
    with ThreadPoolExecutor(max_workers=3) as pool:
        reports = dict(pool.map(evaluate,cases.items()))
    (folder / 'results.json').write_text(json.dumps(reports,indent=2))
    print(json.dumps(reports))


if __name__ == '__main__':
    main()
