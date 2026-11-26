"""Cross-check the raw Stokes contour under fitted-wall Navier–Stokes."""

import argparse
import csv
import json
import os
from pathlib import Path

from .mesh_cases import LEVELS, file_hash, run_case
from .mesh_report import POWER, SIDE, comparison, differences

CASE = 'fitted-stokes'
CASES = {'original': 'fitted-original', 'navier_stokes': 'fitted-raw', 'stokes': CASE}
ROOTS = {
    'original': ('original', 'geometry.json'),
    'navier_stokes': ('selected', 'geometry.json'),
}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def fitted_run(stage, contour, out, spacing):
    out.mkdir(parents=True, exist_ok=False)
    print(f'Starting {CASE} on target h={spacing:g} mm', flush=True)
    try:
        row = run_case(stage, None, out, CASE, spacing, fitted_contour=contour)
    except Exception as exc:
        write_json(out / 'result.json', {
            'case': CASE, 'mesh_target_mm': spacing, 'state_passed': False,
            'failure': f'{type(exc).__name__}: {exc}', 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
        })
        raise
    write_json(out / 'result.json', row)
    print(json.dumps(row, allow_nan=False), flush=True)
    return row['state_passed']


def read_series(case, directory, expected_geometry_hash):
    rows = []
    for spacing in LEVELS:
        file = directory / f'{case}-{spacing:g}' / 'result.json'
        if file.is_file():
            row = json.loads(file.read_text(encoding='utf-8'))
            if row['case'] != case or row['mesh_target_mm'] != spacing:
                raise ValueError(f'Fitted case metadata disagrees with {file}')
            if row['state_passed'] and row['design_fields_sha256'] is not None:
                raise ValueError(f'Fitted case unexpectedly uses immersed density at {file}')
        else:
            row = {'case': case, 'mesh_target_mm': spacing, 'state_passed': False,
                   'failure': 'Missing result. See the Slurm task log and accounting.'}
        rows.append(row)
    hashes = {row['geometry_sha256'] for row in rows if row['state_passed']}
    if len(hashes) > 1 or (expected_geometry_hash and hashes - {expected_geometry_hash}):
        raise ValueError(f'{case} geometry checksum disagrees with its frozen contour')
    return rows


def summarize_pair(navier, stokes, navier_screen, stokes_screen):
    signed = [
        {'mesh_target_mm': s['mesh_target_mm'],
         'stokes_minus_navier_stokes_power_W_per_m': s[POWER] - n[POWER],
         'stokes_minus_navier_stokes_side_flow_pp': 100. * (s[SIDE] - n[SIDE])}
        for n, s in zip(navier, stokes) if n['state_passed'] and s['state_passed']
    ]
    paired = comparison(navier, stokes, navier_screen, stokes_screen, True)
    power = paired.get('paired_power_uncertainty_screen_W_per_m')
    side = paired.get('paired_side_uncertainty_screen_pp')
    energy = paired['by_mesh'][-1]['combined_raw_energy_residual_W_per_m'] if signed else None
    finest = signed[-1] if signed else None
    complete = len(signed) == len(LEVELS)
    power_resolved = bool(complete and paired['screen_passed']
                          and abs(finest['stokes_minus_navier_stokes_power_W_per_m']) > max(power, energy))
    side_resolved = bool(complete and paired['screen_passed']
                         and abs(finest['stokes_minus_navier_stokes_side_flow_pp']) > side)
    return {
        'signed_difference_definition': 'Stokes-optimized contour minus Navier–Stokes-optimized contour, both evaluated under fitted-wall Navier–Stokes',
        'by_mesh': signed,
        'finest_paired': finest,
        'mesh_screen_passed': paired['screen_passed'],
        'power_uncertainty_screen_W_per_m': power,
        'side_flow_uncertainty_screen_pp': side,
        'combined_finest_raw_energy_residual_W_per_m': energy,
        'power_difference_resolved': power_resolved,
        'side_flow_difference_resolved': side_resolved,
        'power_comparison_status': (
            'incomplete' if not complete else
            'unresolved' if not power_resolved else
            'stokes_contour_lower_power' if finest['stokes_minus_navier_stokes_power_W_per_m'] < 0. else
            'navier_stokes_contour_lower_power'
        ),
    }


def fitted_report(stage, archived, stokes_results, out):
    series = {}
    for name, case in CASES.items():
        root = stokes_results if name == 'stokes' else archived
        expected = None if name == 'stokes' else file_hash(stage.joinpath(*ROOTS[name]))
        series[name] = read_series(case, root, expected)
    checks = {name: differences(rows) for name, rows in series.items()}
    screen = {name: all(item['passed'] for item in results) for name, results in checks.items()}
    result = {
        'mesh_levels_mm': LEVELS,
        'physics': 'Fitted-wall Navier–Stokes, Re = 50, common three-port conditions and P2/P1 FEM',
        'archived_fitted_settings': (
            'Reused the original and raw Navier–Stokes contours by exact geometry-file SHA-256. '
            'New cases use the same fitted mesh constructor, nonlinear continuation, boundary data, power, '
            'and residual checks via the archived mesh-audit solver path.'
        ),
        'series': series,
        'successive_comparisons': checks,
        'mesh_screen_passed': screen,
        'comparison': summarize_pair(series['navier_stokes'], series['stokes'],
                                     screen['navier_stokes'], screen['stokes']),
        'stokes_vs_original': comparison(series['original'], series['stokes'],
                                        screen['original'], screen['stokes'], True),
        'uncertainty_definition': (
            'Maximum of the latest two paired differences and corresponding sums of absolute individual '
            'changes. Both contours must pass the existing mesh screen. A resolved power difference must '
            'also exceed the combined raw energy residual. This is a screening estimate.'
        ),
    }
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / 'fitted_comparison.json', result)
    columns = ('design', 'case', 'mesh_target_mm', 'state_passed', 'triangles', 'dofs', POWER, SIDE,
               'hydraulic_resistance_pa_s_per_m2', 'relative_mass_imbalance', 'nonlinear_residual', 'failure')
    with (out / 'fitted_comparison.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, columns, extrasaction='ignore')
        writer.writeheader()
        for name, rows in series.items():
            writer.writerows(dict(row, design=name) for row in rows)
    print(json.dumps(result['comparison'], allow_nan=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    run = commands.add_parser('run')
    run.add_argument('archived_stage', type=Path)
    run.add_argument('stokes_contour', type=Path)
    run.add_argument('output', type=Path)
    run.add_argument('spacing', type=float, choices=LEVELS)
    report = commands.add_parser('report')
    report.add_argument('archived_stage', type=Path)
    report.add_argument('archived_mesh_results', type=Path)
    report.add_argument('stokes_mesh_results', type=Path)
    report.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.command == 'run':
        raise SystemExit(0 if fitted_run(args.archived_stage, args.stokes_contour, args.output, args.spacing) else 1)
    fitted_report(args.archived_stage, args.archived_mesh_results, args.stokes_mesh_results, args.output)


if __name__ == '__main__':
    main()
