"""Aggregate the fixed-pipe mesh study without dropping failed or absent levels."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .mesh_cases import CASES, LEVELS


POWER_TOLERANCE_PERCENT = .5
SIDE_TOLERANCE_PP = .020
PAIRS = {'immersed-raw': 'immersed-original', 'fitted-raw': 'fitted-original',
         'fitted-smoothed': 'fitted-original'}
POWER = 'mechanical_power_W_per_m'
SIDE = 'side_fraction'
COLORS = dict(zip(CASES, ('#336086', '#00786E', '#244574', '#CA6322', '#943A6C')))


def differences(series):
    checks = []
    for coarse, fine in zip(series[:-1], series[1:]):
        if coarse['state_passed'] and fine['state_passed']:
            power = abs(fine[POWER] - coarse[POWER])
            relative = 100. * power / abs(fine[POWER])
            side_pp = 100. * abs(fine[SIDE] - coarse[SIDE])
            checks.append({
                'coarse_mm': coarse['mesh_target_mm'], 'fine_mm': fine['mesh_target_mm'],
                'power_difference_W_per_m': power, 'power_difference_percent': relative,
                'side_difference_pp': side_pp,
                'passed': relative <= POWER_TOLERANCE_PERCENT and side_pp <= SIDE_TOLERANCE_PP,
            })
        else:
            checks.append({'coarse_mm': coarse['mesh_target_mm'], 'fine_mm': fine['mesh_target_mm'],
                           'passed': False, 'incomplete': True})
    return checks


def apparent_order(values):
    """Withhold extrapolation if three values lack a decreasing monotone trend."""
    a, b = np.diff(values[-3:])
    if a * b <= 0. or abs(b) >= abs(a):
        return None
    p = float(np.log2(abs(a / b)))
    return {'apparent_order': p, 'richardson_estimate': float(values[-1] + b / (2. ** p - 1.)),
            'qualification': 'Three-level monotone diagnostic for the nominal 2:1 spacing'}


def comparison(reference, candidate, reference_screen, candidate_screen, fitted):
    rows = []
    for original, optimized in zip(reference, candidate):
        if not (original['state_passed'] and optimized['state_passed']):
            continue
        rows.append({
            'spacing_mm': original['mesh_target_mm'],
            'power_saving_W_per_m': original[POWER] - optimized[POWER],
            'power_saving_percent': 100. * (original[POWER] - optimized[POWER]) / original[POWER],
            'side_increase_pp': 100. * (optimized[SIDE] - original[SIDE]),
            'combined_raw_energy_residual_W_per_m': (abs(original['energy']['raw_closure'])
                                                      + abs(optimized['energy']['raw_closure'])),
        })
    result = {'by_mesh': rows, 'screen_passed': reference_screen and candidate_screen,
              'power_improvement_resolved': False, 'reference_side_flow_maintained': None,
              'fitted_performance_verified': False}
    if len(rows) != len(reference) or len(rows) < 3:
        return result
    last = rows[-1]
    latest_reference, latest_candidate = differences(reference)[-2:], differences(candidate)[-2:]
    power_uncertainty = max(
        np.max(np.abs(np.diff([row['power_saving_W_per_m'] for row in rows[-3:]]))),
        *(r['power_difference_W_per_m'] + c['power_difference_W_per_m']
          for r, c in zip(latest_reference, latest_candidate)),
    )
    side_uncertainty = max(
        np.max(np.abs(np.diff([row['side_increase_pp'] for row in rows[-3:]]))),
        *(r['side_difference_pp'] + c['side_difference_pp'] for r, c in zip(latest_reference, latest_candidate)),
    )
    result['paired_power_uncertainty_screen_W_per_m'] = float(power_uncertainty)
    result['paired_side_uncertainty_screen_pp'] = float(side_uncertainty)
    result['energy_residual_below_power_saving'] = (last['combined_raw_energy_residual_W_per_m']
                                                   < last['power_saving_W_per_m'])
    result['power_improvement_resolved'] = bool(result['screen_passed']
                                                and last['power_saving_W_per_m'] > power_uncertainty
                                                and result['energy_residual_below_power_saving'])
    result['reference_side_flow_maintained'] = bool(result['screen_passed']
                                                    and last['side_increase_pp'] >= side_uncertainty)
    result['fitted_performance_verified'] = bool(fitted and result['power_improvement_resolved']
                                                 and result['reference_side_flow_maintained'])
    result['paired_power_apparent_convergence'] = apparent_order([row['power_saving_W_per_m'] for row in rows])
    return result


def load_results(root, fine):
    series = {}
    for case in CASES:
        series[case] = []
        for h in (*LEVELS, .03125) if fine else LEVELS:
            path = root / f'{case}-{h:g}' / 'result.json'
            if path.is_file():
                row = json.loads(path.read_text(encoding='utf-8'))
                if row['case'] != case or row['mesh_target_mm'] != h:
                    raise ValueError(f'Case metadata does not match directory {path.parent}')
            else:
                row = {'case': case, 'mesh_target_mm': h, 'state_passed': False,
                       'failure': 'No result file. Check the Slurm task log and accounting.'}
            series[case].append(row)
        for field in ('geometry_sha256', 'design_fields_sha256'):
            if len({row[field] for row in series[case] if row['state_passed']}) > 1:
                raise ValueError(f'{case} used different {field} between refinements')
    return series


def save_figure(fig, output, name):
    fig.savefig(output / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(output / f'{name}.png', dpi=220, bbox_inches='tight')
    plt.close(fig)


def convergence_plot(series, paired, output):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4), layout='constrained')
    for key, rows in series.items():
        h = [r['mesh_target_mm'] for r in rows]
        label = key.replace('-', ' ').capitalize()
        for ax, field, scale in zip(axes[:2], (POWER, SIDE), (1e6, 100.)):
            ax.plot(h, [scale * r[field] if r['state_passed'] else np.nan for r in rows], 'o-',
                    color=COLORS[key], label=label)
    for key, result in paired.items():
        rows = result['by_mesh']
        if rows:
            axes[2].plot([r['spacing_mm'] for r in rows], [1e6 * r['power_saving_W_per_m'] for r in rows],
                         'o-', color=COLORS[key])
    axes[0].legend(frameon=False, fontsize=8)
    axes[2].axhline(0., color='.35', linewidth=.7)
    for ax, label in zip(axes, ('Pumping power (µW/m)', 'Side flow (% of inlet)', 'Paired power saving (µW/m)')):
        ax.set(xscale='log', xlabel='Mesh target (mm)', ylabel=label)
        ax.set_xticks(LEVELS, labels=[f'{h:g}' for h in LEVELS])
        ax.tick_params(which='minor', labelbottom=False)
        ax.invert_xaxis()
        ax.grid(alpha=.20)
    save_figure(fig, output, 'mesh_convergence')


def field_plot(family, series, root, output):
    fields = []
    for key, rows in series.items():
        if not key.startswith(family):
            continue
        for row in reversed(rows):
            if row['state_passed']:
                h = row['mesh_target_mm']
                with np.load(root / f'{key}-{h:g}' / 'fields.npz') as data:
                    fields.append((key, h, dict(data)))
                break
    if not fields:
        return
    fig, axes = plt.subplots(len(fields), 2, figsize=(9, 2.5 * len(fields)), squeeze=False, layout='constrained')
    pressures = [d['pressure_pa'] for _, _, d in fields]
    speeds = [np.linalg.norm(d['velocity_m_per_s'], axis=1) for _, _, d in fields]
    p_range = (min(p.min() for p in pressures), max(p.max() for p in pressures))
    speed_max = max(s.max() for s in speeds)
    for i, (key, h, data) in enumerate(fields):
        dofs = data['velocity_dofs']
        subt = np.concatenate([dofs[:, inds] for inds in ((0, 3, 5), (3, 1, 4), (5, 4, 2), (3, 4, 5))])
        p = axes[i, 0].tripcolor(*data['pressure_points_mm'].T, data['triangles'], pressures[i],
                                shading='gouraud', cmap='viridis', vmin=p_range[0], vmax=p_range[1])
        u = axes[i, 1].tripcolor(*data['velocity_points_mm'].T, subt, speeds[i], shading='gouraud',
                                cmap='magma', vmin=0., vmax=speed_max, rasterized=True)
        p.set_rasterized(True)
        axes[i, 0].set_ylabel(f"{key.split('-')[-1].capitalize()}\nh = {h:g} mm")
        for ax in axes[i]:
            ax.set(aspect='equal', xlim=(7., 13.), ylim=(2.5, 7.08))
            ax.set_xticks([8., 10., 12.])
            if i == len(fields) - 1:
                ax.set_xlabel('x (mm)')
    axes[0, 0].set_title('Pressure (Pa)')
    axes[0, 1].set_title('Speed (m/s)')
    fig.colorbar(p, ax=axes[:, 0], shrink=.75)
    fig.colorbar(u, ax=axes[:, 1], shrink=.75)
    save_figure(fig, output, f'{family}_pressure_speed')


def report(root, output, fine):
    output.mkdir(parents=True, exist_ok=False)
    series = load_results(root, fine)
    checks = {key: differences(rows) for key, rows in series.items()}
    screens = {key: all(check['passed'] for check in items[-2:]) for key, items in checks.items()}
    paired = {key: comparison(series[base], series[key], screens[base], screens[key], key.startswith('fitted-'))
              for key, base in PAIRS.items()}
    asymptotic = {key: {field: apparent_order([row[field] for row in rows]) for field in (POWER, SIDE)}
                  for key, rows in series.items() if all(row['state_passed'] for row in rows[-3:])}
    result = {'thresholds': {'power_percent': POWER_TOLERANCE_PERCENT, 'side_pp': SIDE_TOLERANCE_PP},
              'series': series, 'successive_comparisons': checks, 'two_successive_passes': screens,
              'apparent_convergence': asymptotic, 'paired_comparisons': paired,
              'uncertainty_definition': 'Maximum of the two latest paired changes and corresponding sums of '
              'absolute individual changes. It is a screening estimate rather than a rigorous error bound.'}
    (output / 'mesh_report.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    columns = ('case', 'mesh_target_mm', 'state_passed', 'triangles', 'dofs', POWER, SIDE,
               'physical_fluid_area_mm2', 'hydraulic_resistance_pa_s_per_m2', 'nonlinear_residual',
               'relative_mass_imbalance', 'energy_raw_closure_W_per_m', 'failure')
    with (output / 'mesh_table.csv').open('w', encoding='utf-8', newline='') as file:
        writer = csv.DictWriter(file, columns, extrasaction='ignore')
        writer.writeheader()
        for rows in series.values():
            for row in rows:
                writer.writerow(dict(row, energy_raw_closure_W_per_m=row.get('energy', {}).get('raw_closure')))
    convergence_plot(series, paired, output)
    for family in ('immersed', 'fitted'):
        field_plot(family, series, root, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh_case_results', type=Path)
    parser.add_argument('report_output', type=Path)
    parser.add_argument('--additional-level', action='store_true')
    args = parser.parse_args()
    report(args.mesh_case_results, args.report_output, args.additional_level)


if __name__ == '__main__':
    main()
