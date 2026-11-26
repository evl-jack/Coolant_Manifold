"""Replay frozen fixed-pipe designs on immersed and body-fitted analysis meshes."""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import socket
from time import perf_counter

import numpy as np
from shapely.geometry import LineString, box, shape

from ...junction import boundary_data
from ...project43_power import port_power
from ...project43_power_audit import audit as energy_audit
from ..domain import REFERENCE_HEIGHT_MM
from ..evaluator import DesignTransfer, HEIGHT_M, POWER_SCALE, REYNOLDS, VELOCITY_SCALE, VISCOSITY
from ..fixed_pipe_pilot_model import COLLAR_BOUNDS_MM, DESIGN_BOUNDS_MM, MAIN_BOUNDS_MM, build_fixed_pipe
from ..leakage import _edge_flux, _edge_topology
from ..reconstruction import Reconstruction, _PORTS_MM, solve_fitted


CASES = ('immersed-original', 'immersed-raw', 'fitted-original', 'fitted-raw', 'fitted-smoothed')
LEVELS = (.25, .125, .0625)
ALPHA_MAX = 100000.


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def mesh_sizes(fem):
    vertices = fem.points[fem.triangles] * REFERENCE_HEIGHT_MM
    edge_lengths = np.linalg.norm(vertices - np.roll(vertices, 1, axis=1), axis=2)
    diameter = edge_lengths.max(axis=1)
    x, y = vertices.mean(axis=1).T
    junction = (x >= 8.) & (x <= 12.) & (y >= 3.5)
    collar = y >= 6.
    return {
        'minimum_edge_mm': float(edge_lengths.min()),
        'maximum_edge_mm': float(diameter.max()),
        'junction_maximum_edge_mm': float(diameter[junction].max()),
        'collar_maximum_edge_mm': float(diameter[collar].max()),
    }


def solid_transport(domain, projected, values, inlet):
    """Measure finite resistance transport in cells with projected density below 0.5."""
    fem = domain.fem
    indices = domain.cell_index
    editable = indices >= 0
    solid = np.zeros(len(indices), dtype=bool)
    solid[editable] = projected[indices[editable]] < .5
    velocity = values[:2 * fem.n_u].reshape(2, -1).T
    local_velocity = fem.basis @ velocity[fem.velocity_dofs[solid]]
    speed2 = np.sum(local_velocity ** 2, axis=2)
    area = fem.weights[solid].sum()
    rms = np.sqrt(np.sum(fem.weights[solid] * speed2) / area)

    vertices, first, last = _edge_topology(fem)
    interface = np.flatnonzero((first != last) & (solid[first] != solid[last]))
    solid_owner = np.where(solid[first[interface]], first[interface], last[interface])
    a, b = fem.points[vertices[interface]].transpose(1, 0, 2)
    delta = b - a
    normal = np.column_stack((delta[:, 1], -delta[:, 0]))
    normal /= np.linalg.norm(delta, axis=1)[:, None]
    solid_center = fem.points[fem.triangles[solid_owner]].mean(axis=1)
    normal[np.sum(normal * (solid_center - a), axis=1) < 0] *= -1.
    _, directional = _edge_flux(interface, vertices, fem.points, velocity, fem.n_p, normal)
    ingress, egress = directional.sum(axis=0) / inlet
    return {
        'definition': 'Element-aligned projected-density cells below 0.5',
        'solid_area_mm2': float(area * REFERENCE_HEIGHT_MM ** 2),
        'solid_rms_speed_over_inlet_bulk': float(rms),
        'solid_sampled_max_speed_over_inlet_bulk': float(np.sqrt(speed2.max())),
        'gross_ingress_over_inlet': float(ingress),
        'gross_egress_over_inlet': float(egress),
        'net_ingress_over_inlet': float(ingress - egress),
    }


def validation_geometry(path, baseline):
    data = json.loads(path.read_text(encoding='utf-8'))
    geometry = shape(data.get('geometry', data))
    if geometry.geom_type not in ('Polygon', 'MultiPolygon') or not geometry.is_valid or geometry.is_empty:
        raise ValueError(f'Frozen contour is not valid polygonal geometry at {path}')
    if not (geometry.covers(box(*MAIN_BOUNDS_MM)) and geometry.covers(box(*COLLAR_BOUNDS_MM))):
        raise ValueError(f'Frozen contour changed fixed fluid at {path}')
    if geometry.symmetric_difference(baseline).difference(box(*DESIGN_BOUNDS_MM)).area > 1e-11:
        raise ValueError(f'Frozen contour changed the geometry outside the envelope at {path}')
    polygons = (geometry,) if geometry.geom_type == 'Polygon' else tuple(geometry.geoms)
    connected = [i for i, polygon in enumerate(polygons)
                 if all(polygon.boundary.covers(LineString(port)) for port in _PORTS_MM)]
    if len(connected) != 1:
        raise ValueError(f'Frozen contour does not connect the three exact port traces at {path}')
    return Reconstruction(geometry, polygons, (connected[0],) * len(_PORTS_MM), ())


def run_case(stage_dir, smoothed, output, case, spacing, *, fitted_contour=None):
    started = perf_counter()
    name = case.split('-')[1]
    baseline_path = stage_dir / 'original' / 'geometry.json'
    baseline = shape(json.loads(baseline_path.read_text(encoding='utf-8'))['geometry'])
    folder = stage_dir / ('original' if name == 'original' else 'selected')
    if fitted_contour is not None and not case.startswith('fitted-'):
        raise ValueError('An explicit contour requires a fitted-wall case')
    geometry_path = fitted_contour or (smoothed if name == 'smoothed' else folder / 'geometry.json')
    rec = validation_geometry(geometry_path, baseline)
    leakage = None
    archived_fields = None
    if case.startswith('fitted-'):
        _, fem, state = solve_fitted(rec, spacing)
        _, _, ports, flux = boundary_data(fem)
        resistance = 0.
        load = fem.load()
    else:
        domain, pipeline = build_fixed_pipe(max_spacing_mm=spacing)
        fem, ports, flux = domain.fem, domain.ports, domain.flux
        archived_fields = folder / 'fields.npz'
        with np.load(archived_fields) as fields:
            design = fields['design']
            density = pipeline.evaluate(design, q=.10, beta=8., alpha_max=ALPHA_MAX)
            for field in ('raw', 'filtered', 'projected', 'resistance'):
                np.testing.assert_allclose(getattr(density, field), fields[field], rtol=1e-12, atol=1e-12)
            resistance = DesignTransfer(domain.cell_index, pipeline.grid.size).forward(density.resistance)
        load = fem.load()
        state = None
        for reynolds in (0., REYNOLDS / 2, REYNOLDS):
            fem.inertia = reynolds
            state = fem.solve(domain.fixed, domain.prescribed, load, resistance=resistance,
                              initial=None if state is None else state.values)
        leakage = solid_transport(domain, density.projected, state.values, float(flux[0] @ state.values))

    values = state.values
    inlet, main, side = np.asarray(flux @ values, dtype=float)
    mass_imbalance = abs(inlet - main - side) / inlet
    power = POWER_SCALE * sum(port_power(fem, ports, values, REYNOLDS)[0])
    energy = energy_audit(fem, ports, values, resistance, POWER_SCALE)
    np.testing.assert_allclose(energy['port_power'], power, rtol=1e-11)
    operator = fem.viscous_pressure + fem.drag_matrix(resistance)
    residual, _ = fem.residual(values, operator, load, tangent=False)
    weak = float(POWER_SCALE * residual[:2 * fem.n_u] @ values[:2 * fem.n_u])
    kinetic_power = sum(port[1] for port in energy['ports_inlet_main_side_pressure_kinetic_viscous'])
    energy['velocity_weighted_weak_reaction'] = weak
    energy['weak_form_closure'] = (energy['viscous_dissipation'] + energy['brinkman_dissipation']
                                   - energy['weighted_divergence'] - kinetic_power - weak)
    penetration_length_mm = 1000. * np.sqrt(VISCOSITY / ALPHA_MAX)
    sizes = mesh_sizes(fem)
    result = {
        'case': case,
        'mesh_target_mm': spacing,
        'geometry_sha256': file_hash(geometry_path),
        'design_fields_sha256': None if archived_fields is None else file_hash(archived_fields),
        'physical_fluid_area_mm2': rec.area_mm2,
        'port_connected_area_mm2': float(rec.port_component.area),
        'closed_pocket_area_mm2': float(sum(p.area for i, p in enumerate(rec.polygons)
                                             if i != rec.port_components[0])),
        'analysis_area_mm2': float(fem.weights.sum() * REFERENCE_HEIGHT_MM ** 2),
        'triangles': len(fem.triangles),
        'dofs': fem.size,
        'free_dofs': len(state.free),
        'mesh_sizes': sizes,
        'inlet_flow_m2_per_s': inlet * HEIGHT_M * VELOCITY_SCALE,
        'side_fraction': side / inlet,
        'main_fraction': main / inlet,
        'mechanical_power_W_per_m': power,
        'hydraulic_resistance_pa_s_per_m2': power / (inlet * HEIGHT_M * VELOCITY_SCALE) ** 2,
        'relative_mass_imbalance': mass_imbalance,
        'nonlinear_residual': state.residual_history[-1],
        'final_reynolds_residual_history': state.residual_history,
        'state_passed': bool(min(inlet, main, side, power) > 0. and mass_imbalance <= 1e-9
                             and state.residual_history[-1] <= 1e-10),
        'energy': energy,
        'solid_transport': leakage,
        'minimum_penetration_length_mm': None if leakage is None else penetration_length_mm,
        'nominal_steps_per_minimum_penetration_length': None if leakage is None else penetration_length_mm / spacing,
        'elapsed_seconds': perf_counter() - started,
        'hostname': socket.gethostname(),
        'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
    }
    pressure_scale = VISCOSITY * VELOCITY_SCALE / HEIGHT_M
    np.savez_compressed(output / 'fields.npz', triangles=fem.triangles, velocity_dofs=fem.velocity_dofs,
                        pressure_points_mm=fem.points * REFERENCE_HEIGHT_MM,
                        velocity_points_mm=fem.velocity_points * REFERENCE_HEIGHT_MM,
                        pressure_pa=values[2 * fem.n_u:] * pressure_scale,
                        velocity_m_per_s=values[:2 * fem.n_u].reshape(2, -1).T * VELOCITY_SCALE)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage_dir', type=Path)
    parser.add_argument('smoothed_contour', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('case', choices=CASES)
    parser.add_argument('spacing', type=float, choices=(*LEVELS, .03125))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    print(f'Starting {args.case} on target h={args.spacing:g} mm', flush=True)
    try:
        result = run_case(args.stage_dir, args.smoothed_contour, args.output, args.case, args.spacing)
    except Exception as error:
        result = {'case': args.case, 'mesh_target_mm': args.spacing, 'state_passed': False,
                  'failure': f'{type(error).__name__}: {error}', 'slurm_job_id': os.environ.get('SLURM_JOB_ID')}
        (args.output / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
        raise
    (args.output / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result, allow_nan=False), flush=True)
    raise SystemExit(0 if result['state_passed'] else 1)


if __name__ == '__main__':
    main()
