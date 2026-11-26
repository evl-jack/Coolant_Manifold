"""Regenerate web-facing data from frozen numerical evidence; no CFD runs."""
from pathlib import Path
import json
r=Path(__file__).resolve().parent
v=r/'evidence'/'verification'; n=r/'evidence'/'nonlinear'
def read(p): return json.loads(p.read_text())
mesh=read(v/'mesh-evidence/mesh_report.json'); nl=read(n/'nonlinear_comparison.json'); kkt=read(v/'kkt-results/fixed_pipe_kkt.json')
config=read(v/'input/configuration.json')
timing={c:read(v/f'timing-2408534/{c}/benchmark.json') for c in ['pilot-128','scale-32','scale-128','scale-512']}
summary={"schema_version":1,"scope":"archive-derived; planar per-unit-depth physics","primary_design":"raw beta=8 Navier-Stokes","configuration":config,"mesh":mesh,"nonlinear":nl,"kkt":kkt,"adjoint_initial":read(v/'adjoint-evidence/adjoint_merged.json'),"adjoint_followup":read(v/'adjoint-coordinate-supplement.json'),"timing":timing,"sources":{"mesh":"evidence/verification/mesh-evidence/mesh_report.json","nonlinear":"evidence/nonlinear/nonlinear_comparison.json","kkt":"evidence/verification/kkt-results/fixed_pipe_kkt.json","adjoint":"evidence/verification/adjoint-evidence/adjoint_merged.json and adjoint-coordinate-supplement.json","timing":"evidence/verification/timing-2408534/*/benchmark.json","configuration":"evidence/verification/input/configuration.json"},"missing":{"team":None,"public_repository":None,"deployed_site":None,"matched_area_flow_control":None,"graded_solver_bug":None,"validated_3d_cfd":None}}
# All converted units are explicit; preserve raw originals above.
summary['headline']={}
for key in ('fitted-original','fitted-raw','fitted-smoothed'):
 a=mesh['series'][key][-1]; summary['headline'][key]={'power_uW_per_m':a['mechanical_power_W_per_m']*1e6,'side_percent':a['side_fraction']*100,'resistance_pa_s_per_m2':a['hydraulic_resistance_pa_s_per_m2'],'area_mm2':a['physical_fluid_area_mm2'],'mesh_target_mm':a['mesh_target_mm']}
summary['headline']['raw_power_saving_percent']=mesh['paired_comparisons']['fitted-raw']['by_mesh'][-1]['power_saving_percent']
summary['headline']['adjoint_serial_speedup']=timing['pilot-128']['including_common_baseline']['measured_serial_fd_over_adjoint_speedup']
summary['headline']['kkt_stationarity']=kkt['results']['stages'][-1]['stationarity_inf']
(r/'data').mkdir(exist_ok=True);(r/'data/results.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
print('Wrote data/results.json from frozen solver results. No numerical rerun performed.')
