"""Nested bounded optimization of a fitted baffle and fitted branch shoulders."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from .project43_explicit import axis, solve


def run(output):
    output.mkdir()
    history = []
    cache = {}
    def evaluate(h,s):
        key = float(h),float(s)
        if key not in cache:
            assert .15<=h<=4 and 0<=s<=1 and h+4*s<=4+1e-10
            xa = axis([0,9.5-s,9.5,10.5,10.5+s,13,14,20],.25)
            ya = axis([0,h/4,h/2,.75*h,h,4,6,7],.25)
            r = solve(xa,ya,h,s)
            np.testing.assert_allclose(r['area_mm2'],79+h+4*s,atol=1e-8)
            assert r['area_mm2']<=83+1e-8
            r.update(open_height_mm=h,shoulder_width_mm=s)
            cache[key] = r
            history.append(r)
            temp = output/'history.tmp'
            temp.write_text(json.dumps(history,indent=2,allow_nan=False))
            temp.replace(output/'history.json')
        return cache[key]
    baseline = evaluate(4.,0.)
    def feasible_at(s,target):
        lo,hi = .15,4-4*s
        def discrepancy(h):
            return evaluate(h,s)['side_share']-target
        if discrepancy(lo)*discrepancy(hi)>0:
            return None
        height = brentq(discrepancy,lo,hi,xtol=5e-6)
        return evaluate(height,s)
    start = feasible_at(0.,.5)
    if start is None:
        raise RuntimeError('Center split target unbracketed for zero shoulders; inspect history.')
    # Search shoulder width along each band edge and center. Margin protects
    # reported feasibility against numerical root-finder tolerances.
    def objective(s,target):
        design = feasible_at(float(s),target)
        return design['pumping_power_W_per_m'] if design else 1e6*start['pumping_power_W_per_m']
    for target in (.49901,.5,.50099):
        objective(0.,target)
        minimize_scalar(objective,args=(target,),bounds=(0.,.96),method='bounded',
                        options=dict(xatol=.01,maxiter=18))
    best = min((r for r in history if .499<=r['side_share']<=.501),
               key=lambda r:r['pumping_power_W_per_m'])
    report = dict(original=baseline,starting_feasible=start,best_feasible=best,
                  improvement_fraction=1-best['pumping_power_W_per_m']/start['pumping_power_W_per_m'],
                  evaluations=len(history),mesh_spacing_mm=.25,
                  search='Bounded shoulder search at three split targets with bracketed baffle adjustment; best observed feasible geometry.')
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False))


if __name__=='__main__':
    run(Path(sys.argv[1]))
