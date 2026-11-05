"""Persistent robust postural fit: task gains and mirrored neural gains only."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from run import trial
TRAIN_NOISE=.025

def evaluate(parameters,held_out=False):
    neural={k:v for k,v in parameters.items()if k not in('kp','kd')};rows=[]
    if held_out:
        cases=[dict(initial=sign*.08,push=-sign*16.,noise=.05,delay_s=delay,seed=seed)for delay in(0.,.03,.07)for seed in(61210,61211,61212)for sign in(1.,-1.)]
    else:
        cases=[dict(initial=sign*.04,push=sign*8.,noise=0.,delay_s=0.,seed=61200)for sign in(1.,-1.)]
        cases += [dict(initial=sign*.08,push=-sign*16.,noise=TRAIN_NOISE,delay_s=delay,seed=61200 if sign>0 else 61201)for delay in(.03,.07)for sign in(1.,-1.)]
    for case in cases:
        row,_,_=trial('premotor','inhibitory',parameters['kp'],parameters['kd'],parameters=neural,backend='native',**case);rows.append(row)
    costs=[1000.*(not r['accepted'])+min(r['final_rms_deg'],100.)+.2*r['maximum_angle_deg']+.2*(r['settling_after_push_s']or 8.)+.1*r['mean_squared_excitation']for r in rows]
    return dict(ok=True,value=float(np.mean(costs)+.25*np.max(costs)),diagnostics=dict(accepted=all(r['accepted']for r in rows),passed=sum(r['accepted']for r in rows),total=len(rows),conditions=rows,scope='Robust isolated ankle-bench test; uniform added sensor delay, shared antagonistic neural scales, unchanged mechanical model'))

def handle_request(request):
    global TRAIN_NOISE
    if request['op']=='preflight':
        TRAIN_NOISE=float(request.get('training_noise_nA',.025))
        result=evaluate(dict(kp=80.,kd=10.,reciprocal_scale=1.,Ia_scale=1.,II_scale=1.,Ib_scale=1.))
        return dict(ok=True,ready_for_optimization=all(np.isfinite(r['final_rms_deg'])for r in result['diagnostics']['conditions']),baseline=result)
    if request['op']not in('evaluate','held_out'):raise ValueError('Unknown request')
    return evaluate(request['params'],request['op']=='held_out')
