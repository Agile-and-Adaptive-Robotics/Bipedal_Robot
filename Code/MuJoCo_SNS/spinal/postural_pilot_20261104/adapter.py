"""Persistent scientific worker for a small mirrored premotor neural fit."""
from pathlib import Path
import sys,time
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run import trial

def handle_request(request):
    if request['op']=='preflight':
        row,_,_=trial('premotor','inhibitory',80.,10.)
        return dict(ok=True,ready_for_optimization=row['accepted'],baseline=row)
    if request['op']!='evaluate':raise ValueError('Unknown operation')
    parameters=request['params'];rows=[]
    # Training conditions contain both perturbation directions. Larger starts,
    # larger pushes, and neural noise are reserved for held-out evaluation.
    for initial,push in ((.04,8.),(-.04,-8.)):
        row,_,_=trial('premotor','inhibitory',80.,10.,initial=initial,push=push,parameters=parameters)
        rows.append(row)
    value=sum(1000.*(not r['accepted'])+r['final_rms_deg']+ .2*r['maximum_angle_deg']+.2*(r['settling_after_push_s']or 8.)+.1*r['mean_squared_excitation']for r in rows)/len(rows)
    return dict(ok=True,value=value,diagnostics=dict(accepted=all(r['accepted']for r in rows),conditions=rows,scope='Isolated fixed-base ankle bench; not a full walker'))
