"""Compare expected-input neural timing with live-body voltages, without promotion."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('campaign',type=Path)
root=parser.parse_args().campaign.resolve()

def summarize(z,physical=False):
    names=z['names' if physical else 'neuron_names'].tolist()
    t=z['t'];v=z['voltage'];keep=t>=5.;out={}
    for side in ('r','l'):
        cells=['RG_E_'+side,'RG_F_'+side,'PF_KNEE-E_'+side,'PF_KNEE-F_'+side,
               'PF_BI_JOINT_RECT_'+side,'PF_BI_JOINT_HAM_'+side,'MN_rect_fem_'+side,
               'MN_bifemsh_'+side,'MN_vas_int_'+side]
        rows={}
        for cell in cells:
            y=v[keep,names.index(cell)];tt=t[keep]
            peaks,_=find_peaks(y,prominence=.1,distance=max(1,int(.2/np.median(np.diff(t)))))
            rows[cell]=dict(minimum=float(y.min()),maximum=float(y.max()),mean=float(y.mean()),
                            peaks=int(len(peaks)),period_s=float(np.median(np.diff(tt[peaks]))) if len(peaks)>1 else None)
        out[side]=rows
    return out

report={}
for mode in ('afferented','deafferented'):
    path=root/('pf_'+mode)/'pf_teacher.npz'
    if path.exists():report['expected_'+mode]=summarize(np.load(path))
path=root/'air_retry/stage2.mechanics.npz'
if path.exists():report['physical_stage2']=summarize(np.load(path),True)
for mode in ('afferented','deafferented'):
    teacher=np.load(root/('pf_'+mode)/'pf_teacher.npz')
    validation=np.load(root/'weights'/('step_C_'+mode+'.npz'))
    cells=teacher['pf_cells'].tolist();keep=validation['t']>=5.
    report['free_pf_'+mode]={cell:dict(minimum=float(validation['pf'][keep,cells.index(cell)].min()),
        maximum=float(validation['pf'][keep,cells.index(cell)].max())) for cell in cells if cell.startswith(('PF_KNEE','PF_HIP'))}
    target=validation['target'][keep];actual=validation['activation'][keep]
    acts=np.load(root.parent/'targets/targets.npz',allow_pickle=True)['act_names'].tolist()
    order=np.argsort(-target.max(axis=0))[:12]
    report['motor_fit_'+mode]=[dict(muscle=acts[i],target_peak=float(target[:,i].max()),
        predicted_min=float(actual[:,i].min()),predicted_max=float(actual[:,i].max()),
        r_squared=float(1.-np.mean((actual[:,i]-target[:,i])**2)/max(float(np.var(target[:,i])),1e-12))) for i in order]
(root/'air_failure_trace.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
