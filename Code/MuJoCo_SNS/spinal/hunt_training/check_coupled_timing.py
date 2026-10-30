"""Check actual full-network RG phase against the reference used for PF fitting."""
from pathlib import Path
import json,numpy as np
H=Path(__file__).resolve().parent;O=H/'runs/feasible_hunt_20261030_v1';rows={}
for mode,T,duty in [('afferented',1.233,.61),('deafferented',2.466,.5)]:
    z=np.load(O/('pf_'+mode)/'pf_teacher.npz');names=z['neuron_names'].tolist();t=z['t'];v=z['voltage'];keep=t>4*T;by_side={}
    for side,shift in [('r',0.),('l',.5)]:
        contrast=v[:,names.index('RG_E_'+side)]-v[:,names.index('RG_F_'+side)]
        prediction=1./(1.+np.exp(-np.clip(contrast,-30,30)));desired=(t/T+shift)%1<duty
        onset=np.flatnonzero((contrast[1:]>0)&(contrast[:-1]<=0))+1;onset=onset[t[onset]>4*T]
        phase=(t[onset]/T+shift)%1
        by_side[side]=dict(period_s=float(np.median(np.diff(t[onset]))) if len(onset)>1 else None,
            phase_error=float(np.mean(np.minimum(phase,1.-phase))),stance_mse=float(np.mean((prediction[keep]-desired[keep])**2)),cycles=len(onset))
    rows[mode]=by_side
(O/'coupled_timing_check.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
