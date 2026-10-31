import json
import numpy as np
from scipy.signal import find_peaks
from import_virtual_walker import HERE
rows=json.loads((HERE/'runs/CALIBRATION_V1.json').read_text())
for row in rows:
 path=row['assessment']['source'];z=np.load(path+'/trajectory.npz');names=z['neuron_names'].tolist();t=z['t'];mask=t>=1.;dt=row['result']['physics_dt'];metrics={}
 for side in ('L','R'):
  for level,name in [('RG',side+' RG ext'),('PF',side+' Knee PF flx'),('MN',side+' Knee MN flx')]:
   if name in names:
    y=z['neuro'][mask,names.index(name)];pk,_=find_peaks(y,prominence=1,distance=round(.15/dt));metrics[side+level]=dict(range=float(np.ptp(y)),peaks=len(pk),median_period_s=float(np.median(np.diff(t[mask][pk]))) if len(pk)>1 else None)
 if 'R RG ext' in names:metrics['RG_LR_correlation']=float(np.corrcoef(z['neuro'][mask,names.index('L RG ext')],z['neuro'][mask,names.index('R RG ext')])[0,1])
 metrics['knee_LR_correlation']=float(np.corrcoef(z['angles'][mask,1],z['angles'][mask,4])[0,1]);metrics['knee_cycle_counts']=[row['assessment']['legs'][s]['knee_peak_count'] for s in ('L','R')]
 row['phase_diagnostics']=metrics;print(row['name'],json.dumps(metrics))
(HERE/'runs/PHASE_DIAGNOSTICS_V1.json').write_text(json.dumps(rows,indent=2))
