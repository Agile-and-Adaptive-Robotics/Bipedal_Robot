from source_inputs import historical_source
"""Body stability versus coordination, preserving GUI204 source tonic2/0nA."""
import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from run_virtual_walker import run
from import_virtual_walker import HERE,REPO
from assess_guides import analyze
rows=[]
for candidate,stem in [('GUI1nA','W2L_VirtualWalker_Newer_GUI1nA'),('SIContact20nA','W2L_VirtualWalker_Newer_SIContact')]:
 source=historical_source(stem)
 for clearance in (None,.03):
  tag=f'newer_guide_pair_v1_t2_c{clearance}';folder=HERE/'runs'/stem/tag;params=dict(neural_time_scale=2.,guide_clearance=clearance,sagittal_harness=0.,rail_height_offset=-.02,yaw_perturb_nm=0.,symmetric_rg_ext_tonic=None)
  result=run(source=source,duration=10.,tag=tag,**params);assessment=analyze(folder);z=np.load(folder/'trajectory.npz');names=z['neuron_names'].tolist();mask=z['t']>=1.;neural={}
  for name in ('L RG ext','L RG flx','R RG ext','R RG flx'):
   v=z['neuro'][mask,names.index(name)];peaks,_=find_peaks(v,prominence=1.,distance=round(.15/result['physics_dt']));neural[name]=dict(range_mV=float(np.ptp(v)),peaks=int(len(peaks)))
  assessment.update(candidate=candidate,source_left_right_RGext_tonic_nA=[2.,0.],neural_rg=neural,body_height_stability_gate=result['failure'] is None and result['root_y_min_m']>=.6,root_y_min_m=result['root_y_min_m'],interpretation='Body-height retention is not a neural stepping pass; source graph96/204 and tonic asymmetry preserved. No perturbation or new parameter fitting.');(folder/'GUIDE_ASSESSMENT.json').write_text(json.dumps(assessment,indent=2));(folder/'PARAMETERS.json').write_text(json.dumps(params,indent=2));rows.append(dict(name=f'{candidate}_c{clearance}',result=result,assessment=assessment));(HERE/'runs/NEWER_GUIDE_PAIRS_V1.json').write_text(json.dumps(rows,indent=2))
