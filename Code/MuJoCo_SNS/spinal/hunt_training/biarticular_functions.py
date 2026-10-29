"""Function signs from actual sampled muscle moment arms, not muscle names."""
import json
from pathlib import Path
import numpy as np
H=Path(__file__).resolve().parent;z=np.load(H/'runs/targets/targets.npz',allow_pickle=True)
acts=z['act_names'].tolist();joints=z['joints'].tolist();rho=z['rho']
report={}
for base in ('rect_fem','bifemlh','semimem','semiten','grac','sar','tfl','med_gas','lat_gas'):
    act=base+'_r'
    if act not in acts:continue
    cols={}
    for coordinate,sign in [('hip_flexion_r',-1.),('knee_angle_r',-1.),('hip_adduction_r',-1.),('hip_rotation_r',1.)]:
        values=rho[:,joints.index(coordinate),acts.index(act)]*sign*1000.
        cols[coordinate]=dict(min_mm=float(values.min()),median_mm=float(np.median(values)),max_mm=float(values.max()),positive_fraction=float(np.mean(values>1.)))
    report[base]=cols
(H/'runs/biarticular_function_signs.json').write_text(json.dumps(report,indent=2))
print('Positive columns: hip flexion, knee flexion, hip adduction, native/OpenSim hip_rotation coordinate. mm moment arms, sampled experimental gait poses only.')
for muscle,columns in report.items():print(muscle,{key:round(value['median_mm'],2) for key,value in columns.items()})
