from source_inputs import historical_source
"""Small, hypothesis-driven bilateral neural calibration; not a completed Hunt fit."""
import json
from pathlib import Path
from run_virtual_walker import run
from assess_trial import assess
from import_virtual_walker import HERE,REPO
source=str(historical_source('W2L_VirtualWalker_Newer_SIContact'))
cases=[('original_baseline',None,dict()),('original_crossed_ankle',None,dict(ankle_pf_crossed=True)),('original_knee1p1',None,dict(pf_knee_flex_scale=1.1)),('original_crossed_knee1p1',None,dict(ankle_pf_crossed=True,pf_knee_flex_scale=1.1)),('gui204_sagittal',source,dict(sagittal_harness=1.,symmetric_rg_ext_tonic=2.)),('gui204_sagittal_crossed',source,dict(sagittal_harness=1.,symmetric_rg_ext_tonic=2.,ankle_pf_crossed=True)),('gui204_sagittal_crossed_lower20',source,dict(sagittal_harness=1.,symmetric_rg_ext_tonic=2.,ankle_pf_crossed=True,rail_height_offset=-.02)),('gui204_sagittal_lower20',source,dict(sagittal_harness=1.,symmetric_rg_ext_tonic=2.,rail_height_offset=-.02))]
rows=[]
for name,path,params in cases:
 tag='calibration_v1_'+name;folder=HERE/'runs'/('W2L_VirtualWalker' if path is None else 'W2L_VirtualWalker_Newer_SIContact')/tag
 if (folder/'RESULT.json').exists():result=json.loads((folder/'RESULT.json').read_text())
 else:result=run(duration=8.,source=path,tag=tag,**params)
 assessed=assess(folder);config=dict(contact_si=False,pf_knee_flex_scale=1.,ankle_pf_crossed=False,symmetric_rg_ext_tonic=None,sagittal_harness=0.,rail_height_offset=0.);config.update(params);(folder/'PARAMETERS.json').write_text(json.dumps(config,indent=2))
 rows.append(dict(name=name,source=path,parameters=config,result=result,assessment=assessed))
 (HERE/'runs/CALIBRATION_V1.json').write_text(json.dumps(rows,indent=2))
print('COMPLETED',len(rows),'bounded source-controller cases')
