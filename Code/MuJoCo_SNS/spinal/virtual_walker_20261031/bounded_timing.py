from source_inputs import historical_source
"""Timing-first mechanical-scale hypothesis, fixed rails and sagittal guidance."""
import json
from run_virtual_walker import run
from assess_trial import assess
from import_virtual_walker import HERE,REPO
rows=[]
for scale in (1.,2.,3.,4.):
 tag='timing_v1_scale_'+str(int(scale));folder=HERE/'runs/W2L_VirtualWalker'/tag
 result=run(duration=10.,tag=tag,sagittal_harness=1.,rail_height_offset=-.02,neural_time_scale=scale);measures=assess(folder);params=dict(contact_si=False,pf_knee_flex_scale=1.,ankle_pf_crossed=False,symmetric_rg_ext_tonic=None,sagittal_harness=1.,rail_height_offset=-.02,neural_time_scale=scale);(folder/'PARAMETERS.json').write_text(json.dumps(params,indent=2));rows.append(dict(scale=scale,result=result,assessment=measures));(HERE/'runs/TIMING_V1.json').write_text(json.dumps(rows,indent=2))
source=str(historical_source('W2L_VirtualWalker_Newer_SIContact'))
tag='gui204_sagittal_drive0_lower20';result=run(duration=10.,source=source,tag=tag,sagittal_harness=1.,rail_height_offset=-.02,symmetric_rg_ext_tonic=0.);folder=HERE/'runs/W2L_VirtualWalker_Newer_SIContact'/tag;measures=assess(folder);rows.append(dict(scale=1.,result=result,assessment=measures));(HERE/'runs/TIMING_V1.json').write_text(json.dumps(rows,indent=2))
