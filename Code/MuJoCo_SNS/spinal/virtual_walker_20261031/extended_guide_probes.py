import json
from run_virtual_walker import run
from import_virtual_walker import HERE
from assess_guides import analyze
rows=[]
for scale in (.5,6.):
 for clearance in (.01,.03):
  tag=f'guides_v1_t{scale:g}_c{clearance}_h0_p5';folder=HERE/'runs/W2L_VirtualWalker'/tag;params=dict(neural_time_scale=scale,guide_clearance=clearance,sagittal_harness=0.,rail_height_offset=-.02,yaw_perturb_nm=5.)
  result=run(duration=10.,tag=tag,**params);assessment=analyze(folder);(folder/'PARAMETERS.json').write_text(json.dumps(dict(contact_si=False,pf_knee_flex_scale=1.,ankle_pf_crossed=False,symmetric_rg_ext_tonic=None,**params),indent=2));rows.append(dict(name=tag,result=result,assessment=assessment));(HERE/'runs/GUIDE_EXTREMES_V1.json').write_text(json.dumps(rows,indent=2))
rows=json.loads((HERE/'runs/GUIDE_STUDY_V1.json').read_text())
for row in rows:row['assessment']=analyze(row['assessment']['source'])
(HERE/'runs/GUIDE_STUDY_V1.json').write_text(json.dumps(rows,indent=2))
