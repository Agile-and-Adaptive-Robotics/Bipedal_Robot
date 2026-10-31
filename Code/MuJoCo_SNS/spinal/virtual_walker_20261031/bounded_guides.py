"""Same source controller, bounded cadence×physical-guide clearance study."""
import json
from run_virtual_walker import run
from import_virtual_walker import HERE
from assess_guides import analyze
rows=[]
cases=[(s,c,0.,5.) for s in (1.,2.,4.) for c in (.01,.03)]+[(2.,None,0.,5.),(2.,None,1.,5.),(2.,.01,0.,0.)]
for scale,clearance,harness,pulse in cases:
 tag=f'guides_v1_t{scale:g}_c{clearance}_h{harness:g}_p{pulse:g}';folder=HERE/'runs/W2L_VirtualWalker'/tag
 params=dict(neural_time_scale=scale,guide_clearance=clearance,sagittal_harness=harness,rail_height_offset=-.02,yaw_perturb_nm=pulse)
 result=run(duration=10.,tag=tag,**params);assessment=analyze(folder);(folder/'PARAMETERS.json').write_text(json.dumps(dict(contact_si=False,pf_knee_flex_scale=1.,ankle_pf_crossed=False,symmetric_rg_ext_tonic=None,**params),indent=2));rows.append(dict(name=tag,result=result,assessment=assessment));(HERE/'runs/GUIDE_STUDY_V1.json').write_text(json.dumps(rows,indent=2))
