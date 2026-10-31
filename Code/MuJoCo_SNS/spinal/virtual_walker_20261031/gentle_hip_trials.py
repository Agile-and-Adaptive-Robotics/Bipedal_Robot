"""Three declared mirrored hip-only recruitment reductions; baseline unchanged."""
import hashlib,json
from import_virtual_walker import HERE
from run_virtual_walker import run
from assess_trial import assess
path=HERE/'RUSTUNA108_GROUND_WINNER.json';original=path.read_bytes();cfg=json.loads(original);expected=cfg.pop('expected_source_sha256');assert hashlib.sha256(open(cfg['source'],'rb').read()).hexdigest()==expected
rows=[]
for factor in (.8,.6,.4):
 tag='canonical108_ground_gentle_hip_'+str(factor).replace('.','p')+'_v1';r=run(newer=True,duration=10.,tag=tag,hip_pf_mn_scale=factor,**cfg);folder=HERE/'runs/W2L_VirtualWalker_Newer'/tag;a=assess(folder);out=dict(tag=tag,hip_factor=factor,result=r,assessment=a);rows.append(out);(HERE/'runs/GENTLE_HIP_TRIALS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(factor=factor,gate=a['alternating_supported_stepping_gate'],swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],alternation=a['alternating_transition_fraction'],support=r['mean_net_vertical_assistance_BW'],failure=r['failure']),indent=2),flush=True)
assert path.read_bytes()==original,'Original winner must remain unchanged'
