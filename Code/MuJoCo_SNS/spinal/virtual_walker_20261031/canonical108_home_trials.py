"""Same-source108 home posture and explicit four-edge ankle-phase comparison."""
import json
import numpy as np
from import_virtual_walker import HERE
from run_virtual_walker import run
from assess_trial import assess
pose=HERE.parent/'physics_fidelity/runs/canonical_12II/ik_reference/flat_existing_rails.npz'
z=np.load(pose);q=z['qpos'];rows=[]
for crossed in (False,True):
 tag='canonical108_ground_nearflat_'+('ankle_crossed' if crossed else 'source_phase')+'_v1';folder=HERE/'runs/W2L_VirtualWalker_Newer'/tag
 result=run(newer=True,duration=10.,tag=tag,initial_qpos=q,ankle_pf_crossed=crossed);a=assess(folder)
 out=dict(pose_source=str(pose),pose_key='qpos',ankle_pf_crossed=crossed,result=result,assessment=a,scope='Same canonical108 source/limits/rails/guides and SNS/Hill initialization. Only actual initial physical qpos and declared four PF→ankleMN target swaps differ. No force-target allocator or static-tension preload.')
 (folder/'HOME_COMPARISON.json').write_text(json.dumps(out,indent=2));rows.append(out);(HERE/'runs/CANONICAL108_HOME_COMPARISON.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(tag=tag,swings=[a['legs'][s]['unloaded_clear_swings'] for s in ('L','R')],alternation=a['alternating_transition_fraction'],support=result['mean_net_vertical_assistance_BW'],failure=result['failure']),indent=2),flush=True)
