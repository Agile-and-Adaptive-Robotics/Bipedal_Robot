"""Fresh winner repeat, longer run, then contact-only controls from corrected home."""
import json,subprocess,sys
import numpy as np
from import_virtual_walker import HERE
from run_virtual_walker import run
from assess_trial import assess
pose=HERE.parent/'physics_fidelity/runs/canonical_12II/ik_reference/flat_existing_rails.npz'
cfg=dict(initial_qpos=np.load(pose)['qpos'].tolist(),ankle_pf_crossed=True)
config=HERE/'CANONICAL108_NEARFLAT_ANKLE_PHASE.json';config.write_text(json.dumps(cfg,indent=2))
base=HERE/'runs/W2L_VirtualWalker_Newer';rows=[]
tag='canonical108_ground_nearflat_ankle_crossed_fresh_v1'
with (HERE/'runs/canonical108_home_fresh.log').open('w') as log:
 subprocess.run([sys.executable,str(HERE/'run_virtual_walker.py'),'--newer','--config',str(config),'--duration','10','--tag',tag],stdout=log,stderr=subprocess.STDOUT,check=True)
folder=base/tag;a=assess(folder);reference=np.load(base/'canonical108_ground_nearflat_ankle_crossed_v1/trajectory.npz');fresh=np.load(folder/'trajectory.npz')
parity={k:bool(reference[k].shape==fresh[k].shape and reference[k].dtype==fresh[k].dtype and reference[k].tobytes()==fresh[k].tobytes()) for k in reference.files if reference[k].dtype.kind in 'fiu'};assert all(parity.values()),parity
(folder/'FRESH_PARITY.json').write_text(json.dumps(dict(all_numeric_arrays_bit_exact=True,arrays=parity),indent=2))
rows.append(dict(tag=tag,result=json.loads((folder/'RESULT.json').read_text()),assessment=a))
cases=[('canonical108_ground_nearflat_ankle_crossed_20s_v1',dict(duration=20.,belt_motor_end_s=20.)),('canonical108_ground_nearflat_crossed_contact_RG_off_v1',dict(contact_rg_scale=0.)),('canonical108_ground_nearflat_crossed_contact_PF_off_v1',dict(contact_pf_scale=0.)),('canonical108_ground_nearflat_crossed_contact_MN_off_v1',dict(contact_mn_scale=0.)),('canonical108_ground_nearflat_crossed_all_contact_outputs_off_v1',dict(contact_rg_scale=0.,contact_pf_scale=0.,contact_mn_scale=0.))]
for tag,extra in cases:
 args=dict(newer=True,duration=10.,tag=tag,**cfg);args.update(extra);r=run(**args);folder=base/tag;a=assess(folder);out=dict(tag=tag,result=r,assessment=a);rows.append(out)
 (folder/'HOME_VALIDATION.json').write_text(json.dumps(out,indent=2));(HERE/'runs/CANONICAL108_HOME_VALIDATION.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(tag=tag,swings=[a['legs'][s]['unloaded_clear_swings'] for s in ('L','R')],alternation=a['alternating_transition_fraction'],gate=a['alternating_supported_stepping_gate'],support=r['mean_net_vertical_assistance_BW'],failure=r['failure']),indent=2),flush=True)
