"""Standing-pose correction and isolated task-feedback comparisons."""
import sys,json
from pathlib import Path
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/balance_20261029_v2';campaign.O.mkdir(exist_ok=True)
run=campaign.run;base=dict(stage=3,restriction=1.,height_offset=-.02,start_hip_deg=-7.,start_knee_deg=-10.,start_ankle_deg=0.)
records=[]
for name,options in [('inherited',None),('disabled',dict(kp=0,kd=0)),('pd_low',dict(kp=6,kd=2)),('pd_high',dict(kp=20,kd=4)),('pda',dict(kp=6,kd=2,ka=.05))]:
    cfg={**base}
    if options is not None:cfg['balance']={**options,'delay_s':.1}
    r=run(name,cfg)
    if r:records.append((campaign.objective(r),cfg,r))
if not records:raise RuntimeError('No completed neutral-stance tests')
# Record the best active controller separately from both control conditions.
active=[row for row in records if row[1].get('balance',{}).get('kp',0)>0.]
score,winner,r=min(active,key=lambda x:x[0]);(campaign.O/'winner.json').write_text(json.dumps(dict(objective=score,configuration=winner,assessment=r['assessment'],measured=r['measured']),indent=2))
for backend in ('sparse','native'):run('benchmark_'+backend,{**winner,'backend':backend})
for name,extra in [('heldout_forward',dict(push_N=30.,push_time=4.)),('heldout_lateral',dict(lateral_push_N=30.,push_time=4.)),('less_restraint',dict(restriction=.5))]:run(name,{**winner,**extra})
for stage in (1,2,4,5):
    cfg={**winner,'stage':stage}
    if stage in (1,2):cfg.pop('balance',None);cfg['restriction']=0.
    if stage in (4,5):cfg['motor_scales']={'knee_ext_l':.5,'knee_ext_r':.5,'knee_flex_l':2.,'knee_flex_r':2.}
    if stage==5:cfg['contact_pf_scale']=.25
    run('curriculum_stage_'+str(stage),cfg)
