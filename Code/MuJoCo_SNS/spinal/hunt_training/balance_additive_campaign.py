"""Add task feedback while retaining the torso/posture stabilization loop."""
import json
from pathlib import Path
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/balance_20261029_v3';campaign.O.mkdir(exist_ok=True)
run=campaign.run;base=dict(stage=3,restriction=1.,height_offset=-.02,start_hip_deg=-7.,start_knee_deg=-10.,start_ankle_deg=0.)
records=[]
for name,kp,kd,ka in [('disabled',0.,0.,0.),('low',.5,.2,0.),('medium',2.,.5,0.),('high',6.,2.,0.),('acceleration',2.,.5,.03)]:
    cfg={**base,'balance':dict(kp=kp,kd=kd,ka=ka,delay_s=.1,replace_inherited=False)}
    r=run(name,cfg)
    if r:records.append((campaign.objective(r),cfg,r))
score,winner,r=min([row for row in records if row[1]['balance']['kp']>0],key=lambda x:x[0])
(campaign.O/'winner.json').write_text(json.dumps(dict(objective=score,configuration=winner,assessment=r['assessment'],measured=r['measured']),indent=2))
for backend in ('sparse','native'):run('benchmark_'+backend,{**winner,'backend':backend})
for name,extra in [('heldout_forward',dict(push_N=30.,push_time=4.)),('heldout_lateral',dict(lateral_push_N=30.,push_time=4.)),('less_restraint',dict(restriction=.5))]:run(name,{**winner,**extra})
for stage in (1,2,4,5):
    cfg={**winner,'stage':stage}
    if stage in (1,2):
        cfg.pop('balance',None);cfg['restriction']=0.;cfg['height_offset']=0.
        for key in ('start_hip_deg','start_knee_deg','start_ankle_deg'):cfg.pop(key,None)
    if stage in (4,5):cfg['motor_scales']={'knee_ext_l':.5,'knee_ext_r':.5,'knee_flex_l':2.,'knee_flex_r':2.}
    if stage==5:cfg['contact_pf_scale']=.25
    run('curriculum_stage_'+str(stage),cfg)
