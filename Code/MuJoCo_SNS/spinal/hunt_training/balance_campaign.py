"""Resumeable isolated balance tuning; failed trials never graduate stages."""
import os,sys,json,subprocess
from pathlib import Path
H=Path(__file__).resolve().parent;O=H/'runs/balance_20261029_v1';O.mkdir(parents=True,exist_ok=True)
def run(name,cfg):
    p=O/(name+'.json');result=p.with_suffix('.result.json')
    if result.exists():return json.loads(result.read_text())
    p.write_text(json.dumps(cfg,indent=2))
    with p.with_suffix('.log').open('w') as f:r=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=f,stderr=subprocess.STDOUT)
    if r.returncode:
        p.with_suffix('.failure.json').write_text(json.dumps(dict(exit=r.returncode,accepted=False)))
        print(name,'FAILED',r.returncode,flush=True);print(p.with_suffix('.log').read_text()[-1200:],flush=True);return None
    value=json.loads(result.read_text());a=value['assessment'];print(name,'accepted',a['accepted'],'sway',a.get('sway_m'),'force',value['measured']['horizontal_harness_rms_bw'],flush=True);return value
def objective(r):
    a=r['assessment'];m=r['measured']
    if not a.get('finite') or not a.get('complete'):return 1e6
    return 300*a.get('sway_m',1.)+a.get('tilt_deg',100.)+100*m['horizontal_harness_rms_bw']+30*max(.5-sum(m['mean_foot_normal_bw']),0.)+10*m['mean_activation_squared']+100*int(not a.get('frontal_gate',False))
if __name__=='__main__':
    base=dict(stage=3,restriction=1.,height_offset=-.06)
    records=[]
    for i,(kp,kd,ka,delay,weight) in enumerate(((0,0,0,.1,1.),(2,1,0,.1,1.),(6,2,0,.1,1.),(12,3,0,.1,1.),(25,5,0,.1,1.),(6,2,.05,.1,1.),(12,3,0,.05,1.),(12,3,0,.1,.75))):
        cfg={**base,'balance':dict(kp=kp,kd=kd,ka=ka,delay_s=delay,somatic_weight=weight)}
        r=run(f'tune_{i:02}',cfg)
        if r:records.append((objective(r),cfg,r))
    if not records:raise RuntimeError('No completed physical balance tests')
    score,winner,r=min(records,key=lambda row:row[0]);(O/'winner.json').write_text(json.dumps(dict(objective=score,configuration=winner,assessment=r['assessment'],measured=r['measured']),indent=2))
    # Benchmark exactly the same fitted network, model, assistance, and logging.
    for backend in ('sparse','native'):
        run('benchmark_'+backend,{**winner,'backend':backend})
    for name,extra in [('holdout_forward',dict(push_N=30.,push_time=4.)),('holdout_lateral',dict(lateral_push_N=30.,push_time=4.)),('less_restraint',dict(restriction=.5)),('less_vertical_stiffness',dict(vertical_stiffness_scale=.25))]:
        run(name,{**winner,**extra})
    for stage in (1,2,4,5):
        cfg={**winner,'stage':stage}
        if stage in (1,2):cfg.pop('balance',None);cfg['restriction']=0.;cfg['height_offset']=0.
        if stage in (4,5):cfg['motor_scales']={'knee_ext_l':.5,'knee_ext_r':.5,'knee_flex_l':2.,'knee_flex_r':2.}
        if stage==5:cfg['contact_pf_scale']=.25
        run('curriculum_stage_'+str(stage),cfg)
