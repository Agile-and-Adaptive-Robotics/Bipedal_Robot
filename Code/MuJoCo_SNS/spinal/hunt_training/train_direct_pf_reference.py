"""Reference symmetric driven PF circuits, with Hunt A/B fixed and C/D refitted."""
import os,sys,json,shutil,runpy,subprocess
from pathlib import Path
H=Path(__file__).resolve().parent;source=H/'runs/native_joint_pf_mono_20261030_v1/weights';name='direct_pf_reference_20261030_s1_v1';out=H/'runs'/name
if out.exists():raise SystemExit('Existing reference candidate retained')
weights=out/'weights';weights.mkdir(parents=True)
for filename in ('pf_weights.json','timing_weights.json'):shutil.copyfile(source/filename,weights/filename)
z=json.loads((weights/'pf_weights.json').read_text());changes={}
for key,value in z['edges'].items():
    pre,post=key.split(' -> ')
    if post.startswith(('PF_HIP-','PF_KNEE-')) and pre.startswith('HUNT_BANK'):
        z['edges'][key]=0.;changes[key]=0.
    elif pre.startswith(('PF_IN_HIP-','PF_IN_KNEE-')):
        z['edges'][key]=2.749;changes[key]=2.749
z['scope']='Symmetric driven four-cell hip/knee PF reference; bank inputs to those pairs disabled. Other functional PF weights retained.'
(weights/'pf_weights.json').write_text(json.dumps(z,indent=2))
(weights/'circuit_options.json').write_text(json.dumps(dict(joint_pf_inhibition_scale=1.,direct_joint_pf_gain=2.4),indent=2))
(out/'PLAN.json').write_text(json.dumps(dict(changes=changes,scope='Graded driven PF cells with laminated inhibition, not intrinsic NaP PF oscillators; original w2lvar starting conductances. No physical pass assumed.'),indent=2))
shutil.copyfile(Path(__file__),out/Path(__file__).name)
os.environ.update(HUNT_FEASIBLE_CAMPAIGN=name,HUNT_FEASIBLE_TARGET='mono_command',HUNT_FEASIBLE_MODES='afferented')
for key in ('HUNT_FEASIBLE_PF_ONLY','HUNT_JOINT_PF_INHIBITION','HUNT_TIMING_SOURCE','HUNT_PF_INACTIVE'):os.environ.pop(key,None)
runpy.run_path(str(H/'train_feasible_hunt.py'),run_name='__main__')
with (out/'air_pair.log').open('w') as log:
    code=subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(weights),'--output-directory',str(out/'air_pair')],stdout=log,stderr=subprocess.STDOUT).returncode
if code:raise SystemExit(code)
print('DIRECT_PF_REFERENCE_COMPLETE',flush=True)
