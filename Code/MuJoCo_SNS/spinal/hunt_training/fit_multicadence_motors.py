"""Fit Hunt C motor weights at both cadences with the existing A/B/PF/BI frozen.

The normal-only motor fit can reproduce air motion while its slower,
deafferented motor pattern crosses. This isolated candidate reuses the same
native four-cell PF fit and teacher data, fits shared MN conductances at both
cadences, and then fits D at the normal cadence. Physical gates are separate.
"""
import os,sys,json,hashlib,shutil,argparse,runpy
from pathlib import Path

H=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--campaign',default='multicadence_mn_20261030_s1_v1')
a=p.parse_args()
if Path(a.campaign).name!=a.campaign:p.error('Use a campaign directory name, not a path')
source=H/'runs/native_joint_pf_mono_20261030_v1/weights'
out=H/'runs'/a.campaign
if out.exists():p.error('Choose a new campaign; existing results are retained')
weights=out/'weights';weights.mkdir(parents=True)
frozen=('pf_weights.json','bi_weights.json','timing_weights.json','circuit_options.json')
hashes={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in frozen}
for name in frozen:shutil.copyfile(source/name,weights/name)
plan=dict(machine='easteregg2',source=str(source),frozen_sha256=hashes,
    fit_modes=['afferented','deafferented'],target='mono_command',
    hypothesis='Shared MN conductances must fit the native PF outputs at both cadences; separate stage winners do not qualify',
    scope='Offline Hunt C/D fit only; physical stages require separate fresh-process validation')
(out/'EXPERIMENT.json').write_text(json.dumps(plan,indent=2))
shutil.copyfile(Path(__file__),out/Path(__file__).name)
os.environ.update(HUNT_FEASIBLE_CAMPAIGN=a.campaign,HUNT_FEASIBLE_MODES='afferented,deafferented',HUNT_FEASIBLE_TARGET='mono_command')
for key in ('HUNT_FEASIBLE_PF_ONLY','HUNT_JOINT_PF_INHIBITION','HUNT_TIMING_SOURCE','HUNT_PF_INACTIVE'):os.environ.pop(key,None)
runpy.run_path(str(H/'train_feasible_hunt.py'),run_name='__main__')
assert all(hashlib.sha256((weights/name).read_bytes()).hexdigest()==sha for name,sha in hashes.items()),'A/B/PF/BI weights changed during the motor fit'
(out/'FROZEN_INPUTS_VERIFIED.json').write_text(json.dumps(dict(accepted=True,sha256=hashes),indent=2))
