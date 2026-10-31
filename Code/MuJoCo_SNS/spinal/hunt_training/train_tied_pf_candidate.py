"""Isolated tied-inhibition Hunt C/D candidate; A/B and four-cell topology fixed."""
import os,sys,json,shutil,runpy,hashlib
from pathlib import Path
H=Path(__file__).resolve().parent;name='tied_native_pf_20261030_s1_v1';out=H/'runs'/name
if out.exists():raise SystemExit('Existing campaign retained')
os.environ.update(HUNT_NATIVE_PF_CAMPAIGN=name,HUNT_NATIVE_TARGET='mono_command',HUNT_NATIVE_SLOW_WEIGHT='1',HUNT_NATIVE_PF_TIE_INHIBITION='1')
runpy.run_path(str(H/'fit_native_joint_pf.py'),run_name='__main__')
os.environ.update(HUNT_FEASIBLE_CAMPAIGN=name,HUNT_FEASIBLE_TARGET='mono_command',HUNT_FEASIBLE_MODES='afferented')
for key in ('HUNT_FEASIBLE_PF_ONLY','HUNT_JOINT_PF_INHIBITION','HUNT_TIMING_SOURCE','HUNT_PF_INACTIVE'):os.environ.pop(key,None)
shutil.copyfile(Path(__file__),out/Path(__file__).name)
runpy.run_path(str(H/'train_feasible_hunt.py'),run_name='__main__')
with (out/'air_pair.log').open('w') as log:
    import subprocess
    result=subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(out/'weights'),'--output-directory',str(out/'air_pair')],stdout=log,stderr=subprocess.STDOUT)
if result.returncode:raise SystemExit(result.returncode)
print('TIED_CANDIDATE_AIR_COMPLETE',flush=True)
