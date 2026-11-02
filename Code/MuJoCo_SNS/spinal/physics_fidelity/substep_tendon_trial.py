"""Isolated refinement: keep2ms neural clock, subdivide only source mechanics.

The saved2ms task-basis calibration is frozen. Every microstep retains the
constitutive/fiber and zero-hyperextension guards. No production/default edit.
"""
import argparse,sys,runpy,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main():
 p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--substeps',type=int,choices=[1,2,4],default=2);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 import thelen_dynamic as T
 init=T.ConstitutiveBackend.__init__;step=T.ConstitutiveBackend.step
 def refined_init(self,model,data,parameters,*args,**kwargs):
  if abs(float(model.opt.timestep)-.002)>1e-12:raise RuntimeError('Expected frozen2ms neural/physics base protocol')
  model.opt.timestep=.002/a.substeps;init(self,model,data,parameters,*args,**kwargs)
 def refined_step(self):
  for _ in range(a.substeps):step(self)
 T.ConstitutiveBackend.__init__=refined_init;T.ConstitutiveBackend.step=refined_step
 manifest=dict(scope='Integration refinement only; same2ms neural/measurement clock and frozen2ms task recruitment, source-compliant physics microsteps; every microstate checked for hyperextension and fiber domain.',neural_step_s=.002,physics_step_s=.002/a.substeps,substeps=a.substeps,driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),base_driver_sha256=hashlib.sha256((HERE/'live_tendon_trial.py').read_bytes()).hexdigest(),config=str(a.config.resolve()),config_sha256=hashlib.sha256(a.config.read_bytes()).hexdigest())
 (a.out/'SUBSTEP_MANIFEST.json').write_text(json.dumps(manifest,indent=2));sys.argv=[str(HERE/'live_tendon_trial.py'),'--protocol','hunt-ground','--physics','compliant','--source-task-basis','--config',str(a.config),'--out',str(a.out)]
 runpy.run_path(sys.argv[0],run_name='__main__')
if __name__=='__main__':main()
