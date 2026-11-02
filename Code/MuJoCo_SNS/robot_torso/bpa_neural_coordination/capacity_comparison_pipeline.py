"""One heavy process, restartable bounded continuation of registered studies."""
import os,sys,time,json,subprocess
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=r'D:\Anaconda\envs\myo'
HERE=Path(__file__).resolve().parent;DEST=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\pressure_reserve_capacity');DRIVER=HERE.parents[1]/'spinal/virtual_walker_20261031/rustuna_transfer.py'

def command(args,name):
 with (DEST/(name+'.log')).open('a',encoding='utf-8') as log:proc=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
 if proc.returncode:raise RuntimeError('Pipeline stage failed: '+name)

def report():command([sys.executable,str(HERE/'capacity_comparison_report.py')],'report')

def main():
 # The already running pressure search exclusively owns the numerical worker.
 deadline=time.monotonic()+1500.
 while not (DEST/'pressure/rustuna/SUMMARY.json').exists():
  if time.monotonic()>deadline:raise TimeoutError('Pressure study did not finish; no second physics worker started')
  time.sleep(4.)
 report()
 if not (DEST/'ZERO_PHYSICAL_PARITY.json').exists():command([sys.executable,str(HERE/'capacity_estimator_checks.py'),'--physical-zero'],'zero_physical')
 force=DEST/'force/rustuna'
 if not (force/'SUMMARY.json').exists():
  args=[r'D:\Anaconda\python.exe',str(DRIVER),'--manifest',str(DEST/'force_MANIFEST.json')]
  if (force/'TRIALS.jsonl').exists():args.append('--resume')
  command(args,'force_search')
 report()
 for mode in ('pressure','force'):
  if not (DEST/mode/'VALIDATION_STATUS.json').exists():command([sys.executable,str(HERE/'capacity_comparison_validation.py'),'--mode',mode],mode+'_validation')
  report()
 (DEST/'PIPELINE_COMPLETE.json').write_text(json.dumps(dict(completed=True,scope='All preregistered searches complete; each mode independently gated before validation; rejected modes do not receive newvalidation.',machine='easteregg2'),indent=2),encoding='utf-8')

if __name__=='__main__':main()
