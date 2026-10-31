"""Bounded Rustuna orchestration through a persistent myo JSON-lines worker.

Model-specific source binding, measurements, baseline acceptance, and objective
belong in the adapter. This script never invents missing muscle correspondences.
"""
import argparse,sys,os,json,time,hashlib,math,subprocess,threading,queue,platform,importlib.metadata
from pathlib import Path
from rustuna_runtime import load_rustuna
HERE=Path(__file__).resolve().parent
MYO=r'D:\Anaconda\envs\myo\python.exe'

class Worker:
 def __init__(self,adapter,folder):
  self.log=(folder/'worker_stderr.log').open('a',encoding='utf-8');env=os.environ.copy();env['CONDA_PREFIX']=str(Path(MYO).parent)
  for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):env[k]='1'
  self.process=subprocess.Popen([MYO,'-u',str(HERE/'rustuna_myo_worker.py'),'--adapter',str(adapter)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,text=True,encoding='utf-8',env=env,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
  self.lines=queue.Queue()
  def read():
   for line in self.process.stdout:self.lines.put(line)
   self.lines.put(None)
  threading.Thread(target=read,daemon=True).start()
 def request(self,payload,timeout):
  self.process.stdin.write(json.dumps(payload,allow_nan=False)+'\n');self.process.stdin.flush()
  try:line=self.lines.get(timeout=timeout)
  except queue.Empty:raise TimeoutError('Scientific worker exceeded per-request timeout')
  if line is None:raise RuntimeError('Scientific worker exited; inspect worker_stderr.log')
  out=json.loads(line)
  if not out.get('ok'):raise RuntimeError('Scientific worker rejected request: '+json.dumps(out))
  return out
 def close(self):
  if self.process.poll() is None:
   try:self.request({'op':'shutdown'},5.);self.process.wait(timeout=5.)
   except Exception:self.process.terminate();self.process.wait(timeout=5.)
  self.log.close()

def sample_parameters(trial,fixed,space):
 params=dict(fixed)
 for key,bounds in space.items():params[key]=trial.suggest_float(key,float(bounds['low']),float(bounds['high']),log=bool(bounds.get('log',False)))
 return params

def replay_completed(study,rows,fixed,space):
 """Reconstruct sampler RNG state by exact sequential ask/tell replay."""
 for row in rows:
  trial=study.ask();params=sample_parameters(trial,fixed,space)
  if trial.number!=row['number'] or params!=row['params']:raise RuntimeError('Seeded replay diverged; refuse a mixed or incompatible study')
  study.tell(trial.number,float(row['value']))

def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--probe-only',action='store_true');p.add_argument('--resume',action='store_true');a=p.parse_args()
 manifest_path=Path(a.manifest).resolve();raw=manifest_path.read_bytes();cfg=json.loads(raw);adapter=Path(cfg['adapter_path']).resolve();folder=Path(cfg['output_dir']).resolve()
 if not folder.is_relative_to(HERE/'runs'):raise ValueError('Raw study outputs must remain in this package\'s ignored runs folder')
 if not adapter.is_file():raise FileNotFoundError(adapter)
 folder.mkdir(parents=True,exist_ok=True);trialpath=folder/'TRIALS.jsonl'
 if not a.probe_only and trialpath.exists() and not a.resume:raise RuntimeError('Output already has trials; use --resume for exact seeded replay')
 maxtrials=int(cfg.get('max_trials',64));budget=float(cfg.get('wall_budget_s',900.));timeout=float(cfg.get('trial_timeout_s',180.))
 if not 1<=maxtrials<=256 or not 0<budget<=10800 or not 0<timeout<=900:raise ValueError('Bounded study limits exceeded')
 fixed=cfg.get('fixed_parameters',{});space=cfg.get('search_space',{})
 if not space or set(fixed)&set(space):raise ValueError('Explicit non-overlapping fixed_parameters/search_space required')
 for key,bounds in space.items():
  lo,hi=float(bounds['low']),float(bounds['high'])
  if not math.isfinite(lo+hi) or not lo<hi or bounds.get('log',False) and lo<=0:raise ValueError('Invalid search bounds: '+key)
 initial=cfg.get('initial_candidates',[])
 if len(initial)>min(16,maxtrials):raise ValueError('Too many initial candidates for bounded study')
 for candidate in initial:
  if set(candidate)!=set(space):raise ValueError('Initial candidates must name exactly every sampled parameter')
  for key,value in candidate.items():
   if not float(space[key]['low'])<=float(value)<=float(space[key]['high']):raise ValueError('Initial candidate outside bounds: '+key)
 started=time.perf_counter();backend=load_rustuna();rows=[];best=None
 provenance=dict(machine=platform.node(),manifest_sha256=hashlib.sha256(raw).hexdigest(),adapter_sha256=hashlib.sha256(adapter.read_bytes()).hexdigest(),optimizer_python=sys.executable,worker_python=MYO,seed=int(cfg.get('seed',100031)),max_trials=maxtrials,wall_budget_s=budget,trial_timeout_s=timeout,rustuna_version=importlib.metadata.version('rustuna'))
 priorpath=folder/'PROVENANCE.json'
 if a.resume:
  old=json.loads(priorpath.read_text())
  for key,value in provenance.items():
   if old.get(key)!=value:raise RuntimeError('Resume provenance mismatch: '+key)
  provenance=old;rows=[json.loads(line) for line in trialpath.read_text().splitlines() if line.strip()] if trialpath.exists() else []
  if rows:best=min(rows,key=lambda r:r['value'])
 else:provenance['started_epoch_s']=time.time()
 deadline=provenance['started_epoch_s']+budget
 if not a.probe_only and time.time()>=deadline:raise TimeoutError('Original study wall-time budget exhausted; no new evaluations')
 worker=Worker(adapter,folder)
 (folder/'MANIFEST.json').write_bytes(raw);priorpath.write_text(json.dumps(provenance,indent=2))
 try:
  preflight=worker.request({'op':'preflight','manifest':cfg},timeout);(folder/'PREFLIGHT.json').write_text(json.dumps(preflight,indent=2))
  if a.probe_only:print(json.dumps(preflight,indent=2));return
  if preflight.get('ready_for_optimization') is not True:raise RuntimeError('Baseline/mapping gate has not passed: '+json.dumps(preflight))
  study=backend.create_study(direction='minimize',sampler=backend.samplers.TPESampler(seed=provenance['seed'],n_startup_trials=min(10,maxtrials),multivariate=False))
  for candidate in initial:study.enqueue_trial(candidate)
  replay_completed(study,rows,fixed,space)
  while len(rows)<maxtrials and time.time()<deadline:
   remaining=deadline-time.time();trial=study.ask();params=sample_parameters(trial,fixed,space)
   evaluation_id='trial_%04d_%d'%(trial.number,time.time_ns())
   request=dict(op='evaluate',trial_number=trial.number,evaluation_id=evaluation_id,evaluation_seed=provenance['seed']+trial.number,params=params,manifest=cfg);(folder/'LAST_REQUEST.json').write_text(json.dumps(request,indent=2))
   tick=time.perf_counter();out=worker.request(request,min(timeout,remaining));value=float(out['value'])
   if not math.isfinite(value):raise ValueError('Worker objective must be finite and carry explicit failure diagnostics')
   row=dict(number=trial.number,evaluation_id=evaluation_id,evaluation_seed=request['evaluation_seed'],params=params,value=value,evaluation_wall_s=time.perf_counter()-tick,diagnostics=out.get('diagnostics',{}));rows.append(row)
   with trialpath.open('a',encoding='utf-8') as f:f.write(json.dumps(row,allow_nan=False)+'\n');f.flush()
   study.tell(trial.number,value)
   if best is None or value<best['value']:best=row;(folder/'BEST.json').write_text(json.dumps(dict(provenance=provenance,trial=best),indent=2))
  result=dict(provenance=provenance,completed_trials=len(rows),this_process_wall_s=time.perf_counter()-started,original_budget_elapsed_s=time.time()-provenance['started_epoch_s'],best=best,scope='Candidate tuning only. Fresh-process repeat, held-out checks, and physical acceptance remain required.')
  (folder/'SUMMARY.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
 finally:worker.close()

if __name__=='__main__':main()
