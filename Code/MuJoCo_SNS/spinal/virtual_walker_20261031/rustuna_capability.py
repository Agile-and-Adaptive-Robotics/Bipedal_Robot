"""Cheap seeded optimizer checks; never calls a walker or existing study database."""
import argparse,sys,os,time,json,platform,hashlib,importlib.metadata as M
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3]
p=argparse.ArgumentParser();p.add_argument('--backend',choices=['rustuna','optuna'],default='rustuna');p.add_argument('--trials',type=int,default=256);a=p.parse_args()
if a.backend=='rustuna':sys.path.insert(0,str(REPO/'tmp/rustuna_runtime_py311'))
started=time.perf_counter();backend=__import__(a.backend);import_s=time.perf_counter()-started
if hasattr(backend,'logging'):backend.logging.set_verbosity(backend.logging.ERROR)
seed=100031;sampler_config=dict(seed=seed,n_startup_trials=10,multivariate=False)
def exercise():
 rows=[];study=backend.create_study(direction='minimize',sampler=backend.samplers.TPESampler(**sampler_config));clock=time.perf_counter()
 def objective(trial):
  params={f'x{i}':trial.suggest_float(f'x{i}',-2.,2.) for i in range(4)}
  value=sum((params[f'x{i}']-v)**2 for i,v in enumerate((.2,-.4,.6,-.8)))
  rows.append(dict(number=trial.number,params=params,value=value));return value
 study.optimize(objective,n_trials=a.trials);elapsed=time.perf_counter()-clock
 best=min(rows,key=lambda x:x['value']);assert dict(study.best_trial.params)==best['params']
 return rows,elapsed,dict(params=best['params'],value=best['value'])
one,elapsed,best=exercise();two,repeat_s,repeat_best=exercise();assert one==two,'Same-backend fresh-study seed reproducibility failed'
out=HERE/'runs/rustuna_capability';out.mkdir(parents=True,exist_ok=True)
trialpath=out/(a.backend+'_seeded_trials.json');trialpath.write_text(json.dumps(one,indent=2))
summary=dict(machine=platform.node(),python=sys.executable,python_version=platform.python_version(),backend=a.backend,version=M.version(a.backend),backend_file=backend.__file__,
 objective='Four-dimensional known quadratic; optimum[.2,-.4,.6,-.8], value0; no simulations',sampler='TPE',sampler_config=sampler_config,trials=a.trials,import_wall_s=import_s,study_wall_s=elapsed,repeat_wall_s=repeat_s,wall_per_trial_ms=elapsed/a.trials*1000,
 best=best,exact_seed_repeat=True,full_trials=str(trialpath),full_trials_sha256=hashlib.sha256(trialpath.read_bytes()).hexdigest(),
 scope='Optimizer/suggest/objective bookkeeping only. Separate Python versions. Not an end-to-end MuJoCo speed-up. Existing databases untouched.')
(out/(a.backend+'_capability.json')).write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
