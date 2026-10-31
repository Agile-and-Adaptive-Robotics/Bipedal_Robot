"""Known-objective validation of exact JSON checkpoint/sampler-state replay."""
import json,copy
from rustuna_runtime import load_rustuna,HERE
from rustuna_transfer import sample_parameters,replay_completed
R=load_rustuna();fixed={'held':1.};space={'x':{'low':-2.,'high':3.},'y':{'low':.01,'high':10.,'log':True}}
def create():
 study=R.create_study(direction='minimize',sampler=R.samplers.TPESampler(seed=100031,n_startup_trials=10,multivariate=False))
 study.enqueue_trial({'x':0.,'y':1.});return study
def append(study,rows,n):
 for _ in range(n):
  t=study.ask();p=sample_parameters(t,fixed,space);value=(p['x']-.2)**2+(p['y']-2.)**2;study.tell(t.number,value);rows.append(dict(number=t.number,params=p,value=value))
full=[];append(create(),full,64);assert full[0]['params']=={'held':1.,'x':0.,'y':1.}
saved=json.loads(json.dumps(full[:32]));resumed=create();replay_completed(resumed,saved,fixed,space);append(resumed,saved,32);assert saved==full
broken=copy.deepcopy(full[:2]);broken[0]['params']['x']+=.1
try:replay_completed(create(),broken,fixed,space)
except RuntimeError:rejected=True
else:raise AssertionError('Changed checkpoint was not rejected')
out=HERE/'runs/rustuna_capability';out.mkdir(parents=True,exist_ok=True)
result=dict(scope='Cheap known objective only; no simulations',trials=64,checkpoint_after=32,full_history_exact_after_JSON_checkpoint_replay=True,corrupted_checkpoint_rejected=rejected,enqueue_trial_exact_candidate_verified=True,existing_databases_opened=False)
(out/'SEEDED_RESUME_SELFCHECK.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
