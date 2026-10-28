from pathlib import Path
import os,sys,json,subprocess,time,hashlib,concurrent.futures,importlib.util
R=Path('D:/GitHub/Bipedal_Robot');S=R/'Code/MuJoCo_SNS/spinal';O=S/'reports_20261028_publication_flow';P=Path(sys.executable)
for n in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[n]='1'
os.environ['CONDA_PREFIX']=str(P.parent);os.environ['PYTHONIOENCODING']='utf-8'
def child(name):
 out=O/name;out.mkdir(exist_ok=True);sys.path.insert(0,str(S));sys.path.insert(0,str(S/'reports_20261027_codex'))
 if name=='w2l':
  import recover_w2l as m
  m.HERE=out
  jobs=[dict(name='verified',knobs=json.loads((S/'reports_20261027_codex/w2l_verified_config.json').read_text())['knobs'],backend='exact_fast')]
  raw=json.loads((S/'reports_20261027_codex/perturbation_results.json').read_text())
  jobs += [dict(v['job'],backend='exact_fast') for v in raw if v['job']['name'] in ['panel_control','panel_del_all_feedback']]
  for j in jobs:
   inp=out/(j['name']+'_input.json');inp.write_text(json.dumps(j));m.one(inp)
   v=json.loads((out/(j['name']+'.json')).read_text());print(j['name'],json.dumps(v['metrics']),flush=True)
 elif name in ['deng_neural','rybak']:
  script=S/'replication'/('verify_deng.py' if name=='deng_neural' else 'verify_rybak2006.py');subprocess.run([str(P),str(script)],check=True,cwd=out)
 elif name=='li':subprocess.run([str(P),str(S/'w2l_mujoco/test_li_stepping_zdrop.py'),'--drop=0.05','--sink=0.01','--support=0.5'],check=True,cwd=out)
 elif name in ['dirusso_equations','dirusso_model']:
  f=S/'reports_20261027_codex'/('test_dirusso_equations.py' if name=='dirusso_equations' else 'inspect_author_model.py');code=f.read_text();code=code.replace("(Path(__file__).resolve().parent/'dirusso_neural_benchmark.json')","(OUT/'dirusso_neural_benchmark.json')").replace("(HERE/'author_model_inventory.json')","(OUT/'author_model_inventory.json')");exec(compile(code,str(f),'exec'),{'__file__':str(f),'__name__':'__main__','OUT':out})
if __name__=='__main__':
 if len(sys.argv)>1:child(sys.argv[1]);sys.exit()
 if O.exists():raise SystemExit('Output folder exists; refusing overwrite')
 O.mkdir();(O/'PREREG.md').write_text('Fresh-process reruns of existing W2L control/complete-feedback deletion, Deng translated-circuit persistence, Rybak deletion tests, Li MuJoCo contact test, and Di Russo equation parity/model load. No fitting or model edits. Existing gates retained. Native Di Russo gait is a separate unresolved runtime requirement. Stimulus matrix and IK overlays are additional recorded tests, not implicit replication claims.\n')
 def job(n):
  t=time.time()
  with (O/(n+'.log')).open('w',encoding='utf-8') as log:
   try:rc=subprocess.run([str(P),__file__,n],stdout=log,stderr=subprocess.STDOUT,timeout=240).returncode
   except subprocess.TimeoutExpired:rc=124
  v={'name':n,'returncode':rc,'seconds':time.time()-t};print(json.dumps(v),flush=True);return v
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:results=list(ex.map(job,['w2l','deng_neural','rybak','li','dirusso_equations','dirusso_model']))
 (O/'rerun_summary.json').write_text(json.dumps(results,indent=2));print('RESULTS '+str(O),flush=True)
