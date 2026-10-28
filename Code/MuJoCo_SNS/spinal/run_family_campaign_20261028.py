"""Bounded isolated-process screen; existing optimizer databases untouched."""
import sys,subprocess,json,concurrent.futures,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;OUT=HERE/'reports_20261028_family_audit/human_trials';OUT.mkdir(parents=True,exist_ok=True)
CASES={'s3k':['baseline','rg_pf_gate','ib_half','afferent_half'],
       'w2lvar':['baseline','central_half','comm_half','central_quarter'],
       'syn6':['baseline','central05','central15','pf_gate','central_gate']}
def execute(family,case,reference=False):
    label=family+'_'+case+('_reference' if reference else '')
    if (OUT/(label+'.json')).exists():return {'label':label,'existing':True}
    command=[sys.executable,str(HERE/'family_human_trials_20261028.py'),family,case]+(['--reference'] if reference else [])
    with (OUT/(label+'.log')).open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
    if r.returncode:return {'label':label,'error_exit':r.returncode}
    result=json.loads((OUT/(label+'.json')).read_text());m=result['metrics']
    return dict(label=label,kine_score=m['kine_score'],nan=m['nan'],fell=m['bal_fell'],t_end=m['t_end'],wall_s=result['wall_s'])
if __name__=='__main__':
    start=time.perf_counter();results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        jobs=[pool.submit(execute,f,c) for f,cases in CASES.items() for c in cases]
        jobs += [pool.submit(execute,f,'baseline',True) for f in CASES]
        for job in concurrent.futures.as_completed(jobs):
            result=job.result();results.append(result);print(json.dumps(result),flush=True)
            (OUT/'campaign_progress.json').write_text(json.dumps(dict(completed=results,expected=len(jobs),wall_s=time.perf_counter()-start),indent=2))
