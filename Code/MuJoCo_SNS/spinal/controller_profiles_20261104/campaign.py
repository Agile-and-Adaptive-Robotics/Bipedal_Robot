"""Serial bounded campaign, with persistent records and phase gates.

Launch with myo. Each Rustuna block uses Python 3.11 plus one persistent myo
worker; no more than one scientific fit runs in this lane. No body physics is
changed. Existing completed blocks are reused, interrupted blocks preserved.
"""
import argparse,json,os,subprocess,sys,time,traceback
from pathlib import Path
from datetime import datetime
HERE=Path(__file__).resolve().parent
MYO=Path(r'D:\Anaconda\envs\myo\python.exe');OPT=Path(r'D:\Anaconda\python.exe')

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--hours',type=float,default=2.);p.add_argument('--families',nargs='+',default=['hunt','s3k','w2lvar','syn6']);p.add_argument('--rg-seed-file',type=Path);a=p.parse_args()
    if not 0<a.hours<=8.:raise ValueError('Use a bounded campaign of at most eight hours')
    if any(f not in ('hunt','s3k','w2lvar','syn6')for f in a.families):raise ValueError('Unknown family')
    a.output.mkdir(parents=True,exist_ok=True);deadline=time.time()+a.hours*3600
    env=os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):env[key]='1'
    env['CONDA_PREFIX']=str(MYO.parent)
    rows=[]
    def record():
        (a.output/'CAMPAIGN_STATUS.json').write_text(json.dumps(dict(machine=os.environ.get('COMPUTERNAME'),updated=datetime.now().isoformat(),rows=rows,scope='Graded neural replay fitting, not physical curriculum'),indent=2))
    def run(args,log):
        remaining=deadline-time.time()
        if remaining<30:raise TimeoutError('Campaign deadline reached')
        with log.open('a',encoding='utf-8')as stream:
            subprocess.run(list(map(str,args)),stdout=stream,stderr=subprocess.STDOUT,env=env,check=True,timeout=remaining,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    for family in a.families:
        folder=a.output/family
        row=dict(family=family,status='preparing');rows.append(row);record()
        try:
            if not (folder/'STATIC_CHECK.json').exists():
                args=[MYO,HERE/'prepare.py','--family',family,'--output',folder]
                if family=='hunt':args+=['--hunt-config',r'D:\Users\Ben Bolen\NeuroAcceleration_20261104\body\compact_final.json']
                elif a.rg_seed_file:args+=['--rg-seed-file',a.rg_seed_file]
                run(args,a.output/(family+'_prepare.log'))
            row['stages']={}
            for phase,ntrials in (('A',16),('B',32),('C',32),('D',24)):
                final=folder/'training'/phase/'FINAL.json'
                if not final.exists():
                    if final.parent.exists():raise RuntimeError('Interrupted block preserved; use a separate candidate folder')
                    row['status']='training_'+phase;record()
                    run([OPT,HERE/'train.py','--family',family,'--prepared',folder,'--phase',phase,'--trials',ntrials,'--budget-s',min(1500,int(deadline-time.time()-10)),'--execute','--computer-free','--runtime',r'D:\Users\Ben Bolen\ControllerPreparation_20261104\rustuna_runtime_py311'],a.output/(family+'_'+phase+'.log'))
                result=json.loads(final.read_text());row['stages'][phase]=dict(accepted=result['accepted'],value=result['best']['value'],diagnostics=result['best']['diagnostics'],file=str(final));record()
                if not result['accepted']:
                    row['status']='blocked_at_'+phase;break
            else:row['status']='offline_A_D_passed'
        except Exception as exc:
            row['status']='error';row['error']=repr(exc)
        record()
    (a.output/'CAMPAIGN_COMPLETE.json').write_text(json.dumps(dict(completed=datetime.now().isoformat(),rows=rows),indent=2))
if __name__=='__main__':main()
