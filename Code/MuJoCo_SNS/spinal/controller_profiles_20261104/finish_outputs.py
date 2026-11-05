"""Serial bounded C/D fits for families with accepted rhythm and timing fits."""
import argparse,json,subprocess,sys,time,os
from pathlib import Path
HERE=Path(__file__).resolve().parent;MYO=Path(r'D:\Anaconda\envs\myo\python.exe');OPT=Path(r'D:\Anaconda\python.exe')

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);deadline=time.time()+7200
    timing=dict(hunt=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104\recalibrated_v2\hunt\training\B'),s3k=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104\intrinsic_seed_v3\s3k\training\B'),w2lvar=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104\commissural_fit_v4\w2lvar\training\B'),syn6=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104\commissural_fit_v4\syn6\training\B'))
    rows=[]
    for family,source in timing.items():
        row=dict(family=family,source_timing=str(source),status='checking');rows.append(row)
        try:
            if not json.loads((source/'FINAL.json').read_text())['accepted']:raise ValueError('Timing dependency failed')
            for phase in ('C','D'):
                folder=a.output/family/phase;final=folder/'FINAL.json';folder.parent.mkdir(parents=True,exist_ok=True)
                row['status']='fitting_'+phase;(a.output/'STATUS.json').write_text(json.dumps(rows,indent=2))
                if not final.exists():
                    if folder.exists():raise RuntimeError('Interrupted output preserved: '+str(folder))
                    if phase=='C':
                        args=[MYO,HERE/'refit_motor_peaks.py','--timing-block',source,'--output',folder,'--iterations','24','--feedforward-only']
                        if family=='hunt':args+=['--pose-aware']
                    else:args=[OPT,HERE/'refine_reflex.py','--fit',folder.parent/'C','--output',folder,'--trials','24']
                    with(a.output/(family+'_'+phase+'.log')).open('a',encoding='utf-8')as log:subprocess.run(list(map(str,args)),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=min(1800,deadline-time.time()),creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
                result=json.loads(final.read_text());row[phase]=dict(accepted=result['accepted'],file=str(final),diagnostics=result['best']['diagnostics'])
                if not result['accepted']:row['status']='blocked_at_'+phase;break
            else:row['status']='offline_C_D_passed'
        except Exception as exc:row['status']='error';row['error']=repr(exc)
        (a.output/'STATUS.json').write_text(json.dumps(rows,indent=2))
    (a.output/'COMPLETE.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':main()
