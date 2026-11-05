"""Bounded air-expression tests, with no automatic advancement of failures."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    rows=[];deadline=time.time()+3600
    for expression in (.1,.2,.35,.5,.7):
        row=dict(expression=expression,stages={});rows.append(row)
        for stage in (1,2):
            folder=a.output/('expression_'+str(expression))/('stage'+str(stage));folder.parent.mkdir(parents=True,exist_ok=True)
            completed=folder/'body.result.json';stop=folder/'STOP.json'
            if not completed.exists()and not stop.exists():
                if folder.exists():raise RuntimeError('Interrupted probe preserved: '+str(folder))
                with(folder.parent/('stage'+str(stage)+'.log')).open('a',encoding='utf-8')as stream:
                    try:run=subprocess.run([sys.executable,str(HERE/'body_probe.py'),'--reflex-fit',str(a.fit),'--stage',str(stage),'--expression',str(expression),'--output',str(folder)],stdout=stream,stderr=subprocess.STDOUT,timeout=min(600,deadline-time.time()),creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
                    except subprocess.TimeoutExpired:row['stages'][str(stage)]={'accepted':False,'failure':'timeout'};break
            if completed.exists():
                result=json.loads(completed.read_text());row['stages'][str(stage)]=dict(accepted=result['assessment']['accepted'],assessment=result['assessment'],measured=result['measured'],wall_s=result['wall_s'],simulated_s=result['simulated_s'])
            elif stop.exists():row['stages'][str(stage)]=dict(accepted=False,stop=json.loads(stop.read_text()))
            else:row['stages'][str(stage)]=dict(accepted=False,failure='Probe exited without result; inspect log')
            (a.output/'STATUS.json').write_text(json.dumps(rows,indent=2))
            if not row['stages'][str(stage)]['accepted']:break
        if all(row['stages'].get(str(s),{}).get('accepted',False)for s in (1,2)):break
    (a.output/'COMPLETE.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':main()
