"""Actual-body central-feedback calibration after a passed deafferented stage.

This scales mirrored RG sensory projections, not physical damping or muscles.
The multiplier is inert in stage 1 because all afferent inputs/gains are off.
A winner is rerun through both stages in fresh processes before acceptance.
"""
import argparse,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from dominant_cycles import measure,self_check
HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--passed-air',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--pathway',choices=('RG','localII'),default='RG');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    previous=json.loads((a.passed_air/'COMPLETE.json').read_text())
    if not previous[0]['accepted']or previous[0]['stage']!=1:raise ValueError('Verified physical stage 1 required')
    if not self_check()['passed']:raise ValueError('Metric self-check failed')
    rows=[];best=None
    for scale in (.05,.1,.2,.4,.7):
        folder=a.output/('scale_'+str(scale))
        with(a.output/('scale_'+str(scale)+'.log')).open('a',encoding='utf-8')as log:
            run=subprocess.run([sys.executable,str(HERE/'body_probe.py'),'--reflex-fit',str(a.fit),'--stage','2','--expression','.1','--feedback-scale'if a.pathway=='RG'else'--local-ii-scale',str(scale),'--output',str(folder)],stdout=log,stderr=subprocess.STDOUT,timeout=600,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        result_path=folder/'body.result.json'
        row=dict(scale=scale,pathway=a.pathway,accepted=False,folder=str(folder))
        if result_path.exists():
            result=json.loads(result_path.read_text());z=np.load(folder/'capture/aligned_trace.npz');t=z['t_state'];keep=(t>=5.)&(t<14.);names=z['joint_names'].tolist();q=z['qpos'][:,z['joint_qpos_addresses']]
            principal=measure(t[keep],-np.degrees(q[keep][:,[names.index('knee_angle_'+s)for s in('r','l')]]),np.degrees(q[keep][:,[names.index('hip_flexion_'+s)for s in('r','l')]]));old=result['assessment']
            row.update(accepted=bool(principal['accepted']and old['finite']and old['complete']and old['frontal_gate']and old['knee_extension_gate']),principal=principal,legacy=old,measured=result['measured'],wall_s=result['wall_s'])
            if row['accepted']:
                row['cost']=abs(principal['period_s']/1.233-1.)+principal['bilateral_phase_error']
                if best is None or row['cost']<best['cost']:best=row;(a.output/'BEST.json').write_text(json.dumps(best,indent=2))
        elif(folder/'STOP.json').exists():row['stop']=json.loads((folder/'STOP.json').read_text())
        else:row['error']='Exit without complete structured result: '+str(run.returncode)
        rows.append(row);(a.output/'STATUS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(scale=scale,accepted=row['accepted'],principal=row.get('principal'),stop=row.get('stop'))),flush=True)
    if best:
        with(a.output/'fresh_validation.log').open('a',encoding='utf-8')as log:subprocess.run([sys.executable,str(HERE/'air_validation.py'),'--fit',str(a.fit),'--output',str(a.output/'fresh_validation'),'--expression','.1','--feedback-scale'if a.pathway=='RG'else'--local-ii-scale',str(best['scale'])],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1200,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    (a.output/'COMPLETE.json').write_text(json.dumps(dict(rows=rows,best=best),indent=2))
if __name__=='__main__':main()
