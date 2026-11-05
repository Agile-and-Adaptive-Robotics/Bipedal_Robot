"""Repeat the safe air candidate, audit principal cycles, then gate stage 2."""
import argparse,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from dominant_cycles import self_check,measure
HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--expression',type=float,default=.1);p.add_argument('--feedback-scale',type=float,default=1.);p.add_argument('--local-ii-scale',type=float,default=1.);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    checks=self_check();(a.output/'METRIC_SELF_CHECK.json').write_text(json.dumps(checks,indent=2))
    if not checks['passed']:raise ValueError('Physical metric check failed')
    rows=[]
    for stage in (1,2):
        folder=a.output/('stage'+str(stage))
        with(a.output/('stage'+str(stage)+'.log')).open('a',encoding='utf-8')as log:
            result=subprocess.run([sys.executable,str(HERE/'body_probe.py'),'--reflex-fit',str(a.fit),'--stage',str(stage),'--expression',str(a.expression),'--feedback-scale',str(a.feedback_scale),'--local-ii-scale',str(a.local_ii_scale),'--output',str(folder)],stdout=log,stderr=subprocess.STDOUT,timeout=600,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        path=folder/'body.result.json'
        if not path.exists():
            rows.append(dict(stage=stage,accepted=False,stop=json.loads((folder/'STOP.json').read_text())if(folder/'STOP.json').exists()else dict(exit_code=result.returncode)));break
        old=json.loads(path.read_text());z=np.load(folder/'capture/aligned_trace.npz');t=z['t_state'];keep=(t>=5.)&(t<14.);names=z['joint_names'].tolist();q=z['qpos'][:,z['joint_qpos_addresses']]
        principal=measure(t[keep],-np.degrees(q[keep][:,[names.index('knee_angle_'+s)for s in('r','l')]]),np.degrees(q[keep][:,[names.index('hip_flexion_'+s)for s in('r','l')]]))
        gate=old['assessment'];accepted=bool(principal['accepted']and gate['finite']and gate['complete']and gate['frontal_gate']and gate['knee_extension_gate'])
        row=dict(stage=stage,accepted=accepted,legacy_assessment=gate,principal_cycle_audit=principal,measured=old['measured'],wall_s=old['wall_s'],simulated_s=old['simulated_s'],scope='Suspended-air validation only; legacy assessment retained; no ground stepping or body-balance claim')
        rows.append(row);(a.output/'STATUS.json').write_text(json.dumps(rows,indent=2))
        if not accepted:break
    (a.output/'COMPLETE.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':main()
