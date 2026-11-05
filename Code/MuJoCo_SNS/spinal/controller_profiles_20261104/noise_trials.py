"""Fixed-weight air-coordination checks with measured, seeded RG current noise.

These are stochastic current disturbances, not neural parameter perturbations.
The spring-supported air protocol cannot establish standing robustness.
"""
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
from dominant_cycles import measure
HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    rows=[];cases=[('zero',0.,41001)]+[(str(sigma)+'_'+str(seed),sigma,seed)for sigma in(.01,.025)for seed in(41001,41002,41003)]+[('repeat',.025,41001)]
    for label,sigma,seed in cases:
        folder=a.output/label
        with(a.output/(label+'.log')).open('a',encoding='utf-8')as log:
            run=subprocess.run([sys.executable,str(HERE/'body_probe.py'),'--reflex-fit',str(a.fit),'--stage','2','--expression','.1','--noise-rg-sigma',str(sigma),'--noise-seed',str(seed),'--output',str(folder)],stdout=log,stderr=subprocess.STDOUT,timeout=600,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        row=dict(label=label,sigma_nA=sigma,seed=seed,source=str(folder),exit_code=run.returncode,accepted=False)
        path=folder/'body.result.json'
        if path.exists():
            old=json.loads(path.read_text());z=np.load(folder/'capture/aligned_trace.npz');t=z['t_state'];keep=(t>=5.)&(t<14.);names=z['joint_names'].tolist();q=z['qpos'][:,z['joint_qpos_addresses']]
            principal=measure(t[keep],-np.degrees(q[keep][:,[names.index('knee_angle_'+s)for s in('r','l')]]),np.degrees(q[keep][:,[names.index('hip_flexion_'+s)for s in('r','l')]]))
            gate=old['assessment'];row.update(accepted=bool(principal['accepted']and all(gate[k]for k in('finite','complete','frontal_gate','knee_extension_gate'))),principal=principal,maximum_knee_extension_deg=gate['maximum_knee_extension_deg'],measured=old['measured'],wall_s=old['wall_s'])
            if sigma:
                recorded=z['stimulus_nA'];finite=recorded[:,np.any(np.nan_to_num(recorded)!=0,axis=0)];row['observed_noise']=dict(channels=finite.shape[1],mean_nA=np.mean(finite,axis=0).tolist(),std_nA=np.std(finite,axis=0).tolist())
            if label=='repeat':
                first=np.load(a.output/'0.025_41001/capture/aligned_trace.npz');row['fresh_repeat_max_abs_difference']={key:float(np.max(abs(z[key]-first[key])))for key in('qpos','neural_voltage','stimulus_nA')}
        elif(folder/'STOP.json').exists():row['stop']=json.loads((folder/'STOP.json').read_text())
        rows.append(row);(a.output/'STATUS.json').write_text(json.dumps(rows,indent=2))
    summary=dict(cases=rows,passed=sum(row['accepted']for row in rows),total=len(rows),scope='Fixed learned parameters; independent 20 ms OU current at four RG cells; suspended-air motion gates, not unsupported balance or ground stepping')
    (a.output/'COMPLETE.json').write_text(json.dumps(summary,indent=2))
if __name__=='__main__':main()
