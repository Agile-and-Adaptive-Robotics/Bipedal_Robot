"""Seeded delay/noise stress tests with frozen postural parameters and physics."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from run import trial

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    parameters=json.loads((a.fit/'FIT_COMPLETE.json').read_text())['parameters']
    old=np.load(a.fit.parent/'initial_k80_d10/premotor_inhibitory.npz');r,log,neural=trial('premotor','inhibitory',80.,10.)
    difference=max(float(abs(log-old['mechanics']).max()),float(abs(neural-old['neural']).max()))
    (a.output/'ZERO_DELAY_PARITY.json').write_text(json.dumps(dict(maximum_absolute_difference=difference,passed=difference==0.),indent=2))
    if difference!=0.:raise RuntimeError('Zero-delay extension changed the baseline')
    rows=[]
    for label,values in (('baseline',{}),('trained',parameters)):
        for delay in (0.,.03,.07):
            for noise in (.025,.05):
                for seed in (61120,61121,61122):
                    for sign in (1.,-1.):
                        row,_,_=trial('premotor','inhibitory',80.,10.,initial=sign*.08,push=-sign*16.,noise=noise,seed=seed,parameters=values,backend='native',delay_s=delay)
                        row['parameter_set']=label;rows.append(row);(a.output/'RESULTS.json').write_text(json.dumps(rows,indent=2))
                print(json.dumps(dict(parameter_set=label,delay_s=delay,noise_nA=noise,completed=len(rows))),flush=True)
    summary=dict(source_fit=str(a.fit),parameters=parameters,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in (Path(__file__),Path(__file__).parent/'run.py',Path(__file__).parent/'model.py')},groups=[dict(parameter_set=label,delay_s=delay,passed=sum(r['accepted']for r in rows if r['parameter_set']==label and r['added_sensor_delay_s']==delay),total=sum(r['parameter_set']==label and r['added_sensor_delay_s']==delay for r in rows))for label in ('baseline','trained')for delay in (0.,.03,.07)],scope='Uniform added measurement delay stress test on the isolated ankle bench, not calibrated physiological conduction delays')
    (a.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
if __name__=='__main__':main()
