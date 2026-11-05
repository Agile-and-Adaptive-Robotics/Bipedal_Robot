"""Held-out tests and reference-backend checks of a frozen postural bench fit."""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
from run import trial,electrical_check

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    best=json.loads((a.fit/'FIT_COMPLETE.json').read_text());parameters=best['parameters']
    a.output.mkdir(parents=True,exist_ok=False)
    check=electrical_check();(a.output/'ELECTRICAL_CHECK.json').write_text(json.dumps(check,indent=2))
    if not check['passed']:raise ValueError('Neural reference check failed')
    cases=[dict(initial=.04,push=8.,noise=0.,seed=61104),dict(initial=-.04,push=-8.,noise=0.,seed=61104)]
    cases += [dict(initial=q,push=f,noise=noise,seed=seed)for q,f in ((.08,-16.),(-.08,16.),(.08,16.),(-.08,-16.))for noise,seed in ((0.,61105),(.025,61106),(.05,61107))]
    rows=[]
    # Baseline and trained neural parameters are frozen across all conditions.
    for label,params in (('baseline',{}),('trained',parameters)):
        for i,case in enumerate(cases):
            row,log,neur=trial('premotor','inhibitory',80.,10.,parameters=params,**case)
            row.update(condition=i,parameter_set=label,held_out=i>=2);rows.append(row)
            if i in (0,2,4):np.savez_compressed(a.output/(label+'_'+str(i)+'.npz'),mechanics=log,neural=neur)
            (a.output/'RESULTS.json').write_text(json.dumps(rows,indent=2))
            print(json.dumps(dict(parameter_set=label,condition=i,accepted=row['accepted'],rms=row['final_rms_deg'])),flush=True)
    # Fresh-process training repeats and the native dense reference backend use
    # exactly the same physical and neural parameters.
    parity=[]
    for label,params in (('baseline',{}),('trained',parameters)):
        compact,c_log,c_neur=trial('premotor','inhibitory',80.,10.,parameters=params)
        native,n_log,n_neur=trial('premotor','inhibitory',80.,10.,parameters=params,backend='native')
        parity.append(dict(parameter_set=label,maximum_mechanical_difference=float(abs(c_log-n_log).max()),maximum_neural_difference_mV=float(abs(c_neur-n_neur).max()),compact_wall_s=compact['wall_s'],native_wall_s=native['wall_s'],native_over_compact=native['wall_s']/compact['wall_s'],accepted=compact['accepted']and native['accepted']))
    summary=dict(fit=str(a.fit),parameters=parameters,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in Path(__file__).parent.glob('*.py')},held_out={label:dict(passed=sum(r['accepted']for r in rows if r['parameter_set']==label and r['held_out']),total=sum(r['held_out']for r in rows if r['parameter_set']==label))for label in ('baseline','trained')},reference_backend=parity,scope='Fixed-base single-ankle proof of concept. No deployment into full Gait2392 controllers in this experiment.')
    (a.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
if __name__=='__main__':main()
