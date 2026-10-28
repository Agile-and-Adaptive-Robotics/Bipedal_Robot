"""Fresh-process replay and matched-control test for the defense diagnostic."""
import json,subprocess,concurrent.futures,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
options=json.loads((HERE/'runs/abd_v2_control.json').read_text())['run_options']
def execute(label,updates):
    output=HERE/'runs'/label
    settings=dict(options,output=str(output));settings.update(updates)
    command=[sys.executable,str(HERE/'run_neural_robot.py')]
    for key,value in settings.items():
        if value is None:continue
        if isinstance(value,bool):
            if value:command.append('--'+key.replace('_','-'))
        else:command.extend(('--'+key.replace('_','-'),str(value)))
    with output.with_suffix('.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    return output
if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        tasks=[pool.submit(execute,'abd_defense_reference',dict(backend='reference')),
               pool.submit(execute,'abd_defense_exact',dict(backend='exact_fast')),
               pool.submit(execute,'abd_defense_control_off',dict(backend='exact_fast',abduction_gain=0.,abductor_floor=.05))]
        paths=[task.result() for task in tasks]
    import numpy as np
    a,b=[np.load(path.with_suffix('.npz')) for path in paths[:2]]
    differences={key:float(np.max(abs(a[key]-b[key]))) if a[key].shape==b[key].shape else None for key in a.files}
    summaries=[json.loads(path.with_suffix('.json').read_text()) for path in paths]
    result=dict(bit_exact=all(value==0. for value in differences.values()),differences=differences,
        complete=[s['failed'] is None for s in summaries],source_hashes_match=summaries[0]['source_hashes']==summaries[1]['source_hashes'],
        scope='Numerical replay, not a gait or robustness acceptance test')
    (HERE/'abduction_defense_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
