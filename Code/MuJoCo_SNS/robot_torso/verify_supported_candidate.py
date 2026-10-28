"""Fresh-process full-trajectory equivalence check for the supported robot."""
import concurrent.futures,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def replay(backend,options):
    output=HERE/'runs'/('verified65_'+backend)
    command=[sys.executable,str(HERE/'run_neural_robot.py')]
    for key,value in dict(options,duration=40.,backend=backend,output=str(output)).items():
        if value is None:continue
        flag='--'+key.replace('_','-')
        if isinstance(value,bool):
            if value:command.append(flag)
        else:command.extend((flag,str(value)))
    with output.with_suffix('.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,cwd=HERE)
    return output

def main():
    options=json.loads((HERE/'runs/lateral03_cap65.json').read_text())['run_options']
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        paths=list(pool.map(lambda backend:replay(backend,options),('reference','exact_fast')))
    import numpy as np
    a,b=[np.load(path.with_suffix('.npz')) for path in paths]
    differences={key:float(np.max(np.abs(a[key]-b[key]))) if a[key].shape==b[key].shape else None for key in a.files}
    summaries=[json.loads(path.with_suffix('.json').read_text()) for path in paths]
    result=dict(machine='easteregg2',duration_s=40.,completed=[s['failed'] is None for s in summaries],
        bit_exact=all(value==0. for value in differences.values()),maximum_differences=differences,
        source_hashes_match=summaries[0]['source_hashes']==summaries[1]['source_hashes'],
        activation_floor=options['activation_floor'],vertical_support_cap=options['support_cap_final'],
        scope='Backend equivalence and supported survival only; stepping assessed separately')
    (HERE/'supported65_backend_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
