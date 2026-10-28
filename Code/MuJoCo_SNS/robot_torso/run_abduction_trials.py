"""Fresh-process, resumable abductor activation diagnostics on easteregg2."""
import concurrent.futures,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
base=json.loads((HERE/'runs/verified65_exact_fast.json').read_text())['run_options']
cases={
    'abd_stop65':dict(model=str(HERE/'ah001_abduction_diagnostic.xml'),abduction_gain=160.,stance_width=.22,abduction_force_limit=25.,hip_frontal_weight=4.,foot_heading_gain=8.),
    'abd_stop80':dict(model=str(HERE/'ah001_abduction_diagnostic.xml'),abduction_gain=160.,stance_width=.22,abduction_force_limit=25.,hip_frontal_weight=4.,foot_heading_gain=8.,support_cap_final=.8),
    'abd_stop90':dict(model=str(HERE/'ah001_abduction_diagnostic.xml'),abduction_gain=160.,stance_width=.22,abduction_force_limit=25.,hip_frontal_weight=4.,foot_heading_gain=8.,support_cap_final=.9),
}
def execute(item):
    name,updates=item;output=HERE/'runs'/name
    options=dict(base,duration=20.,output=str(output),**updates)
    command=[sys.executable,str(HERE/'run_neural_robot.py')]
    for key,value in options.items():
        if value is None:continue
        if isinstance(value,bool):
            if value:command.append('--'+key.replace('_','-'))
        else:command.extend(('--'+key.replace('_','-'),str(value)))
    with output.with_suffix('.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    result=json.loads(output.with_suffix('.json').read_text())
    return dict(run=name,failed=result['failed'],duration=result['duration_simulated_s'])
if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for result in pool.map(execute,cases.items()):print(json.dumps(result),flush=True)
