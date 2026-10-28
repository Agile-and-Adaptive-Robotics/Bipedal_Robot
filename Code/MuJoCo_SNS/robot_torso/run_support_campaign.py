"""Small, resumable fresh-process robot support campaign on easteregg2."""
import argparse,concurrent.futures,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
CASES={
 'slow0525_cap65':dict(load_compensation=1.,support_cap_final=.65,knee_amplitude=.6,knee_amplitude_right=.85,hip_flexion=.35,hip_extension=-.12,pattern='event',activation_floor=.05,unload_swing=True,lateral_shift=.03,lateral_stiffness=1000,lateral_damping=120,horizontal_stiffness=200,rg_tau_h=.525),
 'slow070_cap65':dict(load_compensation=1.,support_cap_final=.65,knee_amplitude=.6,knee_amplitude_right=.85,hip_flexion=.35,hip_extension=-.12,pattern='event',activation_floor=.05,unload_swing=True,lateral_shift=.03,lateral_stiffness=1000,lateral_damping=120,horizontal_stiffness=200,rg_tau_h=.7),
 'slow105_cap65':dict(load_compensation=1.,support_cap_final=.65,knee_amplitude=.6,knee_amplitude_right=.85,hip_flexion=.35,hip_extension=-.12,pattern='event',activation_floor=.05,unload_swing=True,lateral_shift=.03,lateral_stiffness=1000,lateral_damping=120,horizontal_stiffness=200,rg_tau_h=1.05),
}

def execute(item,duration):
    name,knobs=item;output=HERE/'runs'/name
    if output.with_suffix('.json').exists():return json.loads(output.with_suffix('.json').read_text())
    command=[sys.executable,str(HERE/'run_neural_robot.py'),'--mode','step','--duration',str(duration),'--stand-time','5','--output',str(output)]
    for key,value in knobs.items():
        if isinstance(value,bool):
            if value:command.append('--'+key.replace('_','-'))
        else:command.extend(('--'+key.replace('_','-'),str(value)))
    with output.with_suffix('.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,cwd=HERE)
    result=json.loads(output.with_suffix('.json').read_text())
    return dict(run=name,failed=result['failed'],duration=result['duration_simulated_s'],support=result['mean_support_fraction'],displacement=result['pelvis_forward_displacement_m'])

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--duration',type=float,default=20)
    args=parser.parse_args()
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(execute,item,args.duration) for item in CASES.items()]
        for future in concurrent.futures.as_completed(futures):print(json.dumps(future.result()),flush=True)
