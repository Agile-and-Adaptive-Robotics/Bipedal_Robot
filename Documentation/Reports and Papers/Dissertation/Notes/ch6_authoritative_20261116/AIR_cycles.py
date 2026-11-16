"""Run the archived principal-cycle calculation under the available base SciPy."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np
import scipy

root=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104')
module=root/'source_snapshot_final/controller_profiles_20261104/dominant_cycles.py'
spec=importlib.util.spec_from_file_location('archived_cycles',module)
cycles=importlib.util.module_from_spec(spec);spec.loader.exec_module(cycles)
out={'numpy':np.__version__,'scipy':scipy.__version__,'module':str(module),'module_sha256':hashlib.sha256(module.read_bytes()).hexdigest(),'self_check':cycles.self_check(),'stages':{}}
old=json.loads((root/'hunt_pose_air_validation_v3/COMPLETE.json').read_text())
for stage in (1,2):
    p=root/f'hunt_pose_air_validation_v3/stage{stage}/capture/aligned_trace.npz'
    with np.load(p) as z:
        t=z['t_state'];keep=(t>=5)&(t<14);names=z['joint_names'].tolist();q=z['qpos'][:,z['joint_qpos_addresses']]
        knee=-np.degrees(q[keep][:,[names.index('knee_angle_'+s)for s in ('r','l')]])
        hip=np.degrees(q[keep][:,[names.index('hip_flexion_'+s)for s in ('r','l')]])
        result=cycles.measure(t[keep],knee,hip)
    expected=old[stage-1]['principal_cycle_audit']
    differences={k:float(abs(result[k]-expected[k]))for k in ['period_s','bilateral_phase_error','period_cv']}
    out['stages'][str(stage)]={'source':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'window':'5<=t_state<14 s','recomputed':result,'principal_counts_equal':result['principal_knee_cycles']==expected['principal_knee_cycles'],'scalar_absolute_differences':differences}
target=Path(__file__).with_name('AIR_CYCLE_RESULTS.json');target.write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'self_check_passed':out['self_check']['passed'],'numpy':out['numpy'],'scipy':out['scipy'],'stages':{k:{'counts':v['recomputed']['principal_knee_cycles'],'period':v['recomputed']['period_s'],'phase_error':v['recomputed']['bilateral_phase_error'],'differences':v['scalar_absolute_differences']}for k,v in out['stages'].items()}},indent=2))
