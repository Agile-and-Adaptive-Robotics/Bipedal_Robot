import json
from pathlib import Path
import numpy as np
H=Path(__file__).resolve().parent
lines=['# Stage 3 task-balance implementation and measurements','',
'Machine: easteregg2. Eight added SNS cells: four directional task populations and four inhibitory relays. Controller total: 1,029 neurons, 424 inputs. No runtime MN target cancellation.',
'','CoM displacement/velocity use 100 ms delay and a 40 ms velocity/acceleration filter by default. Acceleration and somatic/orientation cue weighting are separate experimental options. Ground-truth simulation signals are sensory surrogates, not a model of receptor transduction or adaptive sensory reweighting.',
'','Recruitment is calibrated once at the initial pose from finite muscle-activation effects on CoM acceleration under contact. It is held fixed during each trial. Locomotor RG/PF fits are preserved.',
'','| Campaign/case | Stage | Pass | Sway (m) | Horizontal harness RMS / BW | Foot load sum / BW | Wall (s) |',
'|---|---:|---|---:|---:|---:|---:|']
bench=[]
for folder in sorted((H/'runs').glob('balance_20261029_*')):
    for path in sorted(folder.glob('*.result.json')):
        r=json.loads(path.read_text());a=r['assessment'];m=r['measured']
        lines.append(f"| {folder.name}/{path.name[:-12]} | {r['stage']} | {a['accepted']} | {a.get('sway_m',float('nan')):.4f} | {m['horizontal_harness_rms_bw']:.4f} | {sum(m['mean_foot_normal_bw']):.3f} | {r['wall_s']:.2f} |")
    sparse=folder/'benchmark_sparse.result.json';native=folder/'benchmark_native.result.json'
    if sparse.exists() and native.exists():
        s=json.loads(sparse.read_text());n=json.loads(native.read_text());errors={}
        with np.load(folder/'benchmark_sparse.npz',allow_pickle=True) as zs,np.load(folder/'benchmark_native.npz',allow_pickle=True) as zn:
            for name in ('qfull','ctrl','com','neuro'):
                errors[name]=float(np.max(abs(zs[name]-zn[name])))
        row=dict(campaign=folder.name,simulated_s=s['simulated_s'],sparse_wall_s=s['wall_s'],native_wall_s=n['wall_s'],
            end_to_end_ratio=n['wall_s']/s['wall_s'],sparse_neural_step_s=s['compiled']['neural_step_wall_s'],
            native_neural_step_s=n['compiled']['neural_step_wall_s'],maximum_array_differences=errors)
        bench.append(row)
lines+=['','## Matched backend comparisons','',
'Same fitted model, network topology, initialization, assistance, duration, and logging. Native uses the dense SNS fixed-tau equations; accelerated uses the sparse preallocated Numba equations. Background trial load varies, so elapsed times are observations, not controlled hardware benchmarks. Whole-run times include model/network setup, standing solve, calibration, physics, sensing, and file output. Neural-step times also include live feature and conductance updates.']
for row in bench:lines+=['','```json',json.dumps(row,indent=2),'```']
lines+=['','v1 tests replaced all inherited BAL inputs. v2 used a neutral standing pose; the inherited controller passes its initial screen, while replacement task feedback does not. v3 retains the inherited controller and tests additive task feedback against an identical-topology disabled-pathway control. Distributed harness restrictions remain assistance. No independent balance or curriculum graduation follows merely from a stage-3 pass.',
'','The v3 winner file selects the best **nonzero** task-feedback candidate for comparison, not an established improvement over the disabled-task control. The latter has slightly lower harness reaction and CoM displacement in these initial tests. Do not promote the added pathway as beneficial on this evidence.',
'','Prescribed-support cases explicitly remove the vertical rig spring/damper and apply a constant 60% or 40% body-weight force. They differ from stiffness-reduction tests. Both stages 3 and 4 remain failed under those tested support settings; inspect symmetry, positional error, phase, and contact criteria rather than using force support alone.',
'','Air stages named curriculum_air_stage_* retain their calibrated air initialization. Earlier curriculum_stage_1/2 neutral-standing-start probes are retained separately and do not replace the air-stage protocol.',
'','Contact sensing: three forces per foot are already logged. Current neural contact/timing pathways mainly receive scalar contact/load, not the full 3D force vector. COP and ankle-referenced contact moments are useful next measurements, but are not currently neural inputs.',
'','Suggested lineage label: biw2lvar-Hunt-PT experimental. This denotes W2L-variant ancestry plus biarticular integration, 3D joint PF channels, Hunt fitting, and the task-balance pathway. It is not an unchanged W2L model, a literature-established name, or a claim of validated walking.']
(H/'BALANCE_RESULTS.md').write_text('\n'.join(lines)+'\n');(H/'runs/backend_benchmarks.json').write_text(json.dumps(bench,indent=2))
print(json.dumps(bench,indent=2))
