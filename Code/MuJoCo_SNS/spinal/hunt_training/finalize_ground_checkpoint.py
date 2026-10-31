"""Freeze current ground evidence without rewriting historical trial results."""
import os,json,hashlib,shutil
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';folder=root/'aligned_support_ladder_v1';p=folder/'repeat_0p5'
r=json.loads(Path(str(p)+'.result.json').read_text());a=np.load(root/'contact_feasible_start_v1/ground_guide.npz',allow_pickle=True);b=np.load(Path(str(p)+'.npz'),allow_pickle=True)
diff={k:float(np.max(abs(a[k]-b[k]))) for k in ('qfull','ctrl','neuro')};assert max(diff.values())==0. and r['assessment']['accepted']
assert r['assessment']['knee_extension_gate'] and r['assessment']['foot_clearance_gate']
air=json.loads((root/'zero_hyperextension_air_v1/air_pair/SUMMARY.json').read_text());assert air['accepted']
checkpoint=dict(machine='easteregg2',session=1,configuration=str(Path(str(p)+'.json')),result=str(Path(str(p)+'.result.json')),
    measured_assistance_bw=r['ground_contacts']['mean_harness_vertical_bw'],sustained_swings=[f['sustained_swings'] for f in r['foot_clearance']['feet']],
    maximum_knee_extension_deg=r['assessment']['maximum_knee_extension_deg'],fresh_repeat_max_abs_diff=diff,air_pair_summary=str(root/'zero_hyperextension_air_v1/air_pair/SUMMARY.json'),
    scope='Supported stepping in place, corrected ground-aligned initial pose and spring references, protective knee stops, sparse backend. Same neural weights/gains pass air stages1/2; standing notpassed. No completed curriculum, basin robustness, or harness-free claim.',
    supersedes='Do not select earlier43.61%-assistance result: it used a deeply penetrating initial foot pose. Lower-support aligned trials fail geometric clearance.')
snapshot=folder/'source_snapshot';snapshot.mkdir(exist_ok=True)
weights_snapshot=folder/'weights_snapshot'
if not weights_snapshot.exists():shutil.copytree(r['configuration']['weights_directory'],weights_snapshot)
force_snapshot=folder/'shared_force_law_v3.json';shutil.copyfile(r['configuration']['force_law_parameters_file'],force_snapshot)
packaged={**r['configuration'],'weights_directory':str(weights_snapshot),'force_law_parameters_file':str(force_snapshot)}
(folder/'packaged_configuration.json').write_text(json.dumps(packaged,indent=2))
checkpoint['packaged_configuration']=str(folder/'packaged_configuration.json')
for name in ('controller.py','physical_eval.py','balance.py','bilateral.py','contact_events.py','regional_contacts.py','foot_clearance_check.py','verify_bilateral_contract.py'):shutil.copyfile(H/name,snapshot/name)
physics=snapshot/'physics_fidelity';physics.mkdir(exist_ok=True)
for name in ('ground_alignment.py','standing_solver.py','contact_standing_solver.py'):shutil.copyfile(H.parent/'physics_fidelity'/name,physics/name)
checkpoint['current_source_sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in snapshot.iterdir() if f.is_file()}
(folder/'CHECKPOINT.json').write_text(json.dumps(checkpoint,indent=2))
print(json.dumps(checkpoint,indent=2))
note='''

### Corrected-pose checkpoint (supersedes the unaligned lower-support claim)

Fresh repeat: aligned_support_ladder_v1/repeat_0p5, same saved controller as
protected-knee air passers. All current ground gates pass, including actual
foot-mesh clearance and no positive knee angle. Measured assistance63.73%BW;
sustained swings5right/4left; maximum knee angles-1.201/-1.460degrees.
Recorded qfull/ctrl/neuro are bit-identical to the first corrected-pose run.
The initial feet are aligned near the floor before building rig references.
This also rebases the existing ankle spring rest under the runner convention.

The earlier43.61%assistance run used the bad deeply penetrating initial pose;
it is NOT the current checkpoint. After correcting the pose, lower prescribed
support.3/.1 produces measured.556/.598BW assistance and ZERO qualifying
swings per foot. Those runs are rejected; the guide takes more load as the
model collapses. No successful further support wean is claimed.

Corrected-pose standing with the original coupled initializer improves to
12.69degrees tilt and93.6mm sway, load.477BW, but still fails the60mm/.5BW
limits. The contact-aware initializer succeeds in constant-activation holds
but fails in the actual SNS interface (about80mm sway/.17BWfootload), and its
ground test also fails phase criteria. It is an experimental opt-in variant,
not the selected checkpoint. Next investigate delivered MN activation versus
initializer targets, gait-trained reflex activity in standing, and stance
support capacity. Do not sweep insensitive gains or count static holds as SNS
curriculum passes. No stage5 advancement yet.

Full details: SESSION1_PHYSICS_GROUND_CHECKPOINT_20261030.md. Canonical ground
checkpoint: tonic_0p02_exact/aligned_support_ladder_v1/CHECKPOINT.json.
'''
for filename in ('TRAINING_RECORD.md','SESSION1_PHYSICS_GROUND_CHECKPOINT_20261030.md'):
 f=H/filename;s=f.read_text(encoding='utf-8')
 if '### Corrected-pose checkpoint (supersedes the unaligned lower-support claim)' not in s:
  tmp=f.with_suffix('.writing');tmp.write_text(s+note,encoding='utf-8');tmp.replace(f)
