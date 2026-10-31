"""Matched air-stage checks for the new protected-knee ground candidate."""
import sys,json,subprocess
from pathlib import Path
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'zero_hyperextension_air_v1';out.mkdir()
c=json.loads((root/'zero_hyperextension_v1/ground_no_comm.json').read_text())
for key in ('stage','support_fraction','force_controlled_support','vertical_stiffness_scale','vertical_damping_multiplier','harness_lowering_m','harness_lowering_start_s','harness_lowering_duration_s','height_offset','restriction','start_hip_deg','start_knee_deg','start_ankle_deg','diagnostic_scope'):c.pop(key,None)
(out/'shared.json').write_text(json.dumps(c,indent=2))
with (out/'air_pair.log').open('w') as log:subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',c['weights_directory'],'--output-directory',str(out/'air_pair'),'--config',str(out/'shared.json')],stdout=log,stderr=subprocess.STDOUT,check=True)
print(json.loads((out/'air_pair/SUMMARY.json').read_text())['accepted'],flush=True)
