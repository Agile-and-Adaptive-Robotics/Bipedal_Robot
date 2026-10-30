import json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent;OUT=HERE/'runs/audit_20261030'
r=json.loads((OUT/'SHARED_FORCE_SURROGATE_BIARTICULAR_V3.json').read_text());rows={}
for name,v in r['muscles'].items():
 p=[p for p in v['points'] if p['pose']['heldout'] and p['activation']<=.05]
 rmse=lambda key:math.sqrt(sum((q[key]-q['opensim_tension_N'])**2 for q in p)/len(p))
 rows[name]=dict(heldout_low_drive_samples=len(p),activation_levels=[.01,.05],baseline_RMSE_N=rmse('baseline_tension_N'),shared_RMSE_N=rmse('joint_tension_N'),maximum_bilateral_force_difference_N=v['joint_fit']['maximum_bilateral_force_difference_N'])
(OUT/'BIARTICULAR_V3_LOW_DRIVE_DIAGNOSIS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
