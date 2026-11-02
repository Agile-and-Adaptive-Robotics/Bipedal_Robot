"""Read the accepted corrected-pose artifact, never re-execute it."""
from pathlib import Path
import numpy as np,json
R=Path(__file__).resolve().parent.parents[3];H=R/'Code/MuJoCo_SNS/spinal/hunt_training';C=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/aligned_support_ladder_v1/CHECKPOINT.json'
c=json.loads(C.read_text());p=Path(c['configuration']).with_suffix('.npz');z=np.load(p)
print(json.dumps(dict(npz=str(p),fields={k:list(z[k].shape) for k in z.files},first_qfull=z['qfull'][0].tolist() if 'qfull' in z else None),indent=2))
