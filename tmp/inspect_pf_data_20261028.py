import os,json
from pathlib import Path
import numpy as np
root=Path(r'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal')
for p in [root/'bsolve_out.npz',root/'fsa_results/fsa_backsolve.npz',root/'synergy_basis.npz',*list((root/'gait_refs').glob('ong*.npz'))[:1]]:
    if not p.exists():continue
    d=np.load(p,allow_pickle=True)
    print(str(p),json.dumps({k:dict(shape=d[k].shape,dtype=str(d[k].dtype)) for k in d.files}))
    for key in ('joint_names','muscle_names','key_joints','names','t','time','phase_grid','q_names'):
        if key in d:print(key,repr(d[key])[:700])
