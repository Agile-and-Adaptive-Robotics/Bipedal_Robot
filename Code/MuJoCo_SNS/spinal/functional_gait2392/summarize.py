"""Compact staged results and the unchanged production-baseline gate."""
import os,sys,json
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import numpy as np
S=Path(__file__).resolve().parent.parent;O=S/'reports_20261028_functional_pf'
old=np.load(S/'reports_20261028_family_audit/human_trials/s3k_baseline_preallocated.npz',allow_pickle=True)
new=np.load(O/'trials/default_s3k_anchor.npz',allow_pickle=True)
comparison={k:bool(np.array_equal(old[k],new[k])) for k in old.files if k in new.files}
assert all(comparison.values()),comparison
(O/'default_baseline_verification.json').write_text(json.dumps(dict(all_common_saved_arrays_bit_exact=True,arrays=comparison),indent=2))
rows=[]
for p in sorted((O/'trials').glob('*.json')):
    d=json.loads(p.read_text());m=d['metrics'];k=m.get('kine') or {}
    rows.append(dict(run=p.stem,stage=d['stage'],variant=d['variant'],t_end=m['t_end'],nan=m['nan'],
        wall_s=d['wall_s'],kine_score=m.get('kine_score'),dx_m=m.get('dx'),cycles=[k.get('n_cycles_r'),k.get('n_cycles_l')],
        foot_load=m.get('bal_foot_load'),rig_share=m.get('bal_rig_share'),fell=m.get('bal_fell')))
(O/'trial_summary.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows[-8:],indent=2))
