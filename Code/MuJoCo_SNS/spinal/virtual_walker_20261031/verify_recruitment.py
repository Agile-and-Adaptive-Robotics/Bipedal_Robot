"""Constructor-only calibration checks; no trial or optimization."""
import json
import numpy as np
from import_virtual_walker import HERE
from native_kernel import NativeController
from recruitment import apply_recruitment
meta=json.loads((HERE/'runs/W2L_VirtualWalker_Newer/IMPORT_AUDIT.json').read_text())
base=NativeController(meta['source_snapshot'],meta['muscles']);default=NativeController(meta['source_snapshot'],meta['muscles'])
default_record=apply_recruitment(default)
assert np.array_equal(base.edges,default.edges) and np.array_equal(base.v,default.v)
candidate=NativeController(meta['source_snapshot'],meta['muscles']);record=apply_recruitment(candidate,2.,0.)
assert len(record['pf_mn_edges'])==12 and len(record['toe_mn_edges'])==2
changed=np.flatnonzero(np.any(base.edges!=candidate.edges,axis=1)).tolist();assert len(changed)==14
out=dict(source=meta['source_snapshot'],neurons=len(base.v),edges=len(base.edges),default_bit_exact=True,changed_edges=changed,candidate=record,scope='Static constructor checks only; no gait, physics, or optimization pass implied.')
(HERE/'runs/RECRUITMENT_CONSTRUCTOR_CHECK.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
