import os
import sys
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import json
import numpy as np
import mujoco
import py_compile
HERE=Path(__file__).resolve().parent
SPINAL=HERE.parent
for path in list(HERE.glob('*.py'))+[SPINAL/'w2l_cpg/exact_fast_eval.py',
                                  SPINAL/'dirusso_gait2392/controller_numba.py']:
    py_compile.compile(str(path),doraise=True)
baseline=np.load(HERE/'trial458_d30.npz')
keys=('q','com','pel','heelF','toeF','rgE','rigFz')
for name in ('repeat_0','repeat_1','repeat_2','canonical_repeat'):
    d=np.load(HERE/(name+'.npz'))
    assert all(np.array_equal(baseline[k],d[k]) for k in keys),name
basin=json.loads((HERE/'full_basin_results.json').read_text())
assert len(basin)==12 and all(x['metrics']['survive'] for x in basin)
model=mujoco.MjModel.from_xml_path(str(SPINAL/'w2l_mujoco/w2l_ground.xml'))
weight=float(model.body_mass.sum()*-model.opt.gravity[2])
source=json.loads((HERE/'plumbing_audit.json').read_text())
assert source['inert_dimensions']==['ia_ank_rg','ia_ank_pf']
author=json.loads((HERE/'author_parameter_inventory.json').read_text())
assert author['total']==256
out=dict(python_sources_compile=True,identical_replays=4,
          full_saved_scalar_basin='12/12',w2l_weight_newtons=weight,
          prior_support_metric_normalization_newtons=411.,
          effective_conductance_dimensions=56,sensor_dimensions=3,
          author_parameter_count=author['total'])
(HERE/'verification.json').write_text(json.dumps(out,indent=2))
print(out)
