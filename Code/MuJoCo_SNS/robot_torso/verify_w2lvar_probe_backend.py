"""Check accelerated neural stepping on the actual 54-actuator variant."""
import os,sys,copy,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
sys.path[:0]=[str(HERE.parent/'spinal'),str(HERE.parent/'spinal/w2l_cpg')]
import numpy as np
import mujoco
import build_network_w2lvar as variant
from exact_fast_eval import ExactFastForward
model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_abduction_diagnostic_v2.xml'))
acts=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
net=variant.build(acts,dt=.002,interleg=True)
reference=net.compiled;fast=ExactFastForward(copy.deepcopy(reference))
rng=np.random.default_rng(20261028);maximum=0.
for step in range(400):
    currents=net.make_inputs()
    for i,name in enumerate(net.inputs):
        currents[i]=3. if name=='DRIVE' else .15 if name=='POSTURE' else float(rng.uniform(0,.3))
    reference.forward(currents);fast.forward(currents)
    assert np.isfinite(reference.V).all() and np.isfinite(fast.V).all()
    maximum=max(maximum,float(np.max(abs(reference.V-fast.V))))
result=dict(machine='easteregg2',true_neurons=int(reference.V.size),inputs=len(net.inputs),
    steps=400,timestep_s=.002,maximum_voltage_difference=maximum,bit_exact=maximum==0.,
    scope='Backend neural-state parity with seeded inputs; not a full-trajectory replay or gait acceptance test')
assert result['bit_exact']
(HERE/'w2lvar_probe_backend_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
