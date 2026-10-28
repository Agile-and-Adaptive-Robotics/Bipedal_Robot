"""Check primitive inputs, anatomical perturbations, and exact backend parity."""
import os,sys,json,copy
from pathlib import Path
H=Path(__file__).resolve().parent;S=H.parent
sys.path[:0]=[str(S),str(H),str(S/'w2l_cpg')]
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
import bsolve_ik as B
import network as F
from preallocated_exact_eval import PreallocatedExactForward
from mechanics import CoupledArms,JOINTS
O=S/'reports_20261028_functional_pf'
m,d=B.build_air_model(False)
acts=[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(m.nu)]
net=F.build(acts,options=dict(axes=True,biarticular=True,primitive=True))
ref=copy.deepcopy(net.compiled);net.compiled=PreallocatedExactForward(net.compiled);net.bind_mujoco(m,d)
ref.g_max_non[:]=net.compiled.g_max_non
rng=np.random.default_rng(28);error=0.
for i in range(200):
    u=net.make_inputs();u[:]=rng.uniform(0,.2,len(u));u[net.ports['DRIVE']]=2.
    net.step(u);ref.forward(u)
    error=max(error,float(np.max(abs(ref.V-net.compiled.V))))
assert error==0.
u=net.make_inputs();u[net.ports['DRIVE']]=0.;net.step(u)
primitive_ports=[i for key,i in net.ports.items() if key.startswith('FUNC_BI_')]
assert np.max(abs(u[primitive_ports]))==0.
# Confirm frontal-plane signs by foot motion, independently of joint labels.
arm=CoupledArms(m);d.qpos[arm.addresses]=0.;arm.project(d);mujoco.mj_forward(m,d)
base=d.qpos.copy();signs={}
for side in ('r','l'):
    foot=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'calcn_'+side)
    pelvis=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'pelvis')
    values=[]
    for angle in (-.001,.001):
        d.qpos[:]=base;d.qpos[arm.addresses[JOINTS.index('hip_adduction_'+side)]]=angle
        arm.project(d);mujoco.mj_forward(m,d)
        values.append(float(-(d.xpos[foot,1]-d.xpos[pelvis,1])))
    derivative=(values[1]-values[0])/.002
    signs[side]=dict(rightward_foot_displacement_per_positive_radian_m=derivative,
        positive_coordinate_action='abduction' if derivative*(1 if side=='r' else -1)>0 else 'adduction')
result=dict(neurons=len(net.compiled.V),inputs=len(net.inputs),reference_fast_max_voltage_difference=error,
            primitive_drive_zero_at_zero_descending_drive=True,anatomical_perturbations=signs,
            scope='200 neural steps at fixed mechanics, not full-trajectory equivalence for the dynamic controller')
(O/'primitive_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
