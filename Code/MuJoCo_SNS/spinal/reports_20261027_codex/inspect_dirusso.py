import os
import sys
from pathlib import Path
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
import numpy as np
import mujoco
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'dirusso_gait2392'))
import dirusso_run as D
import controller as C
model = mujoco.MjModel.from_xml_path(D.XML)
data = mujoco.MjData(model)
mujoco.mj_resetDataKeyframe(model, data, 0)
mujoco.mj_forward(model, data)
ref = np.load(D.BSOLVE, allow_pickle=True)
print('model', model.nq, model.nv, model.nu, 'qpos0', data.qpos[:12])
print('ref', {k: ref[k].shape for k in ref.files})
print('ref t', ref['t'][[0,1,-1]], 'q0', ref['qpos'][0,:12])
print('com0', data.subtree_com[0], 'ref heights', ref['qpos'][:3,:7])
for j in range(model.njnt):
    name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT,j)
    if model.jnt_type[j] != mujoco.mjtJoint.mjJNT_FREE:
        print(name, 'qadr',model.jnt_qposadr[j], 'axis', data.xaxis[j])
print('contacts', data.ncon)
for i in range(data.ncon):
    f = np.zeros(6)
    mujoco.mj_contactForce(model,data,i,f)
    print('contact',data.contact[i].frame.reshape(3,3),f)
print('fit', D.stage1_fit()[1])
