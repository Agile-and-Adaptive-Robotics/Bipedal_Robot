import sys
from import_virtual_walker import mujoco,HERE
path=HERE/'runs'/'W2L_VirtualWalker'/'W2L_VirtualWalker.xml';m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
print('initial ncon',d.ncon)
for c in d.contact:
 print(mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,c.geom1),mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,c.geom2),'gap',c.dist)
print('init acceleration',d.qacc)
