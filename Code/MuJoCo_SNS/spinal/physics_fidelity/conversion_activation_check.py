"""Test whether the old 5 ms force evaluator actually samples activation transients."""
import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import mujoco,numpy as np
ROOT=Path(__file__).resolve().parents[4]
source=ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
m=mujoco.MjModel.from_xml_path(str(source));m.opt.timestep=.005;m.opt.gravity[:]=0.
d=mujoco.MjData(m);aid=m.actuator('bifemsh_r').id;ad=m.actuator_actadr[aid];d.ctrl[aid]=1.
# Replicate the evaluator's initial step followed by seven sampled-pose steps.
rows=[]
for k in range(8):
 d.qpos[:]=m.qpos0;d.qvel[:]=0.;before=float(d.act[ad]);mujoco.mj_step(m,d)
 rows.append(dict(step=k,activation_used_for_force=before,activation_after_step=float(d.act[ad])))
out=dict(dynprm=m.actuator_dynprm[aid].tolist(),timestep_s=.005,steps=rows,
 scope='Activation-state timing only; reset-to-reference poses do not reproduce every conversion force sample.')
path=Path(__file__).resolve().parent/'runs/audit_20261030/CONVERSION_ACTIVATION_CHECK.json';path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
