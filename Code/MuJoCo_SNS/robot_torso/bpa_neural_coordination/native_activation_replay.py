"""Native activation-only replay from actual release state and full control log."""
import json,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import gait_standing_family as F

def replay(prefix):
 r=json.loads(prefix.with_suffix('.result.json').read_text())
 if r.get('bpa_configuration'):raise ValueError('Native activation is not the BPA force state')
 with np.load(prefix.with_suffix('.release_state.npz')) as z:initial=z['act'].copy()
 with np.load(prefix.with_suffix('.npz')) as z:commands=z['ctrl'][-r['completed_samples']:].copy()
 model=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));force=json.loads(Path(r['force_law']).read_text())
 for row in force['muscles'].values():
  for name in row['actuator_names']:model.actuator_dynprm[model.actuator(name).id]=row['dynprm']
 if model.na!=model.nu or initial.shape!=(model.nu,):raise ValueError('Require one native activation per muscle')
 tree=ET.fromstring('<mujoco><option timestep=".002" gravity="0 0 0"/><worldbody><body><joint name="q"/><geom type="sphere" size=".01" mass="1"/></body></worldbody><actuator/></mujoco>')
 for i in range(model.nu):ET.SubElement(tree.find('actuator'),'general',dict(joint='q',name=model.actuator(i).name,dyntype='muscle',dynprm=' '.join(format(x,'.17g') for x in model.actuator_dynprm[i]),gainprm='0',biastype='none',ctrllimited='true',ctrlrange='0 1'))
 m=F.mujoco.MjModel.from_xml_string(ET.tostring(tree,encoding='unicode'));d=F.mujoco.MjData(m);d.act[:]=initial;values=[]
 for command in commands:
  d.ctrl[:]=command;F.mujoco.mj_forward(m,d);d.act[:]+=m.opt.timestep*d.act_dot;np.clip(d.act,0,1,out=d.act);values.append(d.act.copy())
 values=np.asarray(values);error=None
 if prefix.with_suffix('.motor_state.npz').exists():
  with np.load(prefix.with_suffix('.motor_state.npz')) as z:error=float(np.max(abs(z['activation']-values)))
  if error>1e-12:raise ValueError('Prepared native activation reconstruction does not match actual telemetry: '+str(error))
 return values,dict(method='Actual native activation kernel replay from measured full release activation and actual full excitation controls; no legacy act-as-activation substitution',actual_telemetry_error=error)
