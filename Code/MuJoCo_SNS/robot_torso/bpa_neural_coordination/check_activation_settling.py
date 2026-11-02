"""Native muscle activation Euler update parity, without moving a body."""
import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,mujoco

XML='''<mujoco><option timestep=".002" integrator="implicitfast" gravity="0 0 0"/>
<worldbody><body><joint name="q"/><geom type="sphere" size=".01" mass="1"/></body></worldbody>
<actuator><general joint="q" dyntype="muscle" dynprm=".01 .04" gainprm="0" ctrllimited="true" ctrlrange="0 1"/></actuator></mujoco>'''

def main():
 m=mujoco.MjModel.from_xml_string(XML);a=mujoco.MjData(m);b=mujoco.MjData(m);rng=np.random.default_rng(100031);trace=[]
 for u in np.r_[np.zeros(20),np.ones(30),rng.uniform(0,1,150),np.zeros(50)]:
  a.ctrl[:]=b.ctrl[:]=u;mujoco.mj_forward(m,a);a.act[:]+=m.opt.timestep*a.act_dot;np.clip(a.act,0,1,out=a.act);mujoco.mj_step(m,b);trace.append([a.act[0],b.act[0]])
 values=np.array(trace);error=float(np.max(abs(values[:,0]-values[:,1])));assert error<1e-14
 report=dict(machine='easteregg2',passed=True,steps=len(trace),max_activation_error=error,bit_exact=bool(np.array_equal(values[:,0],values[:,1])),scope='Explicit activation-only update from engine act_dot versus native implicitfast muscle activation integration. Zero actuator gain and zero gravity prevent body work; the optional pinned warmup leaves qpos/qvel unchanged.')
 out=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\activation_settling_parity.json');out.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
