"""Native force/power and callback state-isolation checks; not a walker run."""
import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,mujoco
from native_bpa_force import NativeBPAForce

XML='''<mujoco><option gravity="0 0 0"/><worldbody>
<site name="a" pos="-.1 0 0"/><site name="b" pos=".1 0 0"/>
<body><joint name="q" axis="0 0 1"/><geom type="capsule" fromto="0 0 0 0 -.3 0" size=".01" mass="1"/><site name="c" pos="0 -.3 0"/></body></worldbody>
<tendon><spatial name="ta"><site site="a"/><site site="c"/></spatial><spatial name="tb"><site site="b"/><site site="c"/></spatial></tendon>
<actuator><general name="left" tendon="ta" gainprm="1"/><general name="right" tendon="tb" gainprm="1"/></actuator>
<sensor><actuatorfrc name="force_left" actuator="left"/><actuatorfrc name="force_right" actuator="right"/></sensor></mujoco>'''

def main():
 m=mujoco.MjModel.from_xml_string(XML);d=mujoco.MjData(m)
 reference=mujoco.MjModel.from_xml_string(XML);rd=mujoco.MjData(reference)
 spec=dict(diameter_mm=20,resting_length_m=.34,maximum_contraction_fraction=.25,
 nonmuscle_length_m=.02,series_stiffness_N_m=2e5,parallel_count=1,
 pressure_max_kpa=620.,tau_fill_s=.04,tau_vent_s=.06,delay_s=.004,
 provenance='Synthetic interface fixture; dimensions and pressure constants are declared test assumptions, not mounted robot identification.')
 cfg=dict(command_mode='direct_pressure',actuators={n:spec.copy() for n in ('left','right')})
 cases=[];original=m.actuator_gaintype.copy()
 with NativeBPAForce(m,cfg,.002) as adapter:
  for q in (-.2,0.,.2):
   for velocity in (-.3,.3):
    d.qpos[:]=rd.qpos[:]=q;d.qvel[:]=rd.qvel[:]=velocity;mujoco.mj_forward(m,d)
    adapter.precharge(d,np.array([.8,.6]));adapter.advance(d,np.array([.3,.9]));before=adapter.bank.pressure.copy()
    mujoco.mj_forward(m,d);mujoco.mj_forward(m,d)
    assert np.array_equal(before,adapter.bank.pressure),'Repeated forward calls advanced pressure'
    rd.ctrl[:]=-adapter.bank.force;mujoco.mj_forward(reference,rd)
    force_error=float(np.max(abs(d.actuator_force-rd.actuator_force)))
    qforce_error=float(np.max(abs(d.qfrc_actuator-rd.qfrc_actuator)))
    sensor_error=float(np.max(abs(d.sensordata-d.actuator_force)))
    power=float(d.qfrc_actuator@d.qvel)
    virtual_power=float(-adapter.bank.force@d.actuator_velocity)
    assert max(force_error,qforce_error,sensor_error,abs(power-virtual_power))<1e-10
    assert np.all(d.actuator_force<0.) and np.prod(d.actuator_moment[:,0])<0.
    cases.append(dict(q_rad=q,qvel_rad_s=velocity,pressure_kpa=before.tolist(),tension_N=adapter.bank.force.tolist(),signed_native_force_N=d.actuator_force.tolist(),actuator_moment=d.actuator_moment[:,0].tolist(),mechanical_power_W=power,virtual_work_power_W=virtual_power,force_error_N=force_error,generalized_force_error=qforce_error,sensor_error_N=sensor_error))
 assert np.array_equal(original,m.actuator_gaintype)
 assert mujoco.get_mjcb_act_bias() is None
 result=dict(machine='easteregg2',passed=True,cases=cases,
 scope='Six native interface force/power tests with opposing route signs and both motion directions; no walker, neural coordination, BPA dynamics identification, or standing claim.')
 out=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\native_pressure_interface_parity.json')
 out.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
