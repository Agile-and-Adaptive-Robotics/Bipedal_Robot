"""Native-actuator parity for BPA tensile direction, layout, and mechanical power."""
import os,sys,json,hashlib
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from bpa_mujoco import BPAMuscleSystem

class ConstantTension:
    def __init__(self,tension):self.tension=tension
    def force(self,length,velocity):return self.tension

def compare(coefficients,velocities):
    joints=len(coefficients[0]);count=len(coefficients)
    # Nested bodies ensure distinct DoFs; a fixed tendon tests exact signs.
    bodies=''.join('<body><joint name="j%d" axis="0 0 1"/><geom type="sphere" size=".05" mass="1"/>'%j for j in range(joints))+'</body>'*joints
    tendons=''.join('<fixed name="t%d">'%i+''.join('<joint joint="j%d" coef="%s"/>'%(j,c) for j,c in enumerate(row))+'</fixed>' for i,row in enumerate(coefficients))
    actuators=''.join('<general name="a%d" tendon="t%d" gainprm="-1"/>'%(i,i) for i in range(count))
    xml='<mujoco><option gravity="0 0 0"/><worldbody>'+bodies+'</worldbody><tendon>'+tendons+'</tendon><actuator>'+actuators+'</actuator></mujoco>'
    model=mujoco.MjModel.from_xml_string(xml);data=mujoco.MjData(model)
    tensions=np.arange(count)*50.+70.;data.ctrl[:]=tensions;data.qpos[:]=.3;data.qvel[:]=velocities
    mujoco.mj_forward(model,data)
    expected=data.qfrc_actuator.copy();length_velocity=data.actuator_velocity.copy()
    native_power=float(expected@data.qvel);tensile_power=float(-tensions@length_velocity)
    model.actuator_gainprm[:]=0.;mujoco.mj_forward(model,data)
    assert np.max(np.abs(data.qfrc_actuator))==0.
    bank=BPAMuscleSystem({'a%d'%i:ConstantTension(float(t)) for i,t in enumerate(tensions)})
    bank.attach(model);bank.cb(model,data)
    actual=data.qfrc_applied.copy();bank.cb(model,data)
    assert np.max(np.abs(actual-data.qfrc_applied))==0.,'Callback accumulates force across calls'
    np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-13)
    actual_power=float(actual@data.qvel)
    assert abs(actual_power-native_power)<1e-13 and abs(actual_power-tensile_power)<1e-13
    return dict(coefficients=coefficients,velocity_rad_s=velocities,shape=list(data.actuator_moment.shape),
        native_torque_Nm=expected.tolist(),bpa_glue_torque_Nm=actual.tolist(),
        native_power_W=native_power,bpa_glue_power_W=actual_power,tensile_virtual_work_W=tensile_power)

cases=[]
for sign in (-1.,1.):
    for velocity in (-.7,.7):cases.append(compare([[sign*.03]],[velocity]))
# These catch orientation errors when a first-dimension heuristic is ambiguous.
for velocity in ([.4,-.6],[-.4,.6]):cases.append(compare([[-.03,.01],[.02,-.025]],velocity))
# And a normal nonsquare (nu,nv) case.
cases.append(compare([[-.03,.01]],[.4,-.6]))
out=HERE/'runs/cad_parameter_audit_20261031'
result=dict(status='PASS',cases=cases,checks=len(cases),mujoco=mujoco.__version__,
    glue_sha256=hashlib.sha256((HERE.parents[1]/'bpa_mujoco.py').read_bytes()).hexdigest())
(out/'BPA_GLUE_NATIVE_PARITY.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
