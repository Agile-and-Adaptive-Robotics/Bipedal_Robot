"""Actual CAD derivative/power parity; prescribed tensions are not BPA tests."""
import sys,os,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROBOT=HERE.parent
sys.path.insert(0,str(ROBOT));os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,mujoco
from dynamic_routes import DynamicRoutes
from cad_moving_force import CADMovingForce

def main():
 root=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa');m=mujoco.MjModel.from_xml_path(str(ROBOT/'ah001_neural.xml'));d=mujoco.MjData(m);dynamic=DynamicRoutes(m,json.loads((ROBOT/'robotbody_routes.json').read_text()));q=np.load(ROBOT/'runs/verified65_exact_fast.npz')['qpos'];names=[m.actuator(i).name for i in range(m.nu)];m.actuator_gainprm[:]=0.;m.actuator_biasprm[:]=0.;bridge=CADMovingForce(m,dynamic,names);rng=np.random.default_rng(100031);cases=[]
 for index in (0,len(q)//2,len(q)-1):
  base=q[index].copy();d.qpos[:]=base;d.qvel[:]=rng.normal(0,.1,m.nv);bridge.update_geometry(d);jac=dynamic.length_jacobian[bridge.tids].copy();native=d.ten_J[bridge.tids].copy();masks=[a.copy() for a in dynamic.active];fd=np.zeros_like(jac);valid=np.ones(m.nu,dtype=bool)
  for j in range(m.nv):
   direction=np.zeros(m.nv);direction[j]=1.;lengths=[]
   for sign in (1.,-1.):
    d.qpos[:]=base;mujoco.mj_integratePos(m,d.qpos,direction,sign*1e-6);dynamic.update(d);lengths.append(d.ten_length[bridge.tids].copy())
    switched=np.array([not np.array_equal(a,b) for a,b in zip(masks,dynamic.active)]);changed_tids={dynamic.tendon_ids[i] for i in np.flatnonzero(switched)};valid&=np.array([t not in changed_tids for t in bridge.tids])
   fd[:,j]=(lengths[0]-lengths[1])/2e-6
  d.qpos[:]=base;bridge.update_geometry(d);tension=rng.uniform(10,100,m.nu);d.qfrc_applied[:]=0.;actual=bridge.apply_tensions(d,tension);power=float(actual@d.qvel);expected=float(-tension@(jac@d.qvel));error=float(np.max(abs(jac[valid]-fd[valid])));native_error=float(np.max(abs(native[valid]-fd[valid])))
  cases.append(dict(pose_index=index,dofs=m.nv,actuators=m.nu,smooth_routes=int(sum(valid)),conditional_switch_routes=[names[i] for i in np.flatnonzero(~valid)],full_jacobian_error_m_per_rad_or_m=error,fixed_site_only_error=native_error,power_error_W=abs(power-expected),generalized_force_error=float(np.max(abs(actual+jac.T@tension)))))
  if error>2e-6 or abs(power-expected)>1e-10:raise AssertionError(cases[-1])
 report=dict(machine='easteregg2',passed=True,cases=cases,scope='Actual AH001 geometry, three historical poses, independent generalized-coordinate derivatives. Conditional branch switches excluded explicitly from smooth derivative checks. Prescribed synthetic10–100N tensions test force direction and power only; no pressure dynamics, actuator sizing, mounted hardware, standing, or locomotion claim. Native gains/biases disabled globally; complete Jacobian applied once, so no extra native correction call is allowed.')
 (root/'CAD_MOVING_FORCE_PARITY.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
