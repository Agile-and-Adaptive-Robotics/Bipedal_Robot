"""Saved-pose force-cone diagnosis, no dynamics/newstance/fitting."""
import sys,json,numpy as np
from scipy.optimize import linprog
sys.path.insert(0,r'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination')
import cad_rigid_stance as C
b=C.Native();z=C.P.load(C.OUT/'REFINED_PRELOAD_STATE.npz');m,d=b.m,b.d;d.qpos[:]=z['qpos'];d.qvel[:]=0;b.s.dynamic.update(d);N=C.R.anatomical_knee_tangent(m,d);T=np.zeros((m.nv,6+N.shape[1]));T[:6,:6]=np.eye(6);T[6:,6:]=N;A=T.T@b.s.dynamic.moment_arms[b.s.tids].T;Qrig=np.zeros(m.nv);C.mujoco.mj_applyFT(m,d,[0,.85*b.weight,0],[0,0,0],d.xipos[b.pelvis],b.pelvis,Qrig);points=[];mus=[];jac=[]
for contact in d.contact:
 names=[m.geom(int(g)).name for g in (contact.geom1,contact.geom2)]
 if 'ground' not in names or not any(name.startswith('foot_') for name in names):continue
 body=int(m.geom_bodyid[contact.geom2 if names[0]=='ground' else contact.geom1]);J=np.zeros((3,m.nv));ang=np.zeros_like(J);C.mujoco.mj_jac(m,d,J,ang,contact.pos,body);points.append(contact.pos.copy());jac.append(T.T@J.T);mus.append(float(contact.friction[0])/np.sqrt(2))
out={}
for mode in ('normal_only','full_friction','normal_only_zero_floor_DIAGNOSTIC'):
 normal=mode.startswith('normal_only');columns=[A]+[J[:,1:2] if normal else J for J in jac];E=np.column_stack(columns);bounds=[(0 if 'zero_floor' in mode else .05*cap,cap) for cap in b.s.caps]
 for i in b.s.ids:bounds[i]=(0,float(b.s.caps[i]))
 bounds+=([(0,b.weight)]*len(jac) if normal else [(None,None),(0,b.weight),(None,None)]*len(jac));ineq=[]
 if mode=='full_friction':
  for i,mu in enumerate(mus):
   for axis in (0,2):
    for sign in (-1,1):row=np.zeros(E.shape[1]);row[m.nu+3*i+axis]=sign;row[m.nu+3*i+1]=-mu;ineq.append(row)
 sol=linprog(np.r_[1/b.s.caps,np.zeros(E.shape[1]-m.nu)],A_eq=E,b_eq=T.T@(d.qfrc_bias-Qrig),A_ub=np.array(ineq) if ineq else None,b_ub=np.zeros(len(ineq)) if ineq else None,bounds=bounds,method='highs');row=dict(feasible=bool(sol.success),message=sol.message)
 if sol.success:row.update(residual_max=float(np.max(np.abs(E@sol.x-T.T@(d.qfrc_bias-Qrig)))),contact_forces=sol.x[m.nu:].reshape(-1,1 if normal else 3).tolist())
 out[mode]=row
 if mode=='full_friction' and sol.success:
  n=E.shape[1];count=2*len(jac);augE=np.column_stack((E,np.zeros((E.shape[0],count))));constraints=[np.r_[row,np.zeros(count)] for row in ineq];augbounds=bounds+[(0,None)]*count
  for i in range(len(jac)):
   for k,axis in enumerate((0,2)):
    for sign in (-1,1):row=np.zeros(n+count);row[m.nu+3*i+axis]=sign;row[n+2*i+k]=-1;constraints.append(row)
  minimum=linprog(np.r_[np.zeros(n),np.ones(count)],A_eq=augE,b_eq=T.T@(d.qfrc_bias-Qrig),A_ub=np.array(constraints),b_ub=np.zeros(len(constraints)),bounds=augbounds,method='highs');out['minimum_total_absolute_shear_N']=float(minimum.fun) if minimum.success else None
b.force[:]=z['force_N'];b.actual();native=[]
for i,contact in enumerate(d.contact):
 names=[m.geom(int(g)).name for g in (contact.geom1,contact.geom2)]
 if 'ground' not in names or not any(name.startswith('foot_') for name in names):continue
 force=np.zeros(6);C.mujoco.mj_contactForce(m,d,i,force);world=contact.frame.reshape(3,3).T@force[:3];native.append(dict(geoms=names,world_force_N=world.tolist()))
out['native_contacts']=native;out['native_total_absolute_shear_N']=sum(abs(row['world_force_N'][0])+abs(row['world_force_N'][2]) for row in native);out['floor_ablation_scope']='Capacity diagnostic only; no physical case, floor policy, neural gains, or source model changed.'
out['scope']='Idealrank4capacitydiagnostic at savedrefinedpose, not nativeequilibrium or anotherphysicalrun';(C.OUT/'FORCE_CONE_DIAGNOSTIC.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
