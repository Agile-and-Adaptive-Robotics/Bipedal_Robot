"""Isolated source-preserving newer W2L stance IK and recruitment projection."""
import os, sys, json, time
from pathlib import Path
os.environ.setdefault('CONDA_PREFIX',r'D:\Anaconda\envs\myo')
import numpy as np
import mujoco
from scipy.optimize import least_squares, lsq_linear
BASE=Path(__file__).resolve().parent
VW=BASE.parent/'virtual_walker_20261031'
SOURCE=VW/'runs/W2L_VirtualWalker_Newer/canonical108_ground_baseline_v1'
OUT=BASE/'runs/canonical_12II/ik_reference'
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(VW))
import native_kernel as K
m=mujoco.MjModel.from_xml_path(str(SOURCE/'model.xml')); d=mujoco.MjData(m)
name=lambda kind,s:mujoco.mj_name2id(m,kind,s)
gid=lambda s:name(mujoco.mjtObj.mjOBJ_GEOM,s)
bid=lambda s:name(mujoco.mjtObj.mjOBJ_BODY,s)
jid=lambda s:name(mujoco.mjtObj.mjOBJ_JOINT,s)
jointnames=['hip_L','knee_L','ankle_L','hip_R','knee_R','ankle_R','toe_L','toe_R']
js=np.array([jid(s) for s in jointnames]); qa=m.jnt_qposadr[js]; va=m.jnt_dofadr[js]
feet=[gid('geom_'+s) for s in ['foot_L','toe_L','foot_R','toe_R']]
rail=[gid('geom_'+s+'_support_rail') for s in ['Left','Right']]
wings=[gid('geom_'+s+'_hip_wing') for s in ['Left','Right']]
floor=gid('geom_Treadmill_Belt'); root=bid('Root')
qinitial=m.qpos0.copy()
signs=np.array([[a,b,c] for a in (-1,1) for b in (-1,1) for c in (-1,1)])
def corners(g): return d.geom_xpos[g]+(signs*m.geom_size[g])@d.geom_xmat[g].reshape(3,3).T
def metrics():
 top=corners(floor)[:,1].max(); parts=[]
 for g in feet:
  p=corners(g); R=d.geom_xmat[g].reshape(3,3)
  # physical main-foot underside is local -z; toe underside is local -y.
  local=(signs*m.geom_size[g]); ix=np.where(signs[:,2]==1)[0] if g in feet[::2] else np.where(signs[:,1]==-1)[0]
  # choose the physically lower of the two faces by their mean world y.
  axis=2 if g in feet[::2] else 1
  a=np.where(signs[:,axis]==1)[0]; b=np.where(signs[:,axis]==-1)[0]; ix=a if p[a,1].mean()<p[b,1].mean() else b
  sole=p[ix]
  parts.append(dict(name=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g),min_gap=float(p[:,1].min()-top),sole_y_span=float(np.ptp(sole[:,1])),pitch=float(np.arctan2(R[1,0],R[0,0])),sole=sole.tolist()))
 gaps=[float(corners(w)[:,1].min()-corners(r)[:,1].max()) for w,r in zip(wings,rail)]
 knee=[]
 for s in ('L','R'):
  F=d.xmat[bid('femur_'+s)].reshape(3,3);T=d.xmat[bid('tibia_'+s)].reshape(3,3)
  downf=-F[:,1]; downt=-T[:,2]; knee.append(float(np.degrees(np.arccos(np.clip(downf@downt,-1,1)))))
 return dict(parts=parts,wing_rail_gaps=gaps,knee_unsigned_flex_deg=knee,root_y=float(d.qpos[1]))
mujoco.mj_forward(m,d)
initial=metrics()
print('INITIAL',json.dumps(initial),flush=True)
traj=np.load(SOURCE/'trajectory.npz');print('NPZ',[(k,traj[k].shape) for k in traj.files],flush=True)
print('JOINTLIMITS',dict(zip(jointnames,np.degrees(m.jnt_range[js]).tolist())),flush=True)
print('INITIALLEG',[(s,d.xmat[bid(s)].reshape(3,3).tolist()) for s in ['femur_L','tibia_L','foot_L','femur_R','tibia_R','foot_R']],flush=True)
np.savez(OUT/'initial_data.npz',qpos=qinitial)
lo=np.r_[.85,m.jnt_range[js[:6],0],[-np.pi/6]*2]
hi=np.r_[1.1,m.jnt_range[js[:6],1],[np.pi/6]*2]
# Both source knee coordinates increase physical flexion. Source knee minima
# happen to imply ~3 deg, already stricter than the requested 2 deg.
lo[2]=max(lo[2],np.deg2rad(2.-10.));lo[5]=max(lo[5],np.deg2rad(2.-55.))
def put(x):
 d.qpos[:]=qinitial;d.qpos[1]=x[0];d.qpos[qa]=x[1:];d.qvel[:]=0;mujoco.mj_forward(m,d)
def residual(x,rail_target=None,exact=False):
 put(x);z=metrics();r=[]
 for i,p in enumerate(z['parts']):
  sole=np.array(p['sole'])
  if exact and i%2:r.append(p['min_gap']/.001)
  else:r.extend((sole[:,1]-(corners(floor)[:,1].max()))/.001)
 if rail_target is not None:r.extend((np.array(z['wing_rail_gaps'])-rail_target)/.001)
 # Weak preference for comparable anatomical thigh and knee poses.
 r.extend(.005*np.array([(-10+np.degrees(x[1]))-(20+np.degrees(x[4])),(10+np.degrees(x[2]))-(55+np.degrees(x[5]))]))
 r.extend(.0001*(x[1:]-np.r_[np.deg2rad([20,0,0,-10,-45,0]),[0,0]]))
 return r
results={}
start=np.r_[.98,np.deg2rad([20,0,0,-10,-45,0]),[0,0]]
for label,target in [('flat_free_height',None),('flat_existing_rails',0.),('exact_mainsole_existing_rails',0.)]:
 tic=time.perf_counter();fit=least_squares(lambda x:residual(x,target,label.startswith('exact')),np.clip(start,lo+1e-8,hi-1e-8),bounds=(lo,hi),max_nfev=1000,xtol=1e-12,ftol=1e-12,gtol=1e-10)
 put(fit.x);z=metrics();z.update(qpos=d.qpos.tolist(),joint_native_deg=dict(zip(jointnames,np.degrees(fit.x[1:]).tolist())),root_vertical_change=float(fit.x[0]-qinitial[1]),residual_norm=float(np.linalg.norm(fit.fun)),elapsed_s=time.perf_counter()-tic)
 z['anatomical_runner_angles_deg']=dict(zip(jointnames[:6],(np.array([-10,10,0,20,55,0])+np.degrees(fit.x[1:7])).tolist()))
 z['all_sole_plane_error_max_m']=max(abs(np.array(p['sole'])[:,1]-corners(floor)[:,1].max()).max() for p in z['parts'])
 z['required_rail_vertical_shift_m']=float(np.mean(z['wing_rail_gaps']))
 results[label]=z;np.savez(OUT/(label+'.npz'),qpos=d.qpos.copy())
 print(label,json.dumps(z),flush=True)
(OUT/'IK_RESULTS.json').write_text(json.dumps(dict(source=str(SOURCE),initial=initial,results=results),indent=2))
# Static recruitment: actual MuJoCo pulling actuator Jacobians, passive springs,
# and unknown regional normal reactions. These optimized reactions are not
# substituted into the simulator and do not constitute a controller fit.
meta=json.loads((SOURCE.parent/'IMPORT_AUDIT.json').read_text()); bank=meta['muscles'];params=K.bank_parameters(bank)
muscles=[(i,b,name(mujoco.mjtObj.mjOBJ_ACTUATOR,b['name'].replace(' ','_'))) for i,b in enumerate(bank) if b['type']=='LinearHillMuscle']
print('HILL',[(b['name'],params[i].tolist()) for i,b,a in muscles],flush=True)
fitx=np.r_[results['flat_existing_rails']['root_y'],np.deg2rad(list(results['flat_existing_rails']['joint_native_deg'].values()))]
put(fitx);d.ctrl[:]=0.;mujoco.mj_forward(m,d)
print('TOE_GAP_SCAN',flush=True)
for deg in [-10,-5,-3,-2,-1,0,1,2,5,10,20,25,30]:
 d.qpos[qa[-2:]]=np.deg2rad(deg);mujoco.mj_forward(m,d);z=metrics();print(deg,[z['parts'][i]['min_gap'] for i in [1,3]],flush=True)
put(fitx)
robotdofs=np.arange(14);jac=np.zeros((3,m.nv));jacr=jac.copy();points=[];columns=[]
for g in feet+wings:
 p=corners(g);bottom=p[np.argsort(p[:,1])[:4]]
 if g in wings:
  # Actual wing/rail intersection, not the nonoverlapping inner wing edge.
  rp=corners(rail[wings.index(g)]);boty=bottom[:,1].mean()
  xmin=max(bottom[:,0].min(),rp[:,0].min());xmax=min(bottom[:,0].max(),rp[:,0].max())
  zmin=max(bottom[:,2].min(),rp[:,2].min());zmax=min(bottom[:,2].max(),rp[:,2].max())
  bottom=np.array([[x,boty,z] for x in [xmin,xmax] for z in [zmin,zmax]])
 for pt in bottom:
  mujoco.mj_jac(m,d,jac,jacr,pt,int(m.geom_bodyid[g]));columns.append(jac[1,robotdofs].copy());points.append(dict(geom=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g),position=pt.tolist()))
J=np.array(columns).T
ai=np.array([a for i,b,a in muscles]);Am=d.actuator_moment[ai][:,robotdofs].T.copy()
caps=[];passive=[]
for i,b,a in muscles:
 l=float(d.actuator_length[a]/m.actuator_gear[a,0]); l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=params[i]
 tl=max(0.,1-((l-l0)/width)**2);pa=kpe*(l-l0)/(1+kpe/kse);mx=(kpe*(l-l0)+(B+D)*tl)/(1+kpe/kse)
 caps.append(float(np.clip(mx,0,maximum)));passive.append(float(pa))
caps=np.array(caps);BW=m.body_subtreemass[root]*9.81
target=(d.qfrc_bias-d.qfrc_passive)[robotdofs]
Aeq=np.column_stack((Am*caps,J*BW));reg=np.diag(np.r_[np.full(len(ai),.005),np.full(J.shape[1],.001)])
minimum=np.maximum(0,np.array(passive));lower=np.r_[minimum/np.maximum(caps,1e-10),np.zeros(J.shape[1])]
sol=lsq_linear(np.vstack((Aeq,reg)),np.r_[target,np.zeros(Aeq.shape[1])],bounds=(lower,np.r_[np.ones(len(ai)),np.full(J.shape[1],2)]),tol=1e-12,max_iter=500)
force=sol.x[:len(ai)]*caps; normal=sol.x[len(ai):]*BW
recruit=[]
for k,(i,b,a) in enumerate(muscles):
 l=float(d.actuator_length[a]/m.actuator_gear[a,0]);l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=params[i];tl=max(0.,1-((l-l0)/width)**2)
 required=(force[k]*(1+kpe/kse)-kpe*(l-l0))/max(tl,1e-15);sig=(required-D)/B
 voltage=float(A-np.log(1/sig-1)/C) if 0<sig<1 else None
 recruit.append(dict(name=b['name'],force_N=float(force[k]),static_Hill_capacity_N=float(caps[k]),route_length_m=l,tension_length_factor=tl,parallel_passive_equilibrium_N=float(passive[k]),required_stimulus_V=voltage,stimulus_inverse_interior=bool(0<sig<1),sigmoid_fraction=float(sig)))
for p,f in zip(points,normal):p['normal_force_N']=float(f)
controller=K.NativeController(meta['source_snapshot'],bank)
for row,(i,b,a) in zip(recruit,muscles):
 binds=[(pre,co) for pre,target,co in controller.outputs if target==i]
 row['MN_bindings']=[]
 for pre,co in binds:
  v=row['required_stimulus_V'];mode='interior'
  if force[len([r for r in recruit[:recruit.index(row)]])]<1e-6:mode='upper_bound_for_zero_clipped_tension'
  if row['sigmoid_fraction']<=1e-7:v=-.15;mode='inactive_low_voltage_convention'
  poly=np.array(co,dtype=float).copy();poly[-1]-=v if v is not None else 0
  roots=np.roots(np.trim_zeros(poly,'f')) if v is not None else []
  real=[float(r.real) for r in roots if abs(r.imag)<1e-8 and -.2<=r.real<=.05]
  row['MN_bindings'].append(dict(name=controller.names[pre],adapter_coefficients=co.tolist(),target_MN_mV=[r*1000 for r in real],interpretation=mode))
res=Aeq@sol.x-target
projection=dict(method='Bounded static equilibrium with original spatial-route moment Jacobian, native passive springs, full root+leg coordinates, upward normal reactions at physical foot/toe and wing bottom corners. Normal-only frictionless static projection; not measured MuJoCo reactions or an MN network fit.',robot_mass_kg=float(m.body_subtreemass[root]),BW_N=float(BW),residual_by_dof=res.tolist(),recruitment=recruit,contact_points=points,foot_normal_BW=float(sum(p['normal_force_N'] for p in points if 'foot' in p['geom'] or 'toe' in p['geom'])/BW),rail_normal_BW=float(sum(p['normal_force_N'] for p in points if 'wing' in p['geom'])/BW),source_hill_formula='F_eq=clip((Kpe*(L-L0)+TL*(Bsig/(1+exp(C*(A-V)))+D))/(1+Kpe/Kse),0,Fmax); TL=max(0,1-((L-L0)/width)^2)',scope='Static equilibrium approximation of original Hill law, no dynamic Hill state validation, no pretense of fitting expected MN traces.')
(OUT/'STATIC_RECRUITMENT.json').write_text(json.dumps(projection,indent=2));print('RECRUITMENT',json.dumps(projection),flush=True)
scan=[]
for deg in np.linspace(-30,30,61):
 put(fitx);d.qpos[qa[-2:]]=np.deg2rad(deg);mujoco.mj_forward(m,d)
 scan.append(dict(toe_native_deg=float(deg),passive_torque=d.qfrc_passive[va[-2:]].tolist(),gaps=[metrics()['parts'][i]['min_gap'] for i in [1,3]]))
(OUT/'TOE_PASSIVE_SCAN.json').write_text(json.dumps(scan,indent=2))
put(fitx);d.ctrl[ai]=force;mujoco.mj_forward(m,d)
hold=[];cf=np.zeros(6);tic=time.perf_counter();kneemin=np.full(2,np.inf);kneemax=np.full(2,-np.inf)
for step in range(1000):
 mujoco.mj_step(m,d);mujoco.mj_forward(m,d)
 knees=np.degrees(d.qpos[qa[[1,4]]])+[10,55];kneemin=np.minimum(kneemin,knees);kneemax=np.maximum(kneemax,knees)
 if step%10==0:
  footF=0.;railF=0.
  for ci,c in enumerate(d.contact):
   mujoco.mj_contactForce(m,d,ci,cf);f=c.frame.reshape(3,3).T@cf[:3]
   if floor in [c.geom1,c.geom2] and any(g in feet for g in [c.geom1,c.geom2]):footF+=abs(f[1])
   if any(g in rail for g in [c.geom1,c.geom2]) and any(g in wings for g in [c.geom1,c.geom2]):railF+=abs(f[1])
  R=d.xmat[root].reshape(3,3);z=metrics();hold.append(dict(t=float(d.time),root_y=float(d.qpos[1]),root_pitch_deg=float(np.degrees(np.arctan2(R[1,0],R[0,0]))),foot_BW=footF/BW,rail_BW=railF/BW,knee_signed_source_deg=(np.degrees(d.qpos[qa[[1,4]]])+[10,55]).tolist(),toe_deg=np.degrees(d.qpos[qa[-2:]]).tolist(),foot_min_gap=[p['min_gap'] for p in z['parts']]))
np.savez(OUT/'FIXED_FORCE_HOLD.npz',samples=np.array([[r['t'],r['root_y'],r['root_pitch_deg'],r['foot_BW'],r['rail_BW'],*r['knee_signed_source_deg'],*r['toe_deg']] for r in hold]))
np.savez(OUT/'fixed_force_settled_qpos.npz',qpos=d.qpos.copy(),qvel=d.qvel.copy(),actuator_force=force)
(OUT/'FIXED_FORCE_HOLD.json').write_text(json.dumps(dict(scope='2s fixed recruitment forces, not SNS or original dynamic Hill closed loop; belt held at zero speed.',warmed_wall_s=time.perf_counter()-tic,knee_min_over_all_physics_states_deg=kneemin.tolist(),knee_max_over_all_physics_states_deg=kneemax.tolist(),final_qpos=d.qpos.tolist(),final_metrics=metrics(),samples=hold),indent=2))
print('HOLD_FINAL',hold[-1],flush=True)
