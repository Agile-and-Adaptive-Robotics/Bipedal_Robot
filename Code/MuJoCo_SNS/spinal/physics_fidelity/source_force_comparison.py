"""Source-force audits, native-contact keyframes, and supported transitions.

Counterfactual force comparisons query native passive force on a separate
native-law model. Source components use actual saved fiber state or a checked
reconstruction. The explicit native-transition mode advances the real model
with source-muscle forces and declared support, without kinematic clamps.
"""
import os,sys,json,argparse
from time import perf_counter
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import numpy as np,mujoco
from tendon_replay_campaign import build,SOURCE,OSIM
from thelen_dynamic import DynamicThelenBank,parameters_from_osim,fiber_components

def compare(folder,out):
 z=np.load(folder/'backend_trace.npz');t=z['t'];q=z['qpos'];a=z['activation'];dt=float(np.median(np.diff(t)));names=z['muscle_names'].tolist()
 cfg=json.loads((folder/'controller.json').read_text());m,d=build(q[0],cfg,1.)
 if m.nq!=m.nv or np.any(np.isin(m.jnt_type,[mujoco.mjtJoint.mjJNT_FREE,mujoco.mjtJoint.mjJNT_BALL])):raise ValueError('Backward-difference reconstruction requires scalar slide/hinge coordinates')
 p=parameters_from_osim(OSIM,names,0.);bank=DynamicThelenBank(p,dt);d.qpos[:]=q[0];mujoco.mj_forward(m,d);bank.initialize(d.actuator_length,a[0]);old=bank.state.copy();components=np.zeros((3,m.nu));rows=[];maxerror=0.;Lnorm=[]
 saved='fiber_state' in z and np.any(z['fiber_state'])
 for i in range(1,len(t)):
  if saved and not (5.<=t[i]<=13. and i%4==0):continue
  d.qpos[:]=q[i];d.qvel[:]=(q[i]-q[i-1])/dt;d.act[:]=a[i-1];mujoco.mj_forward(m,d)
  native=-d.actuator_force.copy();d.act[:]=0.;mujoco.mj_fwdActuation(m,d);nativepassive=-d.actuator_force.copy()
  if saved:state=z['fiber_state'][i];old[:]=z['fiber_state'][i-1];fiber_components(state,old,a[i-1],p,dt,components)
  else:
   old[:]=bank.state;bank.commit(d.actuator_length,a[i-1]);state=bank.state;fiber_components(state,old,a[i-1],p,dt,components)
  if 5.<=t[i]<=13. and i%4==0:
   maxerror=max(maxerror,float(np.max(abs(components.sum(axis=0)-z['force'][i]))));Lnorm.append(np.hypot(state,p[:,3])/p[:,1])
   rows.append((i,native,nativepassive,components.copy()))
 if maxerror>1e-5:raise RuntimeError('Fiber reconstruction/activation convention failed measured-force check: '+str(maxerror))
 ids=np.array([r[0] for r in rows]);native=np.array([r[1] for r in rows]);passive=np.array([r[2] for r in rows]);sourcecomp=np.array([r[3] for r in rows]);source=z['force'][ids];Lnorm=np.array(Lnorm)
 output=dict(source_folder=str(folder),scope='Identical physical qpos, reconstructed scalar-joint qvel, and pre-step activation. Native hypothetical force has no closed-loop feedback or effect on recorded motion.',samples=len(rows),component_force_check_max_error_N=maxerror,source_fiber_state='saved' if saved else 'reconstructed and checked against recorded delivered force',muscles={})
 for fam in ('psoas','iliacus','bifemsh','semimem','semiten','bifemlh','rect_fem','vas_lat','soleus','med_gas','lat_gas','tib_ant'):
  records=[]
  for side in ('r','l'):
   j=names.index(fam+'_'+side);adr=int(m.joint('knee_angle_'+side).qposadr[0]);flex=(q[ids,adr]-q[ids-1,adr])/dt<-.1
   records.append(dict(side=side,mean_native_N=float(native[:,j].mean()),mean_source_N=float(source[:,j].mean()),mean_native_passive_N=float(passive[:,j].mean()),mean_source_passive_N=float(sourcecomp[:,1,j].mean()),mean_native_active_during_knee_flexion_N=float((native[flex,j]-passive[flex,j]).mean()) if np.any(flex) else None,mean_source_active_during_knee_flexion_N=float(sourcecomp[flex,0,j].mean()) if np.any(flex) else None,normalized_fiber_p01_p99=np.percentile(Lnorm[:,j],[1,99]).tolist()))
  output['muscles'][fam]=records
 out.write_text(json.dumps(output,indent=2));print(json.dumps(dict(path=str(out),samples=len(rows),force_check=maxerror,muscles={k:v for k,v in output['muscles'].items() if k in ('bifemsh','semimem','rect_fem','soleus','tib_ant')}),indent=2))

def reference_audit(model_folder,out,reference_override=None,passive_sliding=False):
 """Same-source reference feasibility; optimized reactions are never applied."""
 from scipy.optimize import lsq_linear
 from mechanics import CoupledArms,fill_projection,JOINTS
 from contact_standing_solver import _support_regions,_external
 S=HERE.parent;target_path=S/'hunt_training/runs/targets/targets.npz';archive=np.load(target_path);z={key:archive[key] for key in archive.files};metadata=json.loads((target_path.parent/'metadata.json').read_text());cfg=json.loads((model_folder/'controller.json').read_text());m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));snap=np.load(model_folder/'runtime_model_arrays.npz');validation={name:bool(np.array_equal(getattr(m,name),snap[name])) for name in snap.files if not name.startswith('initial_')}
 if reference_override is not None:
  override=np.load(reference_override)
  for key in ('qpos','qvel','qacc'):z[key]=override[key]
 if not all(validation.values()):raise RuntimeError('Binary model does not match captured live arrays')
 out.mkdir(parents=True,exist_ok=True);d=mujoco.MjData(m);arms=CoupledArms(m);acts=[m.actuator(i).name for i in range(m.nu)];p=parameters_from_osim(OSIM,acts,0.);bank=DynamicThelenBank(p,.002);bw=float(m.body_mass.sum()*9.81);dependent=set(arms.dependent.tolist());drivers=np.array([i for i in range(m.nv) if i not in dependent],dtype=np.int64);E=np.zeros((m.nv,len(drivers)));M=np.zeros((m.nv,m.nv));jp=np.zeros((3,m.nv));jr=np.zeros_like(jp);rows=[];indices=np.linspace(0,len(z['t'])-1,17).astype(int);names=[m.joint(int(m.dof_jntid[i])).name for i in drivers];scales=np.array([1./(bw if m.jnt_type[int(m.dof_jntid[i])]==mujoco.mjtJoint.mjJNT_SLIDE else bw*.5) for i in drivers]);root=np.array([j for j,name in enumerate(names) if name.startswith('pelvis_')]);internal=np.array([j for j,name in enumerate(names) if not name.startswith('pelvis_')]);saved_alpha=[];full_alpha=[]
 def force(alpha):
  bank.initialize(d.actuator_length,alpha);return bank.force.copy(),bool(np.any(bank.bounds))
 def fit(rho,tau,weights,contact):
  alpha=np.full(m.nu,.05);f,bound=force(alpha);reaction=np.zeros(contact.shape[1]);domain=not bound
  for _ in range(5):
   lo=np.maximum(.01,alpha-.005);hi=np.minimum(.98,alpha+.005);fl,bl=force(lo);fh,bh=force(hi);domain=domain and not bl and not bh;derivative=(fh-fl)/(hi-lo)
   if np.any(derivative<-1e-6):raise RuntimeError('Nonmonotone source activation tangent')
   A=rho*np.maximum(derivative,0.)[None,:];design=np.column_stack([A,contact]);rhs=tau-rho@f+A@alpha;reg=np.eye(design.shape[1])*.001;sol=lsq_linear(np.vstack([weights[:,None]*design,reg]),np.r_[weights*rhs,np.zeros(design.shape[1])],bounds=(np.r_[np.full(m.nu,.01),np.zeros(contact.shape[1])],np.r_[np.full(m.nu,.98),np.full(contact.shape[1],1.5)]),method='bvls',tol=1e-8,max_iter=150)
   if not sol.success:raise RuntimeError('Source inverse NNLS failed')
   alpha=sol.x[:m.nu];reaction=sol.x[m.nu:];f,hit=force(alpha);domain=domain and not hit
  residual=rho@f+contact@reaction-tau
  return dict(activation=alpha,force=f,reaction=reaction,residual=residual,domain_valid=domain,weighted_residual_norm=float(np.linalg.norm(weights*residual)))
 families=dict(hip_flexor_mono=('psoas','iliacus'),hip_extensor_mono=('glut_max1','glut_max2','glut_max3'),knee_extensor_mono=('vas_lat','vas_med','vas_int'),knee_flexor_mono=('bifemsh',),HAM=('semimem','semiten','bifemlh'),RF=('rect_fem',),GAS=('med_gas','lat_gas'),SAR=('sar',),GRAC=('grac',),TFL=('tfl',))
 for index in indices:
  d.qpos[:]=z['qpos'][index];arms.project(d);d.qvel[:]=z['qvel'][index];acc=z['qacc'][index].copy()
  for dep,ind,c,offset,qa in zip(arms.dependent,arms.independent,arms.coefficients,arms.offsets,arms.qadr):
   x=d.qpos[qa]-offset;first=c[1]+x*(2*c[2]+x*(3*c[3]+4*x*c[4]));second=2*c[2]+x*(6*c[3]+12*x*c[4]);d.qvel[dep]=first*d.qvel[ind];acc[dep]=first*acc[ind]+second*d.qvel[ind]**2
  d.act[:]=.05;d.xfrc_applied[:]=0.;d.xfrc_applied[m.body('pelvis').id,2]=float(cfg['support_fraction'])*bw;mujoco.mj_forward(m,d);rho=arms.update(d).copy();fill_projection(E,drivers,arms.dependent,arms.independent,arms.coefficients,arms.offsets,d.qpos,arms.qadr);mujoco.mj_fullM(m,M,d.qM);gross=E.T@(M@acc+d.qfrc_bias-d.qfrc_passive-_external(m,d));regions,normal=_support_regions(m,d);rays=[];contact_rows=[]
  for ci in range(d.ncon):
   con=d.contact[ci];g1,g2=int(con.geom1),int(con.geom2)
   if m.geom_type[g1]==mujoco.mjtGeom.mjGEOM_PLANE:foot=g2
   elif m.geom_type[g2]==mujoco.mjtGeom.mjGEOM_PLANE:foot=g1
   else:continue
   body=int(m.geom_bodyid[foot]);owner=body
   while owner and m.body(owner).name not in tuple(b+'_'+s for s in ('r','l') for b in ('talus','calcn','toes')):owner=int(m.body_parentid[owner])
   if not owner:continue
   mujoco.mj_jac(m,d,jp,jr,con.pos,body);velocity=jp@d.qvel;mu=float(con.friction[0]);basis=[normal+mu*np.array([1.,0.,0.]),normal-mu*np.array([1.,0.,0.]),normal+mu*np.array([0.,1.,0.]),normal-mu*np.array([0.,1.,0.])]
   if passive_sliding and np.linalg.norm(velocity[:2])>.01:
    axis=int(np.argmax(abs(velocity[:2])));direction=np.zeros(3);direction[axis]=-mu*np.sign(velocity[axis]);basis=[normal+direction]
   rays.extend([E.T@jp.T@ray*bw for ray in basis]);contact_rows.append(dict(body=m.body(owner).name,dist_m=float(con.dist),friction=mu,point_velocity_mps=velocity.tolist()))
  contact=np.column_stack(rays) if rays else np.zeros((len(drivers),0));internal_fit=fit(rho,z['torque'][index],np.full(len(JOINTS),1./(bw*.5)),np.zeros((len(JOINTS),0)));full_fit=fit(E.T@(-d.actuator_moment.T),gross,scales,contact);saved_alpha.append(internal_fit['activation']);full_alpha.append(full_fit['activation']);contributions={}
  for family,members in families.items():
   ids=[i for i,a in enumerate(acts) if a[:-2] in members];contributions[family]={j:float((rho[row,ids]*internal_fit['force'][ids]).sum()) for row,j in enumerate(JOINTS) if j.startswith(('hip_flexion','knee_angle','ankle_angle'))}
  rows.append(dict(reference_index=int(index),time_s=float(z['t'][index]),strict_knees_valid=bool(all(d.qpos[m.joint('knee_angle_'+s).qposadr[0]]<=0. for s in ('r','l'))),regional_mesh_gaps_m={r['label']:r['gap_m'] for r in regions},contacts=contact_rows,source_saved_internal_fit=dict(domain_valid=internal_fit['domain_valid'],residual_Nm=internal_fit['residual'].tolist(),weighted_residual_norm=internal_fit['weighted_residual_norm'],activation_max=float(internal_fit['activation'].max()),group_torque_Nm=contributions),actual_full_fit=dict(domain_valid=full_fit['domain_valid'],weighted_residual_norm=full_fit['weighted_residual_norm'],root_residual=full_fit['residual'][root].tolist(),internal_residual_Nm=full_fit['residual'][internal].tolist(),optimized_contact_normal_total_bw=float(full_fit['reaction'].sum()),contact_scope='Unknown nonnegative actual-contact friction-pyramid rays; upper-bound force feasibility only,not sticking or forward validation. No reaction is applied.')))
 np.savez_compressed(out/'REFERENCE_SOURCE_RECRUITMENT.npz',indices=indices,actuators=acts,source_static_activation=np.array(saved_alpha),actual_contact_static_activation=np.array(full_alpha));report=dict(scope='Same-source source-curve activation feasibility at17 IK poses. Exact captured live model,source Thelen2003 static fiber equilibrium; no force-velocity/activation-history validation yet. Source lengths,and physical constraints unchanged. No runtime trajectory clamp or virtual reaction applied. Saved-internal-torque comparison retains its archived target even when root path is explicitly retargeted; only the full inverse case recomputes dynamics.',reference_override=None if reference_override is None else str(reference_override),passive_sliding_friction_rays=passive_sliding,binary_model_arrays_exact=validation,source_model=str(model_folder/'runtime_model.mjb'),source_parameters=p.tolist(),actuators=acts,independent_coordinates=names,saved_reference_mean_foot_load_bw=metadata['mean_contact_fraction_bw'],actual_protocol_support_fraction=cfg['support_fraction'],saved_reference_root_qpos=z['qpos'][0,:6].tolist(),reference_file=str(target_path),rows=rows)
 (out/'REFERENCE_SOURCE_FEASIBILITY.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(rows=len(rows),all_model_arrays_exact=all(validation.values()),mesh_gap_range_m=[min(v for row in rows for v in row['regional_mesh_gaps_m'].values()),max(v for row in rows for v in row['regional_mesh_gaps_m'].values())],contact_counts=[len(row['contacts']) for row in rows]),indent=2))

def prepare_contact_reference(model_folder,out):
 """Declare offline root-height diagnostics; anatomical reference is immutable."""
 from scipy.signal import savgol_filter
 import bsolve_ik as B
 from mechanics import CoupledArms
 from contact_standing_solver import _support_regions
 target=HERE.parent/'hunt_training/runs/targets/targets.npz';z=np.load(target);m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);arms=CoupledArms(m);tg,gn,gv=B.read_mot(B.GRF_MOT);force=np.array([B.grf_at(z['t'],tg,gn,gv,pref=p)[0][:,2] for p in ('','1_')]).T;anchor=np.argmax(force,axis=1);loaded=force>50.;shifts=[];gaps=[];t=z['t'];dt=float(np.median(np.diff(t)));vertical=m.joint('pelvis_ty').qposadr[0]
 for k in range(len(t)):
  d.qpos[:]=z['qpos'][k];arms.project(d);mujoco.mj_forward(m,d);regions,normal=_support_regions(m,d);pair=np.array([min(regions[2*s]['gap_m'],regions[2*s+1]['gap_m']) for s in (0,1)]);projection=float(d.xaxis[m.joint('pelvis_ty').id]@normal)
  if abs(projection)<.99:raise ValueError('Explicit vertical adapter required')
  gaps.append(pair);shifts.append((-.001-pair[anchor[k]])/projection)
 shifts=np.array(shifts);smooth=savgol_filter(shifts,11,3);rows={};paths={};out.mkdir(parents=True,exist_ok=True)
 for variant,delta in [('raw',np.zeros_like(shifts)),('exact_loaded_foot',shifts),('smooth_loaded_foot',smooth)]:
  q=z['qpos'].copy();v=z['qvel'].copy();a=z['qacc'].copy();q[:,vertical]+=delta
  if variant!='raw':v[:,vertical]=savgol_filter(q[:,vertical],11,3,deriv=1,delta=dt);a[:,vertical]=savgol_filter(q[:,vertical],11,3,deriv=2,delta=dt)
  anatomy=np.delete(np.arange(m.nq),vertical)
  if not np.array_equal(q[:,anatomy],z['qpos'][:,anatomy]):raise RuntimeError('Anatomical reference changed')
  observed=[];active_contacts=[];speed=[]
  for k in range(len(t)):
   d.qpos[:]=q[k];arms.project(d);d.qvel[:]=v[k];mujoco.mj_forward(m,d);regions,normal=_support_regions(m,d);pair=[];velocity=[]
   for s in (0,1):
    region=min(regions[2*s:2*s+2],key=lambda r:r['gap_m']);pair.append(region['gap_m']);jp=np.zeros((3,m.nv));jr=np.zeros_like(jp);mujoco.mj_jac(m,d,jp,jr,region['point'],region['body']);velocity.append((jp@d.qvel).tolist())
   observed.append(pair);speed.append(velocity);active_contacts.append(int(d.ncon))
  observed=np.array(observed);speed=np.array(speed);source_knees=np.array([m.joint('knee_angle_'+s).qposadr[0] for s in ('r','l')]);rows[variant]=dict(root_shift_range_m=[float(delta.min()),float(delta.max())],loaded_mesh_gap_range_m=[float(observed[loaded].min()),float(observed[loaded].max())],loaded_tangential_point_speed_p95_mps=float(np.percentile(np.linalg.norm(speed[:,:,:2],axis=2)[loaded],95)),strict_knee_invalid_frames=int(np.any(q[:,source_knees]>0.,axis=1).sum()),foot_gaps_m=observed.tolist(),contact_counts=active_contacts,foot_point_velocities_mps=speed.tolist(),root_vertical_velocity_range_mps=[float(v[:,vertical].min()),float(v[:,vertical].max())],root_vertical_acceleration_range_mps2=[float(a[:,vertical].min()),float(a[:,vertical].max())]);path=out/(variant+'.npz');np.savez_compressed(path,t=t,qpos=q,qvel=v,qacc=a,expected_loaded=loaded,expected_raw_normal_force_N=force);paths[variant]=str(path)
 report=dict(scope='Offline declared reference preparation diagnostic only. Original saved target and all anatomical angles/source dimensions unchanged. No candidate controller tracks or clamps these trajectories. Exact and Savitzky-Golay11/3 smoothed vertical pelvis translations align the dominant measured-GRF foot mesh to−1mm; they do not guarantee actual contact,passive friction,or dynamics.',source_reference=str(target),reference_variants=paths,source_raw_GRF_normal_force_threshold_N=50.,anchor_side_index=anchor.tolist(),rows=rows);(out/'CONTACT_REFERENCE_VARIANTS.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:{name:row[name] for name in ('root_shift_range_m','loaded_mesh_gap_range_m','loaded_tangential_point_speed_p95_mps','strict_knee_invalid_frames')} for k,row in rows.items()},indent=2))

def coordinate_stop_force(model,data):
 """Generalized coordinate/tendon-stop force,without altering constraints."""
 force=np.zeros(model.nv)
 if data.nefc:
  kinds=np.asarray(data.efc_type[:data.nefc]);selected=(kinds==mujoco.mjtConstraint.mjCNSTR_LIMIT_JOINT)|(kinds==mujoco.mjtConstraint.mjCNSTR_LIMIT_TENDON);multipliers=np.where(selected,data.efc_force[:data.nefc],0.);mujoco.mj_mulJacTVec(model,data,force,multipliers)
 return force

def load_recruitment(folder):
 """Decode immutable keyframes once,not from ZIP at every physics step."""
 with np.load(folder/'source_recruitment.npz') as archive:return {key:archive[key] for key in archive.files}

def full_native_equilibrium(model_folder,folder,spec_path,out,activation_floor=.01,max_nfev=80,resume_from=None,force_priority=1.,feasibility_first=False,forward_acceleration=False,stationary_load_split=False,root_geometry_prefit=False,root_geometry_only=False,initialize_actual_lumbar=False,capacity_aware_lumbar=False):
 """Full native inverse balance; loaded follower states remain real unknowns."""
 from mechanics import CoupledArms
 from contact_standing_solver import _support_regions,_external
 from scipy.optimize import least_squares
 from scipy.spatial.transform import Rotation
 from thelen_dynamic import ConstitutiveBackend
 code_hash=__import__('hashlib').sha256(Path(__file__).read_bytes()).hexdigest()
 if not np.isfinite(force_priority) or force_priority<=0.:raise ValueError('Static force priority must be finite and positive')
 posture_weight=.01 if feasibility_first else 1.;toe_objective_halfwidth=.015 if feasibility_first else .02;toe_objective_sigma=.001 if feasibility_first else .01
 if (out/'SOURCE_NATIVE_STANCE.json').exists():raise ValueError('Preserve existing equilibrium result')
 out.mkdir(parents=True,exist_ok=True);saved=load_recruitment(folder);qideal=saved['qpos'];p=saved['parameters'];names=saved['muscles'].tolist();enabled=saved.get('force_channel_enabled',np.ones(len(names),bool));old=json.loads((folder/'SOURCE_NATIVE_STANCE.json').read_text());oldcfg=json.loads((Path(old['source_pose'])/'MANIFEST.json').read_text());support=float(oldcfg['support_fraction']);m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);ar=CoupledArms(m);spec=json.loads(spec_path.read_text());drivers=np.r_[np.arange(6),ar.rows];addresses=np.r_[np.arange(6),ar.addresses];dep=ar.dependent;depaddr=np.array([m.jnt_qposadr[m.dof_jntid[k]] for k in dep]);bw=float(m.body_mass.sum()*9.81);feet=[m.body('calcn_'+s).id for s in ('r','l')];force6=np.zeros(6)
 if m.nq!=m.nv or len(set(drivers.tolist()+dep.tolist()))!=m.nv:raise ValueError('Full scalar-coordinate coverage required')
 load_reference=oldcfg.get('stationary_load_split_reference') if stationary_load_split else None
 if stationary_load_split and load_reference is None:raise ValueError('Explicit rolling load-split seed manifest required')
 qseed=qideal.copy();alpha_seed=saved['activation'].copy()
 if (folder/'source_hold.npz').exists() and old.get('forward_hold',{}).get('accepted',False):
  with np.load(folder/'source_hold.npz') as hold:qseed=hold['qpos'][-1].copy()
 resumed=None
 if resume_from is not None:
  resumed=load_recruitment(resume_from)
  if not np.array_equal(resumed['parameters'],p) or resumed['muscles'].tolist()!=names or not np.array_equal(resumed.get('force_channel_enabled',np.ones(len(names),bool)),enabled):raise ValueError('Resume changes source parameters,actuator order,or pruning policy')
  qseed=resumed['qpos'].copy();alpha_seed=resumed['activation'].copy()
 d.qpos[:]=qideal;mujoco.mj_fwdPosition(m,d);target_position=d.xpos[feet].copy();target_rotation=d.xmat[feet].reshape(2,3,3).copy();d.qpos[:]=qseed;ar.project(d);offset=qseed[depaddr]-d.qpos[depaddr]
 low=np.r_[qideal[:3]-[.03,.02,.003],np.deg2rad([-10.]*3)];high=np.r_[qideal[:3]+[.03,.02,.003],np.deg2rad([10.]*3)]
 if stationary_load_split:
  anchor=np.asarray(load_reference['root_anchor_native']);low[:6]=np.maximum(low[:6],anchor-np.r_[[.005]*3,np.deg2rad([2.]*3)]);high[:6]=np.minimum(high[:6],anchor+np.r_[[.005]*3,np.deg2rad([2.]*3)])
 for name in [m.joint(int(j)).name for j in ar.jids]:
  lo,hi=spec['joints'][name]['native_bounds_rad'] if name in spec['joints'] else np.deg2rad([-10.,10.])
  if name.startswith('knee_angle_'):hi=min(hi,np.deg2rad(-2.))
  low=np.r_[low,lo];high=np.r_[high,hi]
 ndriver=len(drivers);ndep=len(dep);nalpha=len(names);low=np.r_[low,np.full(ndep,-.005),np.full(nalpha,activation_floor),-.002,np.deg2rad(-.5)];high=np.r_[high,np.full(ndep,.005),np.full(nalpha,.98),.002,np.deg2rad(.5)];x=np.r_[qseed[addresses],offset,alpha_seed,0.,0.]
 if resumed is not None:
  if 'fit_variables' in resumed:x=resumed['fit_variables'].copy()
  else:
   d.qpos[:]=qseed;mujoco.mj_fwdPosition(m,d);x[-2]=np.mean((d.xpos[feet]-target_position)[:,2]);x[-1]=np.mean([Rotation.from_matrix(d.xmat[b].reshape(3,3)@target_rotation[s].T).as_rotvec()[1] for s,b in enumerate(feet)])
 x=np.clip(x,low+1e-10,high-1e-10);xscale=np.r_[np.array([.01]*3+[.1]*(ndriver-3)),np.full(ndep,.001),np.full(nalpha,.1),.001,.01];bank=DynamicThelenBank(p,float(m.opt.timestep));last={};force_scale=np.where(m.jnt_type[m.dof_jntid]==mujoco.mjtJoint.mjJNT_SLIDE,.1,.01)/force_priority
 if forward_acceleration:
  import copy
  acceleration_model=copy.copy(m);acceleration_data=mujoco.MjData(acceleration_model);acceleration_model.actuator_gaintype[:]=mujoco.mjtGain.mjGAIN_FIXED;acceleration_model.actuator_gainprm[:]=0.;acceleration_model.actuator_biastype[:]=mujoco.mjtBias.mjBIAS_AFFINE;acceleration_model.actuator_biasprm[:]=0.;acceleration_scale=np.where(m.jnt_type[m.dof_jntid]==mujoco.mjtJoint.mjJNT_SLIDE,1e-5,1e-4)
 def residual(value):
  d.qpos[:]=qideal;d.qpos[addresses]=value[:ndriver];ar.project(d);d.qpos[depaddr]+=value[ndriver:ndriver+ndep];alpha=value[ndriver+ndep:ndriver+ndep+nalpha];d.act[:]=alpha;d.qvel[:]=0.;d.qacc[:]=0.;d.xfrc_applied[:]=0.;d.xfrc_applied[m.body('pelvis').id,2]=support*bw;mujoco.mj_inverse(m,d);bank.initialize(d.actuator_length,alpha);f=bank.force.copy();f[~enabled]=0.;raw=-d.actuator_moment.T@f-(d.qfrc_inverse-_external(m,d));target=target_position.copy();target[:,2]+=value[-2];R=Rotation.from_rotvec([0.,value[-1],0.]).as_matrix();rotation=np.concatenate([Rotation.from_matrix((R@target_rotation[s])@d.xmat[b].reshape(3,3).T).as_rotvec() for s,b in enumerate(feet)]);normal=toe=0.
  normal_side=np.zeros(2)
  for ci in range(d.ncon):
   c=d.contact[ci];a,b=int(c.geom1),int(c.geom2)
   if m.geom_type[a]==mujoco.mjtGeom.mjGEOM_PLANE:g=b
   elif m.geom_type[b]==mujoco.mjtGeom.mjGEOM_PLANE:g=a
   else:continue
   name=m.body(int(m.geom_bodyid[g])).name
   if not name.startswith(('calcn_','toes_')):continue
   mujoco.mj_contactForce(m,d,ci,force6);normal+=force6[0];toe+=force6[0] if name.startswith('toes_') else 0.;normal_side[0 if name.endswith('_r') else 1]+=force6[0]
  acceleration_residual=np.empty(0)
  if forward_acceleration:
   acceleration_data.qpos[:]=d.qpos;acceleration_data.qvel[:]=0.;acceleration_data.qacc_warmstart[:]=0.;acceleration_data.act[:]=alpha;acceleration_data.ctrl[:]=alpha;acceleration_data.xfrc_applied[:]=d.xfrc_applied;acceleration_model.actuator_biasprm[:,0]=-f;mujoco.mj_forward(acceleration_model,acceleration_data);last['native_forward_acceleration']=acceleration_data.qacc.copy();acceleration_residual=acceleration_data.qacc/acceleration_scale
  regions,_=_support_regions(m,d);gaps=np.array([r['gap_m'] for r in regions]);ratio=toe/max(normal,1e-12);position=d.xpos[feet]-target;last.update(qpos=d.qpos.copy(),activation=alpha.copy(),force=f.copy(),raw=raw.copy(),gaps=gaps.copy(),ratio=ratio,position=position.copy(),rotation=rotation.copy(),normal_side_BW=normal_side/bw,source_domain=bool(not np.any(bank.bounds[enabled])));ratio_error=np.sign(ratio-.08)*max(0.,abs(ratio-.08)-toe_objective_halfwidth);load_residual=(normal_side/bw-np.asarray(load_reference['normal_target_BW']))/.001 if stationary_load_split else np.empty(0);return np.r_[raw/force_scale,acceleration_residual,load_residual,position.ravel()/.0002,rotation/.001,posture_weight*(d.qpos[ar.addresses]-qideal[ar.addresses])/.05,ratio_error/toe_objective_sigma,np.maximum(0.,-.002-gaps)/.00005,1e-4*alpha,1e3*bank.bounds[enabled]]
 residual(x);before_horizontal=last['raw'][[0,2]].copy();before_length=d.actuator_length.copy();neutral=x.copy();neutral[[0,2]]=m.qpos_spring[[0,2]];residual(neutral);neutral_horizontal=last['raw'][[0,2]].copy();length_change=float(np.max(abs(d.actuator_length-before_length)));root_moment=float(np.max(abs(d.actuator_moment[:,[0,2]])));fixed=[0,2] if np.max(abs(neutral_horizontal))<=1e-7 and length_change<=1e-12 and root_moment<=1e-12 and np.all(m.jnt_stiffness[[0,2]]>0) else []
 root_constraints=dict(analytically_eliminated=bool(fixed),coordinates=[m.joint(k).name for k in fixed],spring_reference=m.qpos_spring[[0,2]].tolist(),raw_force_before_N=before_horizontal.tolist(),raw_force_at_reference_N=neutral_horizontal.tolist(),maximum_MTU_length_change_m=length_change,maximum_horizontal_translation_muscle_derivative=root_moment,scope='Static inverse friction is zero at zero velocity/acceleration. With no horizontal external force and translation-invariant muscles,the nonzero rig springs require their neutral horizontal positions. Full105 force residual is still checked; no runtime root constraint is added.')
 if fixed:x=neutral
 geometry_report=None
 if root_geometry_prefit:
  if not stationary_load_split:raise ValueError('Geometry preconditioning requires explicit stationary split')
  residual(x);root_moment_max=float(np.max(abs(d.actuator_moment[:,:6])))
  if root_moment_max>1e-10:raise ValueError('Root cannot be isolated from source muscle actuation')
  lumbar_columns=[6+j for j,jid in enumerate(ar.jids) if m.joint(int(jid)).name.startswith('lumbar_')] if initialize_actual_lumbar else []
  for k in lumbar_columns:x[k]=qideal[ar.addresses[k-6]]
  geometry_columns=np.array([k for k in list(range(ndriver))+[len(x)-2,len(x)-1] if k not in fixed+lumbar_columns]);geometry_base=x.copy();geometry_h=np.array([1e-8 if k<3 or k==len(x)-2 else 1e-6 for k in geometry_columns])
  def lumbar_capacity():
   from scipy.optimize import lsq_linear
   rho=ar.update(d).copy();required=(ar.E.T@(d.qfrc_inverse-_external(m,d)))[-3:];selected=np.array([j for j,n in enumerate(names) if n.startswith(('ercspn_','intobl_','extobl_'))])
   bank.initialize(d.actuator_length,np.full(nalpha,.01));minimum=bank.force[selected].copy();bank.initialize(d.actuator_length,np.full(nalpha,.98));maximum=bank.force[selected].copy();span=maximum-minimum;A=rho[-3:,selected]*span[None,:];target=required-rho[-3:,selected]@minimum
   solution=lsq_linear(A,target,bounds=(np.zeros(len(selected)),np.ones(len(selected))),method='bvls',tol=1e-10,max_iter=50)
   return A@solution.x-target
  def geometry_objective(y):
   value=geometry_base.copy();value[geometry_columns]=y;residual(value);ratio_error=np.sign(last['ratio']-.08)*max(0.,abs(last['ratio']-.08)-toe_objective_halfwidth)
   capacity=lumbar_capacity()/.01 if capacity_aware_lumbar else np.empty(0)
   return np.r_[last['raw'][:6]/force_scale[:6],(last['normal_side_BW']-np.asarray(load_reference['normal_target_BW']))/.001,last['position'].ravel()/.0002,last['rotation']/.001,ratio_error/toe_objective_sigma,np.maximum(0.,-.002-last['gaps'])/.00005,.01*(d.qpos[ar.addresses]-qideal[ar.addresses])/.05,capacity]
  def geometry_jacobian(y):
   base=geometry_objective(y);matrix=np.empty((len(base),len(y)))
   for j,h in enumerate(geometry_h):
    plus=y.copy();minus=y.copy();plus[j]+=h;minus[j]-=h
    if minus[j]>=low[geometry_columns[j]] and plus[j]<=high[geometry_columns[j]]:matrix[:,j]=(geometry_objective(plus)-geometry_objective(minus))/(2*h)
    elif plus[j]<=high[geometry_columns[j]]:matrix[:,j]=(geometry_objective(plus)-base)/h
    else:matrix[:,j]=(base-geometry_objective(minus))/h
   geometry_objective(y);return matrix
  geometry_started=perf_counter();geometry_fit=least_squares(geometry_objective,x[geometry_columns],jac=geometry_jacobian,bounds=(low[geometry_columns],high[geometry_columns]),x_scale='jac',max_nfev=120,xtol=1e-10,ftol=1e-10,gtol=1e-10);x[geometry_columns]=geometry_fit.x;residual(x)
  geometry_report=dict(root_muscle_moment_max=root_moment_max,root_residual_max_N_or_Nm=float(np.max(abs(last['raw'][:6]))),normal_error_max_BW=float(np.max(abs(last['normal_side_BW']-np.asarray(load_reference['normal_target_BW'])))),foot_position_error_m=float(np.max(abs(last['position']))),foot_rotation_error_rad=float(np.max(abs(last['rotation']))),toe_fraction=float(last['ratio']),nfev=int(geometry_fit.nfev),cost=float(geometry_fit.cost),optimality=float(geometry_fit.optimality),wall_s=perf_counter()-geometry_started,scope='Necessary native root/contact geometry preconditioning with original bounds and fixed loaded follower offsets. No muscles can affect global root wrench; complete105-force solve still follows.')
  geometry_report['actual_lumbar_initialization_only']=initialize_actual_lumbar;geometry_report['full_solver_lumbar_bounds_unchanged']=True
  if capacity_aware_lumbar:geometry_report['necessary_lumbar_capacity_residual_max_Nm']=float(np.max(abs(lumbar_capacity())));residual(x)
  (out/'GEOMETRY_PREFIT.json').write_text(json.dumps(geometry_report,indent=2));print(json.dumps(dict(geometry_prefit=geometry_report),indent=2),flush=True)
  if root_geometry_only:
   seed_manifest=dict(format='aarl_root_geometry_reference_only_v1',support_fraction=support,stationary_load_split_reference=load_reference,source_code_sha256=code_hash,root_geometry_prefit=geometry_report,scope='Root/contact geometry only; source-muscle/follower force equilibrium and five-second validation are NOT passed.')
   (out/'MANIFEST.json').write_text(json.dumps(seed_manifest,indent=2));np.savez_compressed(out/'source_recruitment.npz',qpos=last['qpos'],activation=last['activation'],source_force=last['force'],parameters=p,muscles=names,force_channel_enabled=enabled,fit_variables=x)
   result=dict(source_pose=str(out),accepted_static_source_recruitment=False,accepted_stationary_endpoint=False,root_geometry=geometry_report,scope=seed_manifest['scope']);(out/'SOURCE_NATIVE_STANCE.json').write_text(json.dumps(result,indent=2));return result
 initial=residual(x);initial_cost=float(.5*np.dot(initial,initial));free=np.array([k for k in range(len(x)) if k not in fixed]);fixed_x=x.copy()
 def objective(value):
  full=fixed_x.copy();full[free]=value;return residual(full)
 derivative_step=np.r_[[1e-8]*3,[1e-6]*(ndriver-3),np.full(ndep,1e-8),np.full(nalpha,1e-5),1e-8,1e-6][free]
 def jacobian(value):
  base=objective(value);J=np.empty((len(base),len(value)))
  for column,h in enumerate(derivative_step):
   plus=value.copy();minus=value.copy();plus[column]+=h;minus[column]-=h
   if minus[column]>=low[free[column]] and plus[column]<=high[free[column]]:J[:,column]=(objective(plus)-objective(minus))/(2.*h)
   elif plus[column]<=high[free[column]]:J[:,column]=(objective(plus)-base)/h
   else:J[:,column]=(base-objective(minus))/h
  objective(value);return J
 started=perf_counter();fit=least_squares(objective,x[free],jac=jacobian,bounds=(low[free],high[free]),x_scale='jac',max_nfev=max_nfev,xtol=1e-10,ftol=1e-10,gtol=1e-10);fit_x=fixed_x.copy();fit_x[free]=fit.x;residual(fit_x);elapsed=perf_counter()-started
 static=dict(full_generalized_residual_max_N_or_Nm=float(np.max(abs(last['raw']))),follower_force_residual_max_N=float(np.max(abs(last['raw'][dep]))),driver_force_residual_max_N_or_Nm=float(np.max(abs(last['raw'][drivers]))),foot_position_error_m=float(np.max(abs(last['position']))),foot_rotation_error_rad=float(np.max(abs(last['rotation']))),minimum_mesh_gap_m=float(last['gaps'].min()),toe_normal_fraction=float(last['ratio']),source_domain_valid=last['source_domain'],source_equilibrium_error_N=float(np.max(abs(bank.error[enabled]))),initial_cost=initial_cost,final_cost=float(fit.cost),optimality=float(fit.optimality),nfev=int(fit.nfev),optimizer_status=int(fit.status),wall_s=elapsed)
 accepted=bool(static['full_generalized_residual_max_N_or_Nm']<=.01 and static['foot_position_error_m']<=.0002 and static['foot_rotation_error_rad']<=.001 and static['minimum_mesh_gap_m']>=-.002 and abs(static['toe_normal_fraction']-.08)<=.02 and last['source_domain'] and static['source_equilibrium_error_N']<=1e-5);q=last['qpos'].copy();alpha=last['activation'].copy()
 if stationary_load_split:
  static['normal_side_BW']=last['normal_side_BW'].tolist();static['normal_target_BW']=load_reference['normal_target_BW'];static['normal_target_error_max_BW']=float(np.max(abs(last['normal_side_BW']-np.asarray(load_reference['normal_target_BW']))));accepted=accepted and static['normal_target_error_max_BW']<=.002
 manifest=dict(format='aarl_full_native_static_equilibrium_v2',machine='easteregg2',source_folder=str(folder),resume_from=str(resume_from) if resume_from is not None else None,source_model=str(model_folder/'runtime_model.mjb'),source_code_sha256=code_hash,support_fraction=support,activation_floor=activation_floor,maximum_evaluations=max_nfev,static_force_objective_priority=force_priority,analytic_static_root_constraints=root_constraints,toe_load_objective='Hinge loss outside the unchanged admissible6–10% normal-load interval; no competing exact8% preference inside that interval.',scope='All105 native generalized force residuals,including loaded soft-follower coordinates. No projected cancellation,dimension,contact,stiffness,or source-strength changes. Independent phase activation states; bilateral source parameters preserved.',acceptance=dict(maximum_static_raw_force_residual=.01,forward_duration_s=5.,joint_drift_max_deg=.25,last_half_joint_speed_max_deg_s=.05,last_half_fiber_speed_max_l0_s=.001))
 manifest['numerical_jacobian']=dict(method='Bound-aware central finite differences in physical units; one-sided only at variable bounds.',translation_step_m=1e-8,follower_step_m=1e-8,angle_step_rad=1e-6,activation_step=1e-5,trust_region_variable_scale='Inverse Jacobian column norm; numerical conditioning only.')
 manifest['toe_load_objective']='Hinge loss outside the explicit objective interval below; the admissible6–10% gate is unchanged.'
 manifest['objective_regularization']=dict(feasibility_first=feasibility_first,posture_weight=posture_weight,toe_load_objective_interval=[.08-toe_objective_halfwidth,.08+toe_objective_halfwidth],toe_load_objective_sigma=toe_objective_sigma,toe_load_acceptance_interval=[.06,.10],scope='Only objective preferences change; original target,ROM bounds,and all acceptance gates remain unchanged.')
 manifest['native_forward_acceleration_objective']=dict(enabled=forward_acceleration,translation_scale_mps2=1e-5,angular_scale_radps2=1e-4,scope='Separate native source-force model/data at rest,all105 accelerations,zeroed probe warmstart; no live state or physics changes.')
 if stationary_load_split:manifest['stationary_load_split_reference']=load_reference
 if geometry_report is not None:manifest['root_geometry_prefit']=geometry_report
 if forward_acceleration:static['native_forward_acceleration_max_native']=float(np.max(abs(last['native_forward_acceleration'])))
 static.update(njev=int(fit.njev),optimizer_message=str(fit.message),active_lower_bounds=int(np.sum(fit.active_mask==-1)),active_upper_bounds=int(np.sum(fit.active_mask==1)))
 (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));np.savez_compressed(out/'source_recruitment.npz',qpos=q,activation=alpha,source_force=last['force'],parameters=p,muscles=names,force_channel_enabled=enabled,full_generalized_static_residual=last['raw'],fit_variables=fit_x)
 result=dict(source_pose=str(out),accepted_static_source_recruitment=accepted,accepted_stationary_endpoint=False,forward_validation='Five-second constant-control validation pending',actuation_policy='Original source92' if enabled.all() else 'Saved pruning mask preserved',static_equilibrium=static,scope=manifest['scope']);(out/'SOURCE_NATIVE_STANCE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result

def audit_landing_reference(model_folder,folder,out):
 from contact_standing_solver import _support_regions
 with np.load(folder/'training_reference.npz') as archive:z={key:archive[key] for key in archive.files}
 report=json.loads((folder/'TRANSITION_RESULT.json').read_text());side=report['metrics']['swing_side'];s=0 if side=='r' else 1;air=(z['foot_mesh_gaps_m'][:,s]>.001)&(z['normal_side_BW'][:,s]<.025);last=int(np.flatnonzero(air)[-1]);land=int(np.flatnonzero((np.arange(len(air))>last)&(z['normal_side_BW'][:,s]>.025))[0]);stride=max(1,round(.02/np.median(np.diff(z['t']))));indices=np.unique(np.r_[np.arange(max(0,land-5*stride),len(air),stride),land]);m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);bw=float(m.body_mass.sum()*9.81);m.actuator_gaintype[:]=mujoco.mjtGain.mjGAIN_FIXED;m.actuator_gainprm[:]=0.;m.actuator_biastype[:]=mujoco.mjtBias.mjBIAS_AFFINE;m.actuator_biasprm[:]=0.;force=np.zeros(6);rows=[];max_load_error=0.
 for k in indices:
  d.qpos[:]=z['qpos'][k];d.qvel[:]=z['qvel'][k];d.act[:]=z['activation'][k];d.ctrl[:]=z['excitation'][k];d.qacc_warmstart[:]=z['qacc'][k];d.xfrc_applied[m.body('pelvis').id,2]=report['manifest']['support_fraction']*bw;m.actuator_biasprm[:,0]=-z['source_force_N'][k];mujoco.mj_forward(m,d);regions,_=_support_regions(m,d);normal=np.zeros(2);body_load={}
  for ci in range(d.ncon):
   c=d.contact[ci];a,b=int(c.geom1),int(c.geom2)
   if m.geom_type[a]==mujoco.mjtGeom.mjGEOM_PLANE:g=b
   elif m.geom_type[b]==mujoco.mjtGeom.mjGEOM_PLANE:g=a
   else:continue
   name=m.body(int(m.geom_bodyid[g])).name
   if not name.endswith(('_r','_l')):continue
   mujoco.mj_contactForce(m,d,ci,force);normal[0 if name.endswith('_r') else 1]+=force[0]/bw;body_load[name]=body_load.get(name,0.)+float(force[0]/bw)
  max_load_error=max(max_load_error,float(np.max(abs(normal-z['normal_side_BW'][k]))));axis=d.xmat[m.body('calcn_'+side).id].reshape(3,3)[:,0];rows.append(dict(time_s=float(z['t'][k]),heel_gap_m=float(regions[2*s]['gap_m']),toe_gap_m=float(regions[2*s+1]['gap_m']),foot_forward_axis_pitch_deg=float(np.rad2deg(np.arctan2(axis[2],np.hypot(axis[0],axis[1])))),ankle_native_deg=float(np.rad2deg(d.qpos[m.joint('ankle_angle_'+side).qposadr[0]])),mtp_native_deg=float(np.rad2deg(d.qpos[m.joint('mtp_angle_'+side).qposadr[0]])),body_normal_BW=body_load,stored_normal_BW=z['normal_side_BW'][k].tolist()))
 result=dict(source=str(folder),first_loaded_time_s=float(z['t'][land]),swing_side=side,maximum_reconstructed_normal_load_error_BW=max_load_error,accepted_force_reconstruction=bool(max_load_error<=1e-6),samples=rows,scope='Read-only reconstruction from actual saved full state and delivered source forces. Original native contact geometry/solver/support. Toe/heel regional gaps and actual body contact loads are separate.')
 out.mkdir(parents=True,exist_ok=True);target=out/'LANDING_CONTACT_AUDIT.json'
 if target.exists():raise ValueError('Preserve existing landing audit')
 target.write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='samples'},indent=2))

class NativeAccelerationTeacher:
 """Source recruitment against native floating-root/contact acceleration response."""
 def __init__(self,model,arms,enabled,omega=5.,delta_force=.01):
  import copy
  self.live_model=model;self.model=copy.copy(model);self.data=mujoco.MjData(self.model);self.arms=arms;self.enabled=enabled;self.omega=omega;self.delta=delta_force
  if np.shares_memory(model.actuator_biasprm,self.model.actuator_biasprm):raise ValueError('Native probe model shares mutable force parameters')
  self.rows=np.r_[np.arange(6),arms.rows];self.addresses=np.r_[np.arange(6),arms.addresses];self.E=np.zeros((model.nv,len(self.rows)));self.E[:6,:6]=np.eye(6)
  self.mask=mujoco.mjtState.mjSTATE_INTEGRATION;self.state=np.empty(mujoco.mj_stateSize(model,self.mask));self.M=np.zeros((model.nv,model.nv));self.G=np.zeros((len(self.rows),model.nu));self.skip=True;self.skip_checked=False
  self.updates=0;self.max_clone_error=0.;self.max_skip_error=0.;self.max_linearity_error=0.;self.invalid_maps=0;self.wall=0.;self.skip_checks=0;self.max_impedance_identity_error=0.
  families=sorted({model.actuator(i).name[:-2] for i in range(model.nu)});self.test=np.array([.1*delta_force*np.sin(families.index(model.actuator(i).name[:-2])+1.) if enabled[i] else 0. for i in range(model.nu)])

 def update(self,live,gain,activation,qref,vref,aref,alpha_reference,torque_task=None,torque_rhs=None,joint_impedance=None,source_prediction=None):
  from scipy.optimize import lsq_linear
  start=perf_counter();m,d=self.model,self.data;m.actuator_gainprm[:]=self.live_model.actuator_gainprm;m.actuator_biasprm[:]=self.live_model.actuator_biasprm;bias=m.actuator_biasprm[:,0].copy()
  mujoco.mj_getState(self.live_model,live,self.state,self.mask);mujoco.mj_setState(m,d,self.state,self.mask);mujoco.mj_forward(m,d);base=d.qacc[self.rows].copy();clone_error=float(np.max(abs(base-live.qacc[self.rows])));self.max_clone_error=max(self.max_clone_error,clone_error);actual_base=base.copy();linearization_control=activation
  if source_prediction is not None:
   future_force,future_control,current_force=source_prediction;bias-=future_force-current_force;m.actuator_biasprm[:,0]=bias;mujoco.mj_forward(m,d);base=d.qacc[self.rows].copy();linearization_control=future_control
  def evaluate():
   d.qacc_warmstart[:]=live.qacc_warmstart
   if self.skip:mujoco.mj_forwardSkip(m,d,mujoco.mjtStage.mjSTAGE_VEL,1)
   else:mujoco.mj_forward(m,d)
   return d.qacc[self.rows].copy()
  checked=False
  for k in np.flatnonzero(self.enabled):
   m.actuator_biasprm[k,0]=bias[k]-self.delta
   if not checked:
    d.qacc_warmstart[:]=live.qacc_warmstart;mujoco.mj_forwardSkip(m,d,mujoco.mjtStage.mjSTAGE_VEL,1);fast=d.qacc[self.rows].copy();d.qacc_warmstart[:]=live.qacc_warmstart;mujoco.mj_forward(m,d);full=d.qacc[self.rows].copy();error=float(np.max(abs(fast-full)));self.max_skip_error=max(self.max_skip_error,error);self.skip=error<=1e-7;self.skip_checked=True;checked=True;self.skip_checks+=1;response=full
   else:response=evaluate()
   self.G[:,k]=(response-base)/self.delta;m.actuator_biasprm[k,0]=bias[k]
  self.E[:,6:]=self.arms.E;mujoco.mj_fullM(self.live_model,self.M,live.qM);reduced=self.E.T@self.M@self.E;metric=np.linalg.cholesky(reduced).T
  m.actuator_biasprm[:,0]=bias-self.test;observed=evaluate()-base;m.actuator_biasprm[:,0]=bias;error=float(np.linalg.norm(metric@(observed-self.G@self.test)));scale=float(np.linalg.norm(metric@observed));self.max_linearity_error=max(self.max_linearity_error,error)
  valid=clone_error<=1e-7*(1.+np.max(abs(actual_base))) and error<=1e-4+.05*scale;self.invalid_maps+=int(not valid)
  desired=aref[self.rows]+self.omega**2*(qref[self.addresses]-live.qpos[self.addresses])+2.*self.omega*(vref[self.rows]-live.qvel[self.rows]);impedance_error=0.
  if joint_impedance is not None:
   kp,kd=joint_impedance;feedback=kp*(qref[self.arms.addresses]-live.qpos[self.arms.addresses])+kd*(vref[self.arms.rows]-live.qvel[self.arms.rows]);root_correction=desired[:6]-aref[:6];joint_correction=np.linalg.solve(reduced[6:,6:],feedback-reduced[6:,:6]@root_correction);desired[6:]=aref[self.rows[6:]]+joint_correction;impedance_error=float(np.max(abs(reduced[6:,:6]@root_correction+reduced[6:,6:]@joint_correction-feedback)));self.max_impedance_identity_error=max(self.max_impedance_identity_error,impedance_error);impedance_valid=impedance_error<=1e-8*(1.+np.max(abs(feedback)));self.invalid_maps+=int(valid and not impedance_valid);valid=valid and impedance_valid
  A=self.G*gain[None,:];rhs=desired-base+A@linearization_control;length_scale=None
  if torque_task is None:WA=metric@A;wrhs=metric@rhs
  else:
   length_scale=abs(float(qref[m.joint('pelvis_ty').qposadr[0]]));wrench=np.diag([length_scale]*3+[1.]*3)@reduced[:6,:6];WA=np.vstack([torque_task,wrench@A[:6]]);wrhs=np.r_[torque_rhs,wrench@rhs[:6]]
  norm=np.linalg.norm(WA,axis=0);regularization=.02*np.median(norm[norm>1e-5]);fit=lsq_linear(np.vstack([WA,regularization*np.eye(m.nu)]),np.r_[wrhs,regularization*alpha_reference],bounds=(.01,.98),method='bvls',tol=1e-6,max_iter=100)
  self.updates+=1;self.wall+=perf_counter()-start;info=dict(valid_map=bool(valid),mode='native_root_wrench_plus_joint_torque' if torque_task is not None else ('native_mass_consistent_joint_impedance' if joint_impedance is not None else 'native_mass_weighted_acceleration'),predictive_source_force=source_prediction is not None,virtual_joint_torque_identity_error_Nm=impedance_error,root_force_length_scale_m=length_scale,baseline_clone_max_error=clone_error,linearity_mass_weighted_error=error,linearity_mass_weighted_change=scale,residual_norm=float(np.linalg.norm(WA@fit.x-wrhs)),root_reference_error_native=(qref[:6]-live.qpos[:6]).tolist(),desired_acceleration=desired.tolist(),native_baseline_acceleration=actual_base.tolist(),predicted_baseline_acceleration=base.tolist(),saturated_muscles=[m.actuator(k).name for k in np.flatnonzero(fit.x>.979)])
  return fit.x,info

 def diagnostics(self):
  return dict(updates=self.updates,invalid_maps=self.invalid_maps,maximum_clone_acceleration_error=self.max_clone_error,maximum_forward_skip_error=self.max_skip_error,forward_skip_comparisons=self.skip_checks,maximum_mass_weighted_linearity_error=self.max_linearity_error,maximum_virtual_joint_torque_identity_error_Nm=self.max_impedance_identity_error,force_only_forward_skip=self.skip,map_wall_s=self.wall,omega_rad_s=self.omega,force_probe_N=self.delta,scope='Separate native model/data; exact integration state including warmstart/applied forces. Native contact/soft followers retained; skip/full and combined-force checks each update. Optional mass-consistent joint impedance preserves PD feedback conditional on root acceleration. Only source-muscle excitation is applied to live physics.')

def native_transition(model_folder,start_folder,end_folder,spec_path,out,ramp=1.,warm=.5,hold=.5,stride=10,feedback=False,command_filter=0.,feedback_warmup=False,via_folder=None,mtp_impedance_scaled=False,physics_dt=None,roundtrip=False,return_dwell=0.,export_training_reference=False,landing_duration_scale=1.,cartesian_root=False,stance_posture_cartesian=False,native_acceleration=False,native_root_augment=False,compensate_coordinate_stops=True,world_foot_pitch=False,loaded_cartesian_followers=False,cartesian_tangential_clearance=None,mass_consistent_impedance=False,predictive_source_horizon=None,actuation_audit_times=None,loaded_contact_min_BW=.025):
 if not np.isfinite(loaded_contact_min_BW) or loaded_contact_min_BW<0.:raise ValueError('Finite nonnegative contact audit threshold required')
 """One physical transition between existing feasible source activation states."""
 from thelen_dynamic import ConstitutiveBackend
 import hashlib
 if ramp<=0. or warm<=0. or hold<0. or stride<1 or command_filter<0. or return_dwell<0.:raise ValueError('Invalid transition duration,stride,or excitation-filter time constant')
 if physics_dt is not None and physics_dt<=0.:raise ValueError('Physics timestep must be positive')
 if native_acceleration and native_root_augment:raise ValueError('Choose one native-map recruitment objective')
 if native_acceleration and mtp_impedance_scaled and not mass_consistent_impedance:raise ValueError('Torque-impedance option does not apply to acceleration-only teacher')
 if native_acceleration and not compensate_coordinate_stops and not mass_consistent_impedance:raise ValueError('Stop-compensation option applies to joint-torque teacher only')
 if mass_consistent_impedance and not native_acceleration:raise ValueError('Mass-consistent impedance requires native acceleration mode')
 if predictive_source_horizon is not None and (not feedback or mass_consistent_impedance or predictive_source_horizon<=0. or not np.isfinite(predictive_source_horizon)):raise ValueError('Source prediction requires feedback,positive horizon,and a separately validated recruitment mode')
 if actuation_audit_times and (not feedback or len(actuation_audit_times)>3 or any(not np.isfinite(t) or t<0 for t in actuation_audit_times)):raise ValueError('Actuation audit supports up to3 finite nonnegative feedback snapshots')
 if world_foot_pitch and (not feedback or cartesian_root):raise ValueError('World-foot-pitch teacher requires source feedback and a separate target mode')
 if stance_posture_cartesian and not cartesian_root:raise ValueError('Stance-posture Cartesian option requires root-aware Cartesian targets')
 if loaded_cartesian_followers and not cartesian_root:raise ValueError('Loaded Cartesian followers require root-aware Cartesian targets')
 if cartesian_tangential_clearance is not None and not loaded_cartesian_followers:raise ValueError('Tangential clearance requires loaded Cartesian FK')
 if (out/'TRANSITION_RESULT.json').exists():raise ValueError('Completed transition exists; use a new output folder to preserve evidence')
 out.mkdir(parents=True,exist_ok=True);start=load_recruitment(start_folder);end=load_recruitment(end_folder);a0=start['activation'];a1=end['activation'];p=start['parameters'];names=start['muscles'].tolist();enabled=start['force_channel_enabled'] if 'force_channel_enabled' in start else np.ones(len(names),bool);other=end['force_channel_enabled'] if 'force_channel_enabled' in end else np.ones(len(names),bool)
 if not np.array_equal(p,end['parameters']) or names!=end['muscles'].tolist() or not np.array_equal(enabled,other):raise ValueError('Transition endpoints have different source parameters,actuator order,or pruning policy')
 reports=[json.loads((folder/'SOURCE_NATIVE_STANCE.json').read_text()) for folder in (start_folder,end_folder)];configs=[json.loads((Path(r['source_pose'])/'MANIFEST.json').read_text()) for r in reports];supports=[c['support_fraction'] for c in configs]
 for report,config in zip(reports,configs):
  if config.get('format','').startswith(('aarl_full_native_static_equilibrium_','aarl_full_native_stationarity_validation_')) and not report.get('accepted_stationary_endpoint',False):raise ValueError('New full-native endpoints require accepted five-second stationarity before transition training')
 if supports[0]!=supports[1]:raise ValueError('Support changes must be explicitly modeled,not silently interpolated')
 m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
 if physics_dt is not None:m.opt.timestep=physics_dt
 d=mujoco.MjData(m);dt=float(m.opt.timestep);control_stride=max(1,round(.01/dt));bw=float(m.body_mass.sum()*9.81);spec=json.loads(spec_path.read_text());ref=np.load(spec['source_target']);joint_names=list(spec['joints']);addresses=np.array([spec['joints'][n]['qpos_address'] for n in joint_names]);low=np.array([spec['joints'][n]['native_bounds_rad'][0] for n in joint_names]);high=np.array([spec['joints'][n]['native_bounds_rad'][1] for n in joint_names]);rate_limit=np.max(abs(ref['qvel'][:,addresses]),axis=0);acc_limit=np.max(abs(ref['qacc'][:,addresses]),axis=0);disabled=~enabled;policy='Original source92,including6legacy-pruned shortrotators' if enabled.all() else 'Explicit saved force-channel mask preserved';manifest=dict(format='aarl_native_source_transition_v1',source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),start=str(start_folder),end=str(end_folder),support_fraction=supports[0],actuation_policy=policy,disabled_muscles=[n for n,e in zip(names,enabled) if not e],warm_s=warm,ramp_s=ramp,terminal_hold_s=hold,timestep_s=dt,trace_stride=stride,control='C2 muscle-excitation reference with optional original-source feedback allocator; no neural controller or live kinematic force/pose clamp',symmetry_scope='Bilateral source parameters shared; endpoint activations are states,not fitted independent neural weights. Shared role recruitment has not yet been identified.',acceptance=dict(source_ROM='Frozen6-axis bounds at every physics step',velocity_limits_rad_s=rate_limit.tolist(),acceleration_limits_rad_s2=acc_limit.tolist(),swing_peak_gap_m=.01,sustained_gap_m=.001,sustained_unloaded_duration_s=.08,swing_normal_max_BW=.025,penetration_max_m=.002,loaded_contact_tangential_speed_max_mps=.005));mp=out/'MANIFEST.json'
 via=load_recruitment(via_folder) if via_folder is not None else None
 if via is not None:
  via_enabled=via['force_channel_enabled'] if 'force_channel_enabled' in via else np.ones(len(names),bool)
  if not np.array_equal(p,via['parameters']) or names!=via['muscles'].tolist() or not np.array_equal(enabled,via_enabled):raise ValueError('Intermediate keyframe changes source parameters,actuator order,or pruning policy')
  via_report=json.loads((via_folder/'SOURCE_NATIVE_STANCE.json').read_text());via_support=json.loads((Path(via_report['source_pose'])/'MANIFEST.json').read_text())['support_fraction']
  via_config=json.loads((Path(via_report['source_pose'])/'MANIFEST.json').read_text())
  if via_config.get('format','').startswith(('aarl_full_native_static_equilibrium_','aarl_full_native_stationarity_validation_')) and not via_report.get('accepted_stationary_endpoint',False):raise ValueError('New full-native intermediate endpoint requires accepted five-second stationarity before transition training')
  if via_support!=supports[0]:raise ValueError('Intermediate keyframe has different support')
  manifest['unloaded_nearflat_intermediate']=dict(folder=str(via_folder),global_phase=.5,reference='Two equal-duration C2 segments; unload vertically at initial fore-aft position,then lift/advance. Gates unchanged.')
 sequence=[start,via,end] if via is not None else [start,end];sequence_folders=[start_folder,via_folder,end_folder] if via is not None else [start_folder,end_folder]
 if roundtrip:
  sequence+=sequence[-2::-1];sequence_folders+=sequence_folders[-2::-1];manifest['return_to_double_support']=dict(required=True,terminal_min_normal_each_BW=.05,terminal_max_foot_gap_m=.001,terminal_root_translation_error_max_m=.005,terminal_root_rotation_error_max_deg=2.,scope='Continuous stateful physical lift/return; no resetting qpos,qvel,activation,or fiber state between reference segments.')
 durations=np.full(len(sequence)-1,ramp/(len(sequence)-1))
 if return_dwell>0.:
  if not roundtrip or via is None:raise ValueError('Return dwell requires roundtrip and an unloaded intermediate keyframe')
  sequence.insert(len(sequence)-1,via);sequence_folders.insert(len(sequence_folders)-1,via_folder);durations=np.insert(durations,len(durations)-1,return_dwell);manifest['unloaded_return_dwell_s']=return_dwell
 if landing_duration_scale!=1.:
  if not roundtrip or landing_duration_scale<=0:raise ValueError('Landing duration scale requires roundtrip and a positive scale')
  durations[-1]*=landing_duration_scale;manifest['landing_duration_scale']=landing_duration_scale
 reference_boundaries=np.r_[0.,np.cumsum(durations)];motion_duration=float(reference_boundaries[-1]);manifest['reference_sequence']=[str(x) for x in sequence_folders];manifest['segment_durations_s']=durations.tolist()
 manifest['source_feedback_allocator']=dict(enabled=feedback,kp_Nm_rad=180.,kd_Nms_rad=12.,update_stride=control_stride,update_period_s=control_stride*dt,excitation_filter_tau_s=command_filter,feedback_C2_warmup=feedback_warmup,source='Existing hunt_training/motor_teacher.py default shared gains; local original-source static activation tangent,actual dynamic force and contact feedback; no direct torque actuation')
 manifest['source_feedback_allocator']['compensate_coordinate_stop_forces']=compensate_coordinate_stops
 d.qpos[:]=start['qpos'];d.act[:]=a0;d.ctrl[:]=a0;d.xfrc_applied[m.body('pelvis').id,2]=supports[0]*bw;mujoco.mj_forward(m,d);backend=ConstitutiveBackend(m,d,p,mode='compliant')
 if disabled.any():
  class MaskedBank(DynamicThelenBank):
   def evaluate(self,length,activation,dynamic=True):
    super().evaluate(length,activation,dynamic);self.force[disabled]=0.;self.tangent[disabled]=0.;self.error[disabled]=0.;self.bounds[disabled]=0;return self.force
  backend.bank=MaskedBank(p,dt);backend.bank.initialize(d.actuator_length,a0)
 backend.bank.evaluate(d.actuator_length,a0);geometry=[];allocator_log=[]
 stationary_endpoint_source=any(c.get('format','').startswith(('aarl_full_native_static_equilibrium_','aarl_full_native_stationarity_validation_')) for c in configs)
 if (native_acceleration or native_root_augment) and not feedback:raise ValueError('Native acceleration teacher requires source feedback')
 if native_acceleration or native_root_augment or stationary_endpoint_source:
  m.actuator_biasprm[:,0]=-backend.bank.force;mujoco.mj_forward(m,d)
 manifest['initial_derived_force_state']='Validated source static forces refreshed before first control query' if stationary_endpoint_source else 'Historical source initialization preserved'
 if feedback:
  from mechanics import CoupledArms
  from scipy.optimize import lsq_linear
  arms=CoupledArms(m);static_bank=DynamicThelenBank(p,dt);command=a0.copy();filtered=a0.copy();kp=np.full(len(arms.rows),180.);kd=np.full(len(arms.rows),12.)
  if mtp_impedance_scaled:
   arms.update(d);mass=np.zeros((m.nv,m.nv));mujoco.mj_fullM(m,mass,d.qM);reduced=arms.E.T@mass@arms.E;effective_inertia=1./np.diag(np.linalg.inv(reduced));mtp=np.array([k for k,j in enumerate(arms.jids) if m.joint(int(j)).name.startswith('mtp_angle_')]);stiffness=float(np.mean(m.jnt_stiffness[arms.jids[mtp]]));passive_damping=float(np.mean(m.dof_damping[arms.rows[mtp]]));inertia=float(np.mean(effective_inertia[mtp]));kp[mtp]=stiffness;kd[mtp]=max(0.,2.*np.sqrt(2.*stiffness*inertia)-passive_damping);manifest['source_feedback_allocator']['MTP_impedance']=dict(enabled=True,shared_kp_Nm_rad=stiffness,shared_kd_Nms_rad=float(kd[mtp[0]]),native_passive_kp_Nm_rad=stiffness,native_passive_damping_Nms_rad=passive_damping,root_held_equality_reduced_effective_inertia_kgm2=effective_inertia[mtp].tolist(),formula='Common activeKp=mean native passiveKp; activeKd=max(0,2sqrt((activeKp+passiveKp)*mean effective inertia)-mean native damping). Local joint approximation,not whole-contact critical-damping claim. Same gains both sides.')
 elif mtp_impedance_scaled:raise ValueError('MTP impedance requires source feedback allocator')
 acceleration_teacher=NativeAccelerationTeacher(m,arms,enabled) if native_acceleration or native_root_augment else None;last_acceleration_info=None
 source_predictor=None;last_excitation=a0.copy()
 if actuation_audit_times and any(abs(t/(control_stride*dt)-round(t/(control_stride*dt)))>1e-9 for t in actuation_audit_times):raise ValueError('Actuation snapshots must align with the existing control clock')
 if predictive_source_horizon is not None:
  if not enabled.all():raise ValueError('Source forecast currently validated only for the explicit original-source92 variant')
  from predictive_source_recruitment import SourcePredictor
  prediction_steps=round(predictive_source_horizon/dt)
  if prediction_steps<1 or abs(prediction_steps*dt-predictive_source_horizon)>1e-10:raise ValueError('Prediction horizon must be an exact physics-step multiple')
  source_predictor=SourcePredictor(p,m.actuator_dynprm,dt,command_filter);manifest['predictive_source_recruitment']=dict(enabled=True,horizon_s=predictive_source_horizon,steps=prediction_steps,command_derivative_step=.005,active_after_warmup_s=warm,scope='Post-warmup native activation,engineering excitation filter,and source-fiber forecast; causal linear MTU velocity. Original C2 warmup controller preserved continuously. No physical state change.')
 if native_acceleration or native_root_augment:manifest['native_acceleration_teacher']=dict(enabled=True,objective='root_wrench_augments_source_joint_torque' if native_root_augment else 'root_and_joint_mass_weighted_acceleration',omega_rad_s=5.,force_probe_N=.01,source_initialization='Native derived forces refreshed from source static equilibrium before first control query; no physical state change',scope='Source force-to-whole-body acceleration map on separate native model/data; retained contacts/soft followers. Augmentation keeps original joint torque target and converts root response to wrench through reflected root inertia,with source pelvis-height length scaling. Actual-muscle excitation only.')
 if native_acceleration and not mass_consistent_impedance:manifest['source_feedback_allocator']['kp_Nm_rad']=None;manifest['source_feedback_allocator']['kd_Nms_rad']=None
 if mass_consistent_impedance:manifest['mass_consistent_joint_impedance']=dict(enabled=True,root_omega_rad_s=5.,equation='M_JJ delta_a_J + M_JR delta_a_R = Kp error_q + Kd error_v; desired_a=reference_a+delta_a',reference_acceleration_policy='Original endpoint-joint C2 acceleration; not recomputed from online Cartesian IK. Experimental,not a validated inverse-motion fit.',scope='Same shared PD gains; native actual acceleration allocation. No static contact-torque cancellation or injected torques.')
 cartesian_teacher=None
 if cartesian_root:
  if not feedback:raise ValueError('Root-aware Cartesian targets require source feedback')
  from supported_step_teacher import NativeCartesianTeacher
  cartesian_teacher=NativeCartesianTeacher(m,sequence,spec,stance_posture_cartesian,preserve_loaded_followers=loaded_cartesian_followers,tangential_clearance=cartesian_tangential_clearance);manifest['root_aware_cartesian_teacher']=dict(enabled=True,stance_posture=stance_posture_cartesian,preserve_loaded_followers=loaded_cartesian_followers,tangential_clearance_margin_m=cartesian_tangential_clearance,target_position_tolerance_m=.0002,target_rotation_tolerance_rad=.001,knee_target_extension_limit_deg=-2.,stance_contact_blend='C2 based on declared reference whole-foot gap,fully stance at<=0mm and fully swing at>=1mm',scope='Desired targets only; optional clearance margin settles horizontal/orientation motion before vertical touchdown. Full loaded FK/velocity retained. No live state/external force changes; all physical gates unchanged.')
 pitch_teacher=None
 if world_foot_pitch:
  from supported_step_teacher import NativeFootPitchTeacher
  pitch_teacher=NativeFootPitchTeacher(m,sequence,spec);manifest['world_foot_pitch_teacher']=dict(enabled=True,source='Existing native contact-calibrated keyframes; no added dorsiflexion',scope='Bounded ankle desired targets from actual world-foot pitch/rate; original hip/knee targets unchanged. No physical geometry/strength changes.')
 if export_training_reference:manifest['training_reference_export']=dict(file='training_reference.npz',sample_period_s=dt,scope='Actual continuous mechanical state and applied excitation; not a fitted neural network. Activation_for_force is the split-step constitutive activation,activation is the end-step MuJoCo state. Fiber velocity is backward-difference along-fiber m/s,consistent with source constitutive equations.')
 manifest['loaded_contact_slip_threshold_BW']=loaded_contact_min_BW
 if mp.exists() and json.loads(mp.read_text())!=manifest:raise ValueError('Existing transition manifest differs')
 mp.write_text(json.dumps(manifest,indent=2))
 for g in range(m.ngeom):
  name=m.body(int(m.geom_bodyid[g])).name
  if m.geom_type[g]==mujoco.mjtGeom.mjGEOM_MESH and name in ('calcn_r','toes_r','calcn_l','toes_l'):
   mesh=int(m.geom_dataid[g]);first=int(m.mesh_vertadr[mesh]);count=int(m.mesh_vertnum[mesh]);geometry.append((0 if name.endswith('_r') else 1,g,m.mesh_vert[first:first+count].copy()))
 n=round((warm+motion_duration+hold)/dt);keep=np.arange(0,n,max(1,stride));telemetry=np.empty((len(keep),1+m.nq+m.nu+2+2+6));gaps=np.full((n,2),np.inf);loads=np.zeros((n,2));q6=np.empty((n,len(addresses)));v6=np.empty_like(q6);acc6=np.empty_like(q6);root=np.empty((n,6));force6=np.zeros(6);jp=np.zeros((3,m.nv));jr=np.zeros_like(jp);M=np.zeros((m.nv,m.nv));maximum_residual=0.;maximum_slip=0.;slip_event={};maximum_acceleration_ratio=0.;acceleration_event={};guard='';written=0;started=perf_counter()
 training={k:np.empty((n,width)) for k,width in dict(qpos=m.nq,qvel=m.nv,qacc=m.nv,activation=m.nu,excitation=m.nu,activation_for_force=m.nu,mtu_length_m=m.nu,source_force_N=m.nu,fiber_length_m=m.nu,fiber_velocity_mps=m.nu,root_passive_generalized_force=6).items()} if export_training_reference else None
 snapshot_indices={round(t/dt) for t in (actuation_audit_times or [])};snapshots=[];state_mask=mujoco.mjtState.mjSTATE_INTEGRATION
 for i in range(n):
  if i in snapshot_indices:
   state=np.empty(mujoco.mj_stateSize(m,state_mask));mujoco.mj_getState(m,d,state,state_mask);snapshots.append(dict(step=i,time=float(d.time),state=state,bias=m.actuator_biasprm.copy(),fiber=backend.bank.state.copy(),force=backend.bank.force.copy(),excitation=d.ctrl.copy(),command=command.copy(),qacc=d.qacc.copy(),window=[]))
  if return_dwell>0. or landing_duration_scale!=1.:
   elapsed=np.clip(d.time-warm,0.,motion_duration);index=min(np.searchsorted(reference_boundaries,elapsed,side='right')-1,len(durations)-1);phase=np.clip((elapsed-reference_boundaries[index])/durations[index],0.,1.);segment_start,segment_end=sequence[index:index+2];segment_duration=durations[index]
  else:phase=np.clip((d.time-warm)/ramp,0.,1.);segments=len(sequence)-1;index=min(int(phase*segments),segments-1);segment_start,segment_end=sequence[index:index+2];phase=np.clip(phase*segments-index,0.,1.);segment_duration=ramp/segments
  blend=phase**3*(10.+phase*(-15.+6.*phase));reference_alpha=segment_start['activation']+blend*(segment_end['activation']-segment_start['activation']);d.ctrl[:]=reference_alpha
  if feedback and (d.time>=warm or feedback_warmup):
   if i%control_stride==0:
    rho=arms.update(d).copy();L=d.actuator_length.copy();lo=np.maximum(.01,d.act-.005);hi=np.minimum(.98,d.act+.005);static_bank.initialize(L,lo);fl=static_bank.force.copy();static_bank.initialize(L,hi);gain=np.maximum(0.,(static_bank.force-fl)/(hi-lo));force_intercept=backend.bank.force.copy();control_intercept=d.act.copy();prediction=None
    if source_predictor is not None and d.time>=warm:
     forecast=source_predictor.linearize_command(backend.bank.state,d.act,last_excitation,command,L,d.actuator_velocity,prediction_steps)
     if forecast is None:guard='Source recruitment forecast invalid'
     else:force_intercept,gain,control_intercept=forecast;prediction=(force_intercept,control_intercept,backend.bank.force.copy())
    gain[disabled]=0.;force_intercept[disabled]=0.;A=rho*gain[None,:]
    desired_q=segment_start['qpos']+blend*(segment_end['qpos']-segment_start['qpos']);rate=(30*phase**2*(1-phase)**2/segment_duration) if 0.<phase<1. else 0.;desired_v=rate*(segment_end['qpos']-segment_start['qpos']);root_reference_q=desired_q[:6].copy();root_reference_v=desired_v[:6].copy();acc_rate=60.*phase*(1.-phase)*(1.-2.*phase)/segment_duration**2 if 0.<phase<1. else 0.;desired_a=acc_rate*(segment_end['qpos']-segment_start['qpos']);constraint_for_control=d.qfrc_constraint if compensate_coordinate_stops else d.qfrc_constraint-coordinate_stop_force(m,d)
    if cartesian_teacher is not None:
     desired_q,desired_v=cartesian_teacher.solve(d,index,blend,rate,desired_q,desired_v)
     if not cartesian_teacher.diagnostics()['accepted_target_IK']:guard='Root-aware bounded target IK did not meet pose criterion'
    if pitch_teacher is not None:desired_q,desired_v=pitch_teacher.solve(d,index,blend,rate,desired_q,desired_v)
    if acceleration_teacher is not None:
     desired_q[:6]=root_reference_q;desired_v[:6]=root_reference_v
     if native_root_augment:
      tau=arms.E.T@(d.qfrc_bias-d.qfrc_passive-constraint_for_control)+kp*(desired_q[arms.addresses]-d.qpos[arms.addresses])+kd*(desired_v[arms.rows]-d.qvel[arms.rows]);rhs=tau-rho@force_intercept+A@control_intercept;command,last_acceleration_info=acceleration_teacher.update(d,gain,d.act,desired_q,desired_v,desired_a,reference_alpha,A,rhs,source_prediction=prediction);last_acceleration_info['linearized_joint_torque_residual_Nm']=(A@command-rhs).tolist()
     else:command,last_acceleration_info=acceleration_teacher.update(d,gain,d.act,desired_q,desired_v,desired_a,reference_alpha,joint_impedance=(kp,kd) if mass_consistent_impedance else None,source_prediction=prediction)
     residual_norm=last_acceleration_info['residual_norm']
     if not last_acceleration_info['valid_map']:guard='Native acceleration map validation failed'
    else:
     tau=arms.E.T@(d.qfrc_bias-d.qfrc_passive-constraint_for_control)+kp*(desired_q[arms.addresses]-d.qpos[arms.addresses])+kd*(desired_v[arms.rows]-d.qvel[arms.rows]);rhs=tau-rho@force_intercept+A@control_intercept;norm=np.linalg.norm(A,axis=0);regularization=.02*np.median(norm[norm>1e-5]);fit=lsq_linear(np.vstack([A,regularization*np.eye(m.nu)]),np.r_[rhs,regularization*reference_alpha],bounds=(.01,.98),method='bvls',tol=1e-6,max_iter=100);command=fit.x;residual_norm=float(np.linalg.norm(A@command-rhs))
    allocator_log.append((d.time,residual_norm,float(np.sqrt(np.mean((desired_q[arms.addresses]-d.qpos[arms.addresses])**2)))))
    for snapshot in snapshots:
     if snapshot['step']==i:snapshot.update(reference_q=desired_q.copy(),reference_v=desired_v.copy(),reference_a=desired_a.copy(),kp=kp.copy(),kd=kd.copy())
   warm_phase=np.clip(d.time/warm,0.,1.);warm_blend=warm_phase**3*(10.+warm_phase*(-15.+6.*warm_phase)) if feedback_warmup else 1.;target_command=reference_alpha+warm_blend*(command-reference_alpha)
   if command_filter>0.:filtered+=(1.-np.exp(-dt/command_filter))*(target_command-filtered);d.ctrl[:]=filtered
   else:d.ctrl[:]=target_command
  if source_predictor is not None:last_excitation[:]=d.ctrl
  try:backend.step()
  except RuntimeError as error:guard=str(error)
  for snapshot in snapshots:
   if 0<=i-snapshot['step']<control_stride:snapshot['window'].append(dict(excitation=d.ctrl.copy(),qpos=d.qpos.copy(),qvel=d.qvel.copy(),activation=d.act.copy(),force=backend.bank.force.copy(),qacc=d.qacc.copy()))
  q6[i]=d.qpos[addresses];v6[i]=d.qvel[addresses];acc6[i]=d.qacc[addresses];root[i]=d.qpos[:6]
  ratio=float(np.max(np.divide(abs(acc6[i]),acc_limit,out=np.zeros_like(acc_limit),where=acc_limit>0)))
  if d.time>=warm and ratio>maximum_acceleration_ratio:
   maximum_acceleration_ratio=ratio;stop=coordinate_stop_force(m,d);acceleration_event=dict(time_s=float(d.time),source_acceleration_ratio=ratio,joint_names=joint_names,angles_deg=np.rad2deg(q6[i]).tolist(),acceleration_rad_s2=acc6[i].tolist(),coordinate_stop_torque_Nm=[float(stop[m.joint(n).dofadr[0]]) for n in joint_names],active_joint_stops=[dict(joint=m.joint(int(d.efc_id[k])).name,multiplier=float(d.efc_force[k]),position=float(d.efc_pos[k]),margin=float(d.efc_margin[k])) for k in range(d.nefc) if d.efc_type[k]==mujoco.mjtConstraint.mjCNSTR_LIMIT_JOINT and abs(d.efc_force[k])>1e-8])
  if training is not None:
   length=np.hypot(backend.bank.state,p[:,3]);velocity=(backend.bank.state-backend.old_state)*np.divide(backend.bank.state,length)/(dt)
   values=dict(qpos=d.qpos,qvel=d.qvel,qacc=d.qacc,activation=d.act,excitation=d.ctrl,activation_for_force=backend.activation,mtu_length_m=d.actuator_length,source_force_N=backend.bank.force,fiber_length_m=length,fiber_velocity_mps=velocity,root_passive_generalized_force=d.qfrc_passive[:6])
   for key,value in values.items():training[key][i]=value
  for side,g,vertices in geometry:gaps[i,side]=min(gaps[i,side],float(np.min(vertices@d.geom_xmat[g].reshape(3,3)[2]+d.geom_xpos[g,2])))
  for ci in range(d.ncon):
   c=d.contact[ci];g1,g2=int(c.geom1),int(c.geom2)
   if m.geom_type[g1]==mujoco.mjtGeom.mjGEOM_PLANE:g=g2
   elif m.geom_type[g2]==mujoco.mjtGeom.mjGEOM_PLANE:g=g1
   else:continue
   body=int(m.geom_bodyid[g]);name=m.body(body).name
   if not name.endswith(('_r','_l')):continue
   mujoco.mj_contactForce(m,d,ci,force6);side=0 if name.endswith('_r') else 1;loads[i,side]+=force6[0]/bw
   if d.time>=warm and force6[0]/bw>loaded_contact_min_BW:
    mujoco.mj_jac(m,d,jp,jr,c.pos,body);slip=float(np.linalg.norm((jp@d.qvel)[:2]))
    if slip>maximum_slip:
     maximum_slip=slip;slip_event=dict(time_s=float(d.time),body=name,normal_force_BW=float(force6[0]/bw),speed_mps=slip,contact_world_m=c.pos.tolist())
     if acceleration_teacher is not None:slip_event['native_acceleration_diagnostic']=last_acceleration_info
     elif feedback:
      slip_event['allocator_diagnostic']=dict(last_update_time_s=float(allocator_log[-1][0]),joint_names=[m.joint(int(j)).name for j in arms.jids],linearized_allocation_residual_Nm=(A@command-rhs).tolist(),actual_source_torque_minus_requested_Nm=(rho@backend.bank.force-tau).tolist(),joint_error_rad=(desired_q[arms.addresses]-d.qpos[arms.addresses]).tolist(),joint_velocity_rad_s=d.qvel[arms.rows].tolist(),saturated_commands=[names[k] for k in np.flatnonzero(command>.979)],toe_muscles=[dict(name=n,activation=float(d.act[k]),command=float(command[k]),force_N=float(backend.bank.force[k])) for k,n in enumerate(names) if n.startswith(('flex_dig','flex_hal','ext_dig','ext_hal'))],scope='Last allocator update and live source force; Jacobian/request may precede this contact by up to10ms. Diagnostic only,no control changes.')
  mujoco.mj_fullM(m,M,d.qM);maximum_residual=max(maximum_residual,float(np.max(abs(M@d.qacc-d.qfrc_smooth-d.qfrc_constraint))))
  if i%max(1,stride)==0:telemetry[written]=np.r_[d.time,d.qpos,d.act,loads[i],gaps[i],d.qfrc_passive[:6]];written+=1
  if np.any(d.qpos[[m.joint('knee_angle_'+s).qposadr[0] for s in ('r','l')]]>0.):guard=guard or 'Strict knee extension'
  if np.any(q6[i]<low-1e-10) or np.any(q6[i]>high+1e-10):guard=guard or 'Frozen source ROM'
  if guard:break
 count=i+1;times=dt*np.arange(1,count+1);gaps=gaps[:count];loads=loads[:count];q6=q6[:count];v6=v6[:count];acc6=acc6[:count];root=root[:count];window=times>=warm;terminal=times>=warm+motion_duration
 # Infer airborne side from the saved endpoint geometry,not a naming/phase guess.
 probe=mujoco.MjData(m);probe.qpos[:]=end['qpos'];mujoco.mj_fwdPosition(m,probe);endgap=np.full(2,np.inf)
 for side,g,vertices in geometry:endgap[side]=min(endgap[side],float(np.min(vertices@probe.geom_xmat[g].reshape(3,3)[2]+probe.geom_xpos[g,2])))
 swing=int(np.argmax(endgap));qualified=(gaps[:,swing]>.001)&(loads[:,swing]<.025)&window;longest=current=0
 for value in qualified:current=current+1 if value else 0;longest=max(longest,current)
 rate=np.max(abs(v6[window]),axis=0) if window.any() else np.zeros(len(addresses));acc=np.max(abs(acc6[window]),axis=0) if window.any() else np.zeros(len(addresses));metrics=dict(duration_s=float(d.time),guard=guard,swing_side=('r','l')[swing],maximum_swing_clearance_m=float(gaps[window,swing].max()) if window.any() else None,sustained_unloaded_clearance_s=float(longest*dt),terminal_min_swing_clearance_m=float(gaps[terminal,swing].min()) if terminal.any() else None,terminal_max_swing_normal_BW=float(loads[terminal,swing].max()) if terminal.any() else None,minimum_foot_mesh_gap_m=gaps.min(axis=0).tolist(),maximum_loaded_contact_tangential_speed_mps=maximum_slip,maximum_slip_event=slip_event,maximum_generalized_dynamics_residual=maximum_residual,maximum_joint_velocity_rad_s=rate.tolist(),maximum_joint_acceleration_rad_s2=acc.tolist(),velocity_reference_ratio=np.divide(rate,rate_limit,out=np.zeros_like(rate),where=rate_limit>0).tolist(),acceleration_reference_ratio=np.divide(acc,acc_limit,out=np.zeros_like(acc),where=acc_limit>0).tolist(),maximum_knee_extension_deg=float(np.rad2deg(q6[:,[k for k,n in enumerate(joint_names) if n.startswith('knee_angle_')]].max())),wall_s=perf_counter()-started);accepted=not guard and longest*dt>=.08 and metrics['maximum_swing_clearance_m']>=.01 and gaps.min()>=-.002 and maximum_slip<=.005 and np.all(rate<=rate_limit+1e-10) and np.all(acc<=acc_limit+1e-10)
 metrics['maximum_acceleration_event']=acceleration_event
 if roundtrip:
  recovery=dict(terminal_min_normal_each_BW=loads[terminal].min(axis=0).tolist() if terminal.any() else None,terminal_max_foot_gap_m=gaps[terminal].max(axis=0).tolist() if terminal.any() else None,terminal_root_translation_error_max_m=float(np.max(abs(root[terminal,:3]-start['qpos'][:3]))) if terminal.any() else None,terminal_root_rotation_error_max_deg=float(np.rad2deg(np.max(abs(root[terminal,3:]-start['qpos'][3:6])))) if terminal.any() else None);recovery['accepted']=bool(terminal.any() and np.min(recovery['terminal_min_normal_each_BW'])>=.05 and max(recovery['terminal_max_foot_gap_m'])<=.001 and recovery['terminal_root_translation_error_max_m']<=.005 and recovery['terminal_root_rotation_error_max_deg']<=2.);metrics['double_support_recovery']=recovery;accepted=accepted and recovery['accepted']
 if cartesian_teacher is not None:metrics['root_aware_cartesian_teacher']=cartesian_teacher.diagnostics();accepted=accepted and metrics['root_aware_cartesian_teacher']['accepted_target_IK']
 if pitch_teacher is not None:metrics['world_foot_pitch_teacher']=pitch_teacher.diagnostics()
 if acceleration_teacher is not None:metrics['native_acceleration_teacher']=acceleration_teacher.diagnostics();accepted=accepted and acceleration_teacher.invalid_maps==0
 if source_predictor is not None:metrics['predictive_source_recruitment']=source_predictor.diagnostics();accepted=accepted and source_predictor.invalid==0
 metrics['constitutive_backend']=backend.diagnostics()
 if snapshots:
  complete=[s for s in snapshots if len(s['window'])==control_stride];metrics['actuation_audit']=dict(requested_times_s=actuation_audit_times,complete_windows=len(complete),window_steps=control_stride,file='actuation_windows.npz')
  if complete:
   payload=dict(times=np.array([s['time'] for s in complete]),state_spec=int(state_mask),integration_state=np.stack([s['state'] for s in complete]),actuator_biasprm=np.stack([s['bias'] for s in complete]),fiber_state=np.stack([s['fiber'] for s in complete]),source_force=np.stack([s['force'] for s in complete]),excitation=np.stack([s['excitation'] for s in complete]),command=np.stack([s['command'] for s in complete]),initial_qacc=np.stack([s['qacc'] for s in complete]),parameters=p,muscles=names,force_channel_enabled=enabled,timestep_s=dt,filter_tau_s=command_filter,window_steps=control_stride)
   for key in ('excitation','qpos','qvel','activation','force','qacc'):payload['window_'+key]=np.array([[step[key] for step in s['window']] for s in complete])
   for key in ('reference_q','reference_v','reference_a','kp','kd'):payload[key]=np.array([s[key] for s in complete])
   payload['compensate_coordinate_stops']=compensate_coordinate_stops
   np.savez_compressed(out/'actuation_windows.npz',**payload)
 if training is not None:
  np.savez_compressed(out/'training_reference.npz',t=times,muscles=names,parameters=p,force_channel_enabled=enabled,foot_mesh_gaps_m=gaps,normal_side_BW=loads,**{k:v[:count] for k,v in training.items()});metrics['training_reference_compressed_bytes']=(out/'training_reference.npz').stat().st_size
 report=dict(accepted_dynamic_teacher=bool(accepted),metrics=metrics,manifest=manifest);(out/'TRANSITION_RESULT.json').write_text(json.dumps(report,indent=2));np.savez_compressed(out/'transition_sparse.npz',samples=telemetry[:written],nq=m.nq,nu=m.nu,joint_names=joint_names,physics_t=times,angles_sagittal=q6,foot_mesh_gaps_m=gaps,normal_side_BW=loads,allocator_samples=allocator_log);print(json.dumps(dict(accepted_dynamic_teacher=bool(accepted),metrics=metrics),indent=2))

def native_stance_recruitment(model_folder,folder,out,hold_dt=None,hold_mode='compliant',preserve_pruned=False,activation_floor=.01):
 """Original-source recruitment against native inverse forces,then unclamped hold."""
 from mechanics import CoupledArms,fill_projection
 from contact_standing_solver import _external
 from scipy.optimize import lsq_linear
 from thelen_dynamic import ConstitutiveBackend
 if not json.loads((folder/'NATIVE_STANCE_RESULT.json').read_text())['accepted_native_static_pose']:raise ValueError('Native contact pose did not pass')
 out.mkdir(parents=True,exist_ok=True);m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);q=np.load(folder/'native_stance.npz')['qpos'];d.qpos[:]=q;d.qvel[:]=0.;d.qacc[:]=0.;bw=float(m.body_mass.sum()*9.81);support=json.loads((folder/'MANIFEST.json').read_text())['support_fraction'];d.xfrc_applied[m.body('pelvis').id,2]=support*bw;mujoco.mj_inverse(m,d);ar=CoupledArms(m);drivers=np.array([i for i in range(m.nv) if i not in set(ar.dependent.tolist())],dtype=np.int64);E=np.zeros((m.nv,len(drivers)));fill_projection(E,drivers,ar.dependent,ar.independent,ar.coefficients,ar.offsets,d.qpos,ar.qadr);rho=E.T@(-d.actuator_moment.T);full_required=d.qfrc_inverse-_external(m,d);target=E.T@full_required;names=[m.actuator(i).name for i in range(m.nu)];p=parameters_from_osim(OSIM,names,0.);bank=DynamicThelenBank(p,.002);alpha=np.full(m.nu,.05);scale=np.array([1./bw if m.jnt_type[int(m.dof_jntid[i])]==mujoco.mjtJoint.mjJNT_SLIDE else 1./(bw*.5) for i in drivers]);history=[]
 disabled=np.zeros(m.nu,dtype=bool)
 if preserve_pruned:
  from runner import PRUNE_MUSCLES
  disabled=np.array([name in PRUNE_MUSCLES for name in names])
 class MaskedBank(DynamicThelenBank):
  def evaluate(self,length,activation,dynamic=True):
   super().evaluate(length,activation,dynamic);self.force[disabled]=0.;self.tangent[disabled]=0.;self.error[disabled]=0.;self.bounds[disabled]=0;return self.force
 if preserve_pruned:bank=MaskedBank(p,.002)
 def force(a):bank.initialize(d.actuator_length,a);return bank.force.copy()
 def cost(a,f):return float(np.sum((scale*(rho@f-target))**2)+1e-8*np.sum(a*a))
 f=force(alpha)
 for iteration in range(20):
  oldcost=cost(alpha,f);low=np.maximum(activation_floor,alpha-.005);high=np.minimum(.98,alpha+.005);derivative=np.maximum(0.,(force(high)-force(low))/(high-low));A=rho*derivative[None,:];rhs=target-rho@f+A@alpha;fit=lsq_linear(np.vstack([scale[:,None]*A,np.eye(m.nu)*.0001]),np.r_[scale*rhs,np.zeros(m.nu)],bounds=(activation_floor,.98),method='bvls',tol=1e-9,max_iter=200);best=(oldcost,alpha.copy(),f.copy())
  for step in (1.,.5,.25,.125,.0625):
   candidate=alpha+step*(fit.x-alpha);ff=force(candidate);value=cost(candidate,ff)
   if value<best[0]:best=(value,candidate.copy(),ff.copy())
  alpha,f=best[1:];history.append(dict(iteration=iteration,cost=best[0],maximum_anatomical_error=float(np.max(abs(rho@f-target))),max_activation=float(alpha.max())))
  if oldcost-best[0]<1e-13:break
 f=force(alpha);residual=rho@f-target;domain=not bool(np.any(bank.bounds));accepted=domain and np.max(abs(residual))<=1.;np.savez_compressed(out/'source_recruitment.npz',qpos=q,activation=alpha,muscles=names,source_force=f,parameters=p,force_channel_enabled=~disabled);report=dict(accepted_static_source_recruitment=bool(accepted),source_domain_valid=domain,controller_activation_floor=activation_floor,maximum_anatomical_residual_N_or_Nm=float(np.max(abs(residual))),raw_generalized_residual_max=float(np.max(abs(-d.actuator_moment.T@f-full_required))),source_model=str(model_folder/'runtime_model.mjb'),source_pose=str(folder),actuation_policy='Legacy production short-rotator pruning preserved' if preserve_pruned else 'Original source92,including6shortrotators disabled by legacy runner',zero_force_muscles=[name for i,name in enumerate(names) if disabled[i]],history=history,scope='Native inverse contact law and original source Thelen static equilibrium,with declared actuation policy. Raw follower residual is reported; soft equality settling is tested by forward hold. No neural controller or trajectory clamp.')
 if accepted:
  m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
  if hold_dt is not None:m.opt.timestep=hold_dt
  d=mujoco.MjData(m);d.qpos[:]=q;d.act[:]=alpha;d.ctrl[:]=alpha;d.xfrc_applied[m.body('pelvis').id,2]=support*bw;mujoco.mj_forward(m,d);backend=ConstitutiveBackend(m,d,p,mode=hold_mode)
  if preserve_pruned:
   if hold_mode!='compliant':raise ValueError('Explicit masked rigid adapter not validated')
   backend.bank=MaskedBank(p,m.opt.timestep);backend.bank.initialize(d.actuator_length,alpha)
  backend.bank.evaluate(d.actuator_length,alpha);n=round(2./m.opt.timestep);qq=np.empty((n,m.nq));normal=np.empty(n);normal_side=np.zeros((n,2));ff=np.empty((n,m.nu));meshgap=np.full((n,2),np.inf);guard='';kneeaddr=[int(m.joint('knee_angle_'+s).qposadr[0]) for s in ('r','l')];force6=np.zeros(6);footgeoms=[]
  for g in range(m.ngeom):
   name=m.body(int(m.geom_bodyid[g])).name
   if m.geom_type[g]==mujoco.mjtGeom.mjGEOM_MESH and name in ('calcn_r','toes_r','calcn_l','toes_l'):
    mesh=int(m.geom_dataid[g]);first=int(m.mesh_vertadr[mesh]);count=int(m.mesh_vertnum[mesh]);footgeoms.append((0 if name.endswith('_r') else 1,g,m.mesh_vert[first:first+count].copy()))
  initialgap=np.full(2,np.inf)
  for side,g,vertices in footgeoms:initialgap[side]=min(initialgap[side],float(np.min(vertices@d.geom_xmat[g].reshape(3,3)[2]+d.geom_xpos[g,2])))
  started=perf_counter()
  for i in range(n):
   backend.step();qq[i]=d.qpos;ff[i]=backend.bank.force;total=0.
   for ci in range(d.ncon):
    con=d.contact[ci]
    if m.geom_type[int(con.geom1)]==mujoco.mjtGeom.mjGEOM_PLANE or m.geom_type[int(con.geom2)]==mujoco.mjtGeom.mjGEOM_PLANE:
     mujoco.mj_contactForce(m,d,ci,force6);total+=force6[0];g=int(con.geom2) if m.geom_type[int(con.geom1)]==mujoco.mjtGeom.mjGEOM_PLANE else int(con.geom1);name=m.body(int(m.geom_bodyid[g])).name
     if name.endswith(('_r','_l')):normal_side[i,0 if name.endswith('_r') else 1]+=force6[0]/bw
   normal[i]=total/bw
   for side,g,vertices in footgeoms:meshgap[i,side]=min(meshgap[i,side],float(np.min(vertices@d.geom_xmat[g].reshape(3,3)[2]+d.geom_xpos[g,2])))
   if np.any(d.qpos[kneeaddr]>0.):guard='Strict knee extension';break
   if np.any(backend.bank.bounds):guard='Source fiber bound';break
   if not np.all(np.isfinite(d.qpos)):guard='Nonfinite';break
  elapsed=perf_counter()-started;qq=qq[:i+1];normal=normal[:i+1];normal_side=normal_side[:i+1];ff=ff[:i+1];meshgap=meshgap[:i+1];hold=dict(duration_s=float(d.time),physics=hold_mode,timestep_s=float(m.opt.timestep),wall_s=elapsed,physical_seconds_per_wall_second=float(d.time/elapsed),guard=guard,root_translation_drift_max_m=float(np.max(abs(qq[:,:3]-q[:3]))),root_rotation_drift_max_deg=float(np.rad2deg(np.max(abs(qq[:,3:6]-q[3:6])))),max_knee_extension_deg=float(np.rad2deg(qq[:,kneeaddr].max())),normal_load_mean_BW=float(np.mean(normal)),normal_load_range_BW=[float(np.min(normal)),float(np.max(normal))],initial_foot_mesh_gaps_m=initialgap.tolist(),minimum_foot_mesh_gaps_m=meshgap.min(axis=0).tolist(),maximum_foot_normal_load_BW=normal_side.max(axis=0).tolist(),scope='Two-second matched constant source activation hold,original model/guide/support,no neural controller and no pose/velocity clamps. Wall timing excludes model/recruitment setup and warms the constitutive kernel; includes contact-force and mesh-gap logging. Rigid uses the same compliant-calibrated activation,not a separately refitted controller.');hold['accepted']=not guard and hold['root_translation_drift_max_m']<=.005 and hold['root_rotation_drift_max_deg']<=2.;hold['accepted_contact_envelope']=bool(np.min(meshgap)>=-.002 and all(float(meshgap[:,s].min())>=.01 and float(normal_side[:,s].max())<=.025 for s in range(2) if initialgap[s]>.005));report['forward_hold']=hold;np.savez_compressed(out/'source_hold.npz',qpos=qq,normal_load_BW=normal,normal_side_BW=normal_side,muscle_force=ff,foot_mesh_gaps_m=meshgap)
 (out/'SOURCE_NATIVE_STANCE.json').write_text(json.dumps(report,indent=2));print(json.dumps({key:value for key,value in report.items() if key!='history'},indent=2))

def teacher_inverse(model_folder,folder,out):
 """Necessary full-root static source-muscle inverse gate for a new teacher."""
 from mechanics import CoupledArms,fill_projection
 from contact_standing_solver import _external
 from scipy.optimize import lsq_linear
 result=json.loads((folder/'KINEMATIC_RESULT.json').read_text())
 if not result['accepted_kinematics']:raise ValueError('Teacher kinematics have not passed; do not fit recruitment')
 out.mkdir(parents=True,exist_ok=True);m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);ar=CoupledArms(m);z=np.load(folder/'teacher.npz');cfg=json.loads((folder/'MANIFEST.json').read_text());names=[m.actuator(i).name for i in range(m.nu)];p=parameters_from_osim(OSIM,names,0.);bank=DynamicThelenBank(p,.002);bw=float(m.body_mass.sum()*9.81);drivers=np.array([i for i in range(m.nv) if i not in set(ar.dependent.tolist())],dtype=np.int64);E=np.zeros((m.nv,len(drivers)));M=np.zeros((m.nv,m.nv));jp=np.zeros((3,m.nv));jr=np.zeros_like(jp);root=np.flatnonzero(drivers<6);internal=np.flatnonzero(drivers>=6);scale=np.array([1./bw if m.jnt_type[int(m.dof_jntid[i])]==mujoco.mjtJoint.mjJNT_SLIDE else 1./(bw*.5) for i in drivers]);rows=[];alphas=[];prior=np.full(m.nu,.05);indices=np.linspace(0,len(z['t'])-1,33).astype(int)
 manifest=dict(kind='Full-root necessary source-static recruitment inverse with reference dynamics',source_teacher=str(folder),source_teacher_sha256=__import__('hashlib').sha256((folder/'teacher.npz').read_bytes()).hexdigest(),support_fraction=cfg['support_fraction'],root_linear_error_max_BW=.01,root_angular_error_max_BWm=.01,internal_rms_error_max_Nm=1.,internal_peak_error_max_Nm=5.,scope='Original source Thelen static fiber equilibria; actual model inertia,guide spring/damping,original source dimensions,and contact geometry. Optimized contact reactions are feasibility variables,not applied forces or forward validation. Dynamic fiber/activation-history checks remain required.')
 (out/'INVERSE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
 for index in indices:
  d.qpos[:]=z['qpos'][index];ar.project(d);d.qvel[:]=z['qvel'][index];d.xfrc_applied[:]=0.;d.xfrc_applied[m.body('pelvis').id,2]=cfg['support_fraction']*bw;mujoco.mj_forward(m,d);fill_projection(E,drivers,ar.dependent,ar.independent,ar.coefficients,ar.offsets,d.qpos,ar.qadr);mujoco.mj_fullM(m,M,d.qM);rho=E.T@(-d.actuator_moment.T);tau=E.T@(M@z['qacc'][index]+d.qfrc_bias-d.qfrc_passive-_external(m,d));rays=[];contacts=[]
  for ci in range(d.ncon):
   con=d.contact[ci];g1,g2=int(con.geom1),int(con.geom2)
   if m.geom_type[g1]==mujoco.mjtGeom.mjGEOM_PLANE:geom=g2
   elif m.geom_type[g2]==mujoco.mjtGeom.mjGEOM_PLANE:geom=g1
   else:continue
   body=int(m.geom_bodyid[geom]);owner=body
   while owner and m.body(owner).name not in tuple(n+'_'+s for s in ('r','l') for n in ('talus','calcn','toes')):owner=int(m.body_parentid[owner])
   if not owner:continue
   mujoco.mj_jac(m,d,jp,jr,con.pos,body);velocity=jp@d.qvel;mu=float(con.friction[0]);normal=np.array([0.,0.,1.]);basis=[normal+mu*np.array(v) for v in ([1.,0.,0.],[-1.,0.,0.],[0.,1.,0.],[0.,-1.,0.])]
   if np.linalg.norm(velocity[:2])>.005:
    axis=int(np.argmax(abs(velocity[:2])));v=np.zeros(3);v[axis]=-mu*np.sign(velocity[axis]);basis=[normal+v]
   rays.extend([E.T@jp.T@v*bw for v in basis]);contacts.append(dict(body=m.body(owner).name,velocity_mps=velocity.tolist(),gap_m=float(con.dist)))
  d.qacc[:]=z['qacc'][index];mujoco.mj_inverse(m,d);native_required=E.T@(d.qfrc_inverse-_external(m,d));linear=[k for k in root if m.jnt_type[int(m.dof_jntid[drivers[k]])]==mujoco.mjtJoint.mjJNT_SLIDE];angular=[k for k in root if k not in linear];native_contact_root=dict(linear_max_BW=float(np.max(abs(native_required[linear]))/bw),angular_max_BWm=float(np.max(abs(native_required[angular]))/bw),unmet_root_generalized_force=native_required[root].tolist(),scope='Native MuJoCo inverse constraint/contact forces at prescribed q,v,a; no optimized virtual reaction. Internal muscle fit still required.');native_contact_root['accepted']=native_contact_root['linear_max_BW']<=.01 and native_contact_root['angular_max_BWm']<=.01
  C=np.column_stack(rays) if rays else np.zeros((len(drivers),0));alpha=prior.copy();reaction=np.zeros(C.shape[1]);success=True
  def force(a):bank.initialize(d.actuator_length,a);return bank.force.copy()
  f=force(alpha)
  for _ in range(8):
   lo=np.maximum(.01,alpha-.005);hi=np.minimum(.98,alpha+.005);derivative=np.maximum(0.,(force(hi)-force(lo))/(hi-lo));A=rho*derivative[None,:];design=np.column_stack([A,C]);rhs=tau-rho@f+A@alpha;reg=np.diag(np.r_[np.full(m.nu,.0005),np.full(C.shape[1],.00001)]);fit=lsq_linear(np.vstack([scale[:,None]*design,reg]),np.r_[scale*rhs,np.zeros(design.shape[1])],bounds=(np.r_[np.full(m.nu,.01),np.zeros(C.shape[1])],np.r_[np.full(m.nu,.98),np.full(C.shape[1],1.5)]),method='bvls',tol=1e-8,max_iter=200);success=success and fit.success;alpha=fit.x[:m.nu];reaction=fit.x[m.nu:];f=force(alpha)
  residual=rho@f+C@reaction-tau;error=dict(root_linear_max_BW=float(np.max(abs(residual[linear]))/bw),root_angular_max_BWm=float(np.max(abs(residual[angular]))/bw),internal_rms_Nm=float(np.sqrt(np.mean(residual[internal]**2))),internal_max_Nm=float(np.max(abs(residual[internal]))));accepted=success and not np.any(bank.bounds) and error['root_linear_max_BW']<=.01 and error['root_angular_max_BWm']<=.01 and error['internal_rms_Nm']<=1. and error['internal_max_Nm']<=5.;rows.append(dict(index=int(index),t=float(z['t'][index]),accepted=bool(accepted),native_contact_root=native_contact_root,source_domain_valid=not bool(np.any(bank.bounds)),solver_success=bool(success),errors=error,residual=residual.tolist(),contacts=contacts,inferred_contact_normal_bw=float(reaction.sum()),max_activation=float(alpha.max()),native_guide_generalized_passive=d.qfrc_passive[:6].tolist()));alphas.append(alpha);prior=alpha
 report=dict(accepted_static_inverse_upper_bound=bool(all(r['accepted'] for r in rows)),accepted_native_contact_root=bool(all(r['native_contact_root']['accepted'] for r in rows)),manifest=manifest,independent_coordinates=[m.joint(int(m.dof_jntid[i])).name for i in drivers],rows=rows);(out/'SOURCE_STATIC_INVERSE.json').write_text(json.dumps(report,indent=2));np.savez_compressed(out/'source_static_activation.npz',t=z['t'][indices],indices=indices,activation=alphas,muscles=names);print(json.dumps(dict(accepted_static_inverse_upper_bound=report['accepted_static_inverse_upper_bound'],accepted_native_contact_root=report['accepted_native_contact_root'],accepted_samples=sum(r['accepted'] for r in rows),samples=len(rows),worst_errors={k:max(r['errors'][k] for r in rows) for k in rows[0]['errors']},source_domain_failures=sum(not r['source_domain_valid'] for r in rows)),indent=2))

def actual_torque_audit(model_folder,folder,out):
 """Observed loaded dynamics plus conditional source-curve recruitment capacity."""
 from mechanics import CoupledArms,JOINTS
 b=np.load(folder/'backend_trace.npz');z={k:b[k] for k in b.files};c=np.load(folder/'controller.npz');offset=len(c['t'])-len(z['t'])
 if not np.array_equal(c['qfull'][offset:],z['qpos']):raise RuntimeError('Controller/physical trace alignment failed')
 logical=c['t'][offset:];m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m);arms=CoupledArms(m);names=z['muscle_names'].tolist();p=parameters_from_osim(OSIM,names,0.);dt=float(np.median(np.diff(z['t'])));components=np.zeros((3,m.nu));M=np.zeros((m.nv,m.nv));rows=[];maxerr=dict(source_force_N=0.,actuator_projection_Nm=0.,generalized_balance=0.,constraint_split=0.);bw=float(m.body_mass.sum()*9.81)
 families=dict(hip_flexor_mono=('psoas','iliacus'),hip_extensor_mono=('glut_max1','glut_max2','glut_max3'),knee_extensor_mono=('vas_lat','vas_med','vas_int'),knee_flexor_mono=('bifemsh',),HAM=('semimem','semiten','bifemlh'),RF=('rect_fem',),GAS=('med_gas','lat_gas'),SAR=('sar',),GRAC=('grac',),TFL=('tfl',))
 masks={(family,s):np.array([i for i,n in enumerate(names) if n.endswith('_'+s) and n[:-2] in members],dtype=int) for family,members in families.items() for s in ('r','l')};window=np.flatnonzero((logical>=5.)&(logical<13.));indices=window[np.linspace(0,len(window)-1,101).astype(int)]
 for i in indices:
  d.qpos[:]=z['qpos'][i];d.qvel[:]=z['qvel'][i];mujoco.mj_forward(m,d);rho=arms.update(d).copy();E=arms.E;alpha=z['activation'][i-1];fiber_components(z['fiber_state'][i],z['fiber_state'][i-1],alpha,p,dt,components);maxerr['source_force_N']=max(maxerr['source_force_N'],float(np.max(abs(components.sum(axis=0)-z['force'][i]))));mujoco.mj_fullM(m,M,d.qM)
  external=z['qfrc_smooth'][i]-z['qfrc_passive'][i]+z['qfrc_bias'][i]-z['qfrc_actuator'][i];inertial=M@z['qacc'][i];balance=inertial-z['qfrc_smooth'][i]-z['qfrc_constraint'][i];maxerr['generalized_balance']=max(maxerr['generalized_balance'],float(np.max(abs(balance))));maxerr['constraint_split']=max(maxerr['constraint_split'],float(np.max(abs(z['qfrc_constraint'][i]-z['qfrc_contact'][i]-z['qfrc_limit'][i]-z['qfrc_equality'][i]))));maxerr['actuator_projection_Nm']=max(maxerr['actuator_projection_Nm'],float(np.max(abs(rho@z['force'][i]-E.T@z['qfrc_actuator'][i]))))
  projected={key:E.T@z['qfrc_'+key][i] for key in ('bias','passive','contact','limit','equality','actuator')};projected['external']=E.T@external;projected['inertial']=E.T@inertial;required=projected['bias']-projected['passive']-projected['contact']-projected['limit']-projected['equality']-projected['external'];record=dict(logical_time_s=float(logical[i]),physical_time_s=float(z['t'][i]),sides={})
  for s in ('r','l'):
   hr=JOINTS.index('hip_flexion_'+s);kr=JOINTS.index('knee_angle_'+s);ar=JOINTS.index('ankle_angle_'+s);group={family:dict(hip_Nm=float((rho[hr,ids]*z['force'][i,ids]).sum()),knee_Nm=float((rho[kr,ids]*z['force'][i,ids]).sum()),ankle_Nm=float((rho[ar,ids]*z['force'][i,ids]).sum()),passive_hip_Nm=float((rho[hr,ids]*components[1,ids]).sum()),passive_knee_Nm=float((rho[kr,ids]*components[1,ids]).sum()),mean_activation=float(alpha[ids].mean())) for (family,side),ids in masks.items() if side==s}
   fiber_length=np.hypot(z['fiber_state'][i],p[:,3]);cosine=z['fiber_state'][i]/fiber_length;velocity=(z['fiber_state'][i]-z['fiber_state'][i-1])*cosine/(dt*p[:,1]);isometric=p[:,0]*np.maximum(.01,alpha)*np.exp(-((fiber_length/p[:,1])-1.)**2/.5)*cosine
   for (family,side),ids in masks.items():
    if side==s:group[family].update(mean_fiber_velocity_l0_s=float(velocity[ids].mean()),mean_normalized_fiber_length=float((fiber_length[ids]/p[ids,1]).mean()),mean_active_FV_multiplier=float(np.mean(components[0,ids]/isometric[ids])))
   ids=masks['knee_extensor_mono',s];bank=DynamicThelenBank(p[ids],dt);L=d.actuator_length[ids].copy();a=alpha[ids].copy();cap=max(0.,.98-float(a.max()));actual_quad=float(rho[kr,ids]@z['force'][i,ids]);desired=float(required[kr]-(projected['actuator'][kr]-actual_quad))
   def quad(delta):
    bank.initialize(L,a+delta);return float(rho[kr,ids]@bank.force),bool(np.any(bank.bounds))
   base,basebound=quad(0.);maximum,maxbound=quad(cap);lo,hi=0.,cap;feasible=desired<=maximum
   if feasible and desired>base:
    for _ in range(24):
     mid=(lo+hi)*.5;value,bound=quad(mid)
     if value<desired:lo=mid
     else:hi=mid
   delta=hi if feasible and desired>base else (0. if desired<=base else None)
   record['sides'][s]=dict(native_hip_knee_ankle_deg=[float(np.rad2deg(d.qpos[m.joint(j+'_'+s).qposadr[0]])) for j in ('hip_flexion','knee_angle','ankle_angle')],native_knee_velocity_rad_s=float(d.qvel[m.joint('knee_angle_'+s).dofadr[0]]),knee_torque_Nm={key:float(v[kr]) for key,v in projected.items()},hip_torque_Nm={key:float(v[hr]) for key,v in projected.items()},required_zero_accel_muscle_knee_Nm=float(required[kr]),required_zero_accel_muscle_hip_Nm=float(required[hr]),groups=group,conditional_quadriceps=dict(actual_dynamic_Nm=actual_quad,required_Nm=desired,source_static_current_activation_Nm=base,source_static_max_activation_Nm=maximum,extra_uniform_activation=delta,within_capacity=bool(feasible),source_domain_valid=not(basebound or maxbound),scope='Observed contact/other-muscle torques held fixed; original source static vasti equilibria,uniform added activation up to.98. Conditional capacity only,not a forward controller fit.'))
  rows.append(record)
 if maxerr['source_force_N']>1e-5 or maxerr['actuator_projection_Nm']>1e-4:raise RuntimeError('Source dynamics consistency failed '+str(maxerr))
 summary={}
 for s in ('r','l'):
  rr=[row['sides'][s] for row in rows];loaded=[r for r in rr if abs(r['knee_torque_Nm']['contact'])>5.];summary[s]=dict(samples=len(rr),loaded_contact_torque_samples=len(loaded),loaded_knee_contact_Nm_percentiles=np.percentile([r['knee_torque_Nm']['contact'] for r in loaded],[5,50,95]).tolist(),loaded_actual_muscle_knee_Nm_percentiles=np.percentile([r['knee_torque_Nm']['actuator'] for r in loaded],[5,50,95]).tolist(),loaded_static_needed_muscle_knee_Nm_percentiles=np.percentile([r['required_zero_accel_muscle_knee_Nm'] for r in loaded],[5,50,95]).tolist(),loaded_quad_capacity_failures=sum(not r['conditional_quadriceps']['within_capacity'] for r in loaded),loaded_extra_quad_activation_percentiles=np.percentile([r['conditional_quadriceps']['extra_uniform_activation'] for r in loaded if r['conditional_quadriceps']['extra_uniform_activation'] is not None],[5,50,95]).tolist(),loaded_group_median={f:{key:float(np.median([r['groups'][f][key] for r in loaded])) for key in ('hip_Nm','knee_Nm','passive_hip_Nm','passive_knee_Nm','mean_activation')} for f in families})
 os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
 import opensim as osim
 om=osim.Model(str(OSIM));om.setUseVisualizer(False);state=om.initSystem();coords=om.getCoordinateSet();muscles=om.getMuscles();signs={r['coordinate']:r['sign'] for r in json.loads((HERE.parent/'hunt_training/runs/ik_audit/fk_signs.json').read_text())};static_bank=DynamicThelenBank(p,dt);source_rows=[]
 for i in indices[np.linspace(0,len(indices)-1,9).astype(int)]:
  d.qpos[:]=z['qpos'][i];mujoco.mj_forward(m,d);alpha=np.maximum(.01,z['activation'][i-1])
  for name in JOINTS:coords.get(name).setValue(state,float(d.qpos[m.joint(name).qposadr[0]]/signs[name]),False)
  for k,name in enumerate(names):muscles.get(name).setActivation(state,float(alpha[k]))
  om.equilibrateMuscles(state);om.realizeDynamics(state);source_lengths=np.array([muscles.get(name).getLength(state) for name in names]);source_force=np.array([muscles.get(name).getTendonForce(state) for name in names]);static_bank.initialize(d.actuator_length,alpha)
  source_rows.append(dict(logical_time_s=float(logical[i]),mujoco_minus_opensim_Lmt_m=(d.actuator_length-source_lengths).tolist(),mujoco_length_source_static_minus_original_opensim_force_N=(static_bank.force-source_force).tolist(),opensim_static_fiber_length_l0=[float(muscles.get(name).getFiberLength(state)/p[k,1]) for k,name in enumerate(names)],source_domain_valid=not bool(np.any(static_bank.bounds))))
 source_check=dict(scope='Nine actual full17-coordinate poses mapped using previously validated FK signs. Original OpenSim model static equilibrium at identical activation versus source constitutive law on converted absolute route lengths; no trajectory or source-length correction.',muscles=names,rows=source_rows)
 out.mkdir(parents=True,exist_ok=True);report=dict(scope='Frozen incumbent diagnostic replay. Exact captured physical model and equality-aware17 anatomical coordinates. Positive hip torque is native extension; positive knee torque is extension. No contact or trajectory forces applied by this audit; conditional static capacity is not gait validation.',source_folder=str(folder),logical_warmup_frames=int(offset),consistency_errors=maxerr,loaded_definition='abs(actual equality-aware knee contact torque)>5Nm,not a GRF acceptance gate',summary=summary,original_opensim_pose_check=source_check,rows=rows);(out/'ACTUAL_STANCE_TORQUE_AUDIT.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(consistency_errors=maxerr,summary=summary),indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--folder',type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--reference-model',type=Path);p.add_argument('--prepare-contact-reference',action='store_true');p.add_argument('--actual-torque-audit',action='store_true');p.add_argument('--teacher-inverse',action='store_true');p.add_argument('--native-stance-recruitment',action='store_true');p.add_argument('--hold-dt',type=float);p.add_argument('--hold-physics',choices=('rigid','compliant'),default='compliant');p.add_argument('--preserve-pruned',action='store_true');p.add_argument('--activation-floor',type=float,default=.01);p.add_argument('--transition-start',type=Path);p.add_argument('--transition-end',type=Path);p.add_argument('--transition-via',type=Path);p.add_argument('--transition-roundtrip',action='store_true');p.add_argument('--return-dwell',type=float,default=0.);p.add_argument('--landing-duration-scale',type=float,default=1.);p.add_argument('--export-training-reference',action='store_true');p.add_argument('--root-aware-cartesian',action='store_true');p.add_argument('--spec',type=Path);p.add_argument('--transition-ramp',type=float,default=1.);p.add_argument('--transition-feedback',action='store_true');p.add_argument('--command-filter',type=float,default=0.);p.add_argument('--feedback-warmup',action='store_true');p.add_argument('--mtp-impedance-scaled',action='store_true');p.add_argument('--trace-stride',type=int,default=10);p.add_argument('--reference-override',type=Path);p.add_argument('--passive-sliding',action='store_true')
 p.add_argument('--feasibility-first',action='store_true')
 p.add_argument('--static-forward-acceleration',action='store_true')
 p.add_argument('--loaded-cartesian-followers',action='store_true')
 p.add_argument('--cartesian-tangential-clearance',type=float)
 p.add_argument('--native-impedance-acceleration',action='store_true')
 p.add_argument('--stationary-load-split',action='store_true')
 p.add_argument('--root-geometry-prefit',action='store_true')
 p.add_argument('--root-geometry-only',action='store_true')
 p.add_argument('--initialize-actual-lumbar',action='store_true')
 p.add_argument('--capacity-aware-lumbar',action='store_true')
 p.add_argument('--predictive-source-horizon',type=float)
 p.add_argument('--actuation-audit-time',type=float,action='append')
 p.add_argument('--stance-posture-cartesian',action='store_true');p.add_argument('--native-acceleration-teacher',action='store_true');p.add_argument('--native-root-augmented-torque',action='store_true');p.add_argument('--no-stop-compensation',action='store_true');p.add_argument('--landing-audit',action='store_true');p.add_argument('--world-foot-pitch',action='store_true');p.add_argument('--full-native-equilibrium',action='store_true');p.add_argument('--max-static-evaluations',type=int,default=80);p.add_argument('--resume-equilibrium',type=Path);p.add_argument('--static-force-priority',type=float,default=1.);a=p.parse_args()
 if a.transition_start:
  if a.hold_physics!='compliant':p.error('Transition mode is source-compliant only; --hold-physics applies to static holds')
  if a.activation_floor!=.01:p.error('Transition allocator currently requires the explicit source1% floor; no5% result is implied')
  if a.preserve_pruned:
   saved=np.load(a.transition_start/'source_recruitment.npz');mask=saved['force_channel_enabled'] if 'force_channel_enabled' in saved else np.ones(len(saved['muscles']),bool);expected={family+'_'+side for family in ('quad_fem','gem','peri') for side in ('r','l')};actual={str(n) for n,e in zip(saved['muscles'],mask) if not e}
   if actual!=expected:p.error('--preserve-pruned requires saved endpoints with exactly the six legacy force channels disabled; source92 endpoints cannot silently satisfy this request')
  native_transition(a.reference_model,a.transition_start,a.transition_end,a.spec,a.out,a.transition_ramp,stride=a.trace_stride,feedback=a.transition_feedback,command_filter=a.command_filter,feedback_warmup=a.feedback_warmup,via_folder=a.transition_via,mtp_impedance_scaled=a.mtp_impedance_scaled,physics_dt=a.hold_dt,roundtrip=a.transition_roundtrip,return_dwell=a.return_dwell,export_training_reference=a.export_training_reference,landing_duration_scale=a.landing_duration_scale,cartesian_root=a.root_aware_cartesian,stance_posture_cartesian=a.stance_posture_cartesian,native_acceleration=a.native_acceleration_teacher or a.native_impedance_acceleration,native_root_augment=a.native_root_augmented_torque,compensate_coordinate_stops=not a.no_stop_compensation,world_foot_pitch=a.world_foot_pitch,loaded_cartesian_followers=a.loaded_cartesian_followers,cartesian_tangential_clearance=a.cartesian_tangential_clearance,mass_consistent_impedance=a.native_impedance_acceleration,predictive_source_horizon=a.predictive_source_horizon,actuation_audit_times=a.actuation_audit_time)
 elif a.full_native_equilibrium:full_native_equilibrium(a.reference_model,a.folder,a.spec,a.out,a.activation_floor,a.max_static_evaluations,a.resume_equilibrium,a.static_force_priority,a.feasibility_first,a.static_forward_acceleration,a.stationary_load_split,a.root_geometry_prefit,a.root_geometry_only,a.initialize_actual_lumbar,a.capacity_aware_lumbar)
 elif a.landing_audit:audit_landing_reference(a.reference_model,a.folder,a.out)
 elif a.native_stance_recruitment:native_stance_recruitment(a.reference_model,a.folder,a.out,a.hold_dt,a.hold_physics,a.preserve_pruned,a.activation_floor)
 elif a.teacher_inverse:teacher_inverse(a.reference_model,a.folder,a.out)
 elif a.actual_torque_audit:actual_torque_audit(a.reference_model,a.folder,a.out)
 elif a.reference_model and a.prepare_contact_reference:prepare_contact_reference(a.reference_model,a.out)
 elif a.reference_model:reference_audit(a.reference_model,a.out,a.reference_override,a.passive_sliding)
 elif a.folder:compare(a.folder,a.out)
 else:p.error('Provide --folder or --reference-model')
