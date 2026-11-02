"""Existing-trace geometry/contact-solve audit; no time integration or fitting."""
import json,argparse,sys,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull,QhullError
import cad_pair_support_step as S
C=S.N.C;R=C.R;mujoco=C.mujoco
OUT=S.ROOT.parent/'posture_balance_audit'
def load(p):
 with np.load(p) as z:return {k:z[k] for k in z.files}
def model():
 m=mujoco.MjModel.from_xml_path(str(C.ROBOT/'ah001_neural.xml'));m.actuator_gainprm[:]=0;m.actuator_biasprm[:]=0;m.actuator_gaintype[:]=0;m.actuator_biastype[:]=0
 return m,C.DynamicRoutes(m,json.loads((C.ROBOT/'robotbody_routes.json').read_text()))
def probe():
 OUT.mkdir(exist_ok=True);root=S.ROOT/'mixed/mixed';a=load(root/'ACTUATOR_BRIDGE.npz');b=load(root/'cad.npz');m,dyn=model();errors={name:[] for name in ('previous','current','zero')};rows=[]
 for j in np.linspace(20,len(b['time'])-1,16,dtype=int):
  k=round(float(b['time'][j])/.002)
  for variant in errors:
   d=mujoco.MjData(m);d.qpos[:]=a['qpos_full'][k];d.qvel[:]=a['qvel_full'][k];d.time=float(a['t'][k]);p=max(k-1,0)
   if variant=='previous':d.qfrc_applied[:]=a['qfrc'][p];d.xfrc_applied[:]=a['rig'][p]
   if variant=='current':d.qfrc_applied[:]=a['qfrc'][k];d.xfrc_applied[:]=a['rig'][k]
   if variant=='zero':d.xfrc_applied[:]=a['rig'][k]
   dyn.update(d);computed=R.contact_loads(m,d);error=float(np.max(np.abs(computed-b['contact_loads'][j])));errors[variant].append(error)
   rows.append(dict(t=float(b['time'][j]),variant=variant,max_channel_error_N=error,reported=b['contact_loads'][j].tolist(),reconstructed=computed.tolist()))
 result=dict(machine='easteregg2',scope='Read-only mj_forward snapshots. No mj_step or newtrajectory. Original percontactwrenches/positions were not saved; parity checked against allfour archivedheel/toe loadchannels.',error_max_N={k:max(v) for k,v in errors.items()},rows=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (OUT/'CONTACT_RECONSTRUCTION_PROBE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result['error_max_N'],indent=2))

def archived(folder,name):
 full=folder/name
 if full.exists():return load(full),None
 manifest=json.loads((folder.parent/'TELEMETRY_MANIFEST.json').read_text())['archives'];index=next(i for i,e in enumerate(manifest) if Path(e['original_path']).name==name);z=load(folder.parent/'COMPACT_TELEMETRY.npz');prefix=f'a{index}_';data={k[len(prefix):]:v for k,v in z.items() if k.startswith(prefix) and not k.endswith('_sample_indices')};return data,'Sparse/startup/end retainedwitnesses only'

def inside_margin(points,p):
 try:
  hull=ConvexHull(points);vertices=points[hull.vertices];edges=np.roll(vertices,-1,axis=0)-vertices;cross=edges[:,0]*(p[1]-vertices[:,1])-edges[:,1]*(p[0]-vertices[:,0]);return float(np.min(cross/np.linalg.norm(edges,axis=1)))
 except QhullError:return np.nan

def geometry(m,d,dyn):
 pelvis=m.body('pelvis').id;trunk=m.body('lumbar_torso').id;jac=np.zeros((3,m.nv))
 if hasattr(mujoco,'mj_jacSubtreeCom'):mujoco.mj_jacSubtreeCom(m,d,jac,pelvis)
 else:
  for bid in range(1,m.nbody):
   temp=np.zeros_like(jac);angular=np.zeros_like(jac);mujoco.mj_jacBodyCom(m,d,temp,angular,bid);jac+=m.body_mass[bid]*temp
  jac/=m.body_mass.sum()
 com=d.subtree_com[pelvis].copy();vel=jac@d.qvel;omega=np.sqrt(9.81/max(com[1],1e-9));xcom=com[[0,2]]+vel[[0,2]]/omega
 def lean(bid):axis=d.xmat[bid].reshape(3,3)[:,1];return np.rad2deg(np.arctan2(axis[0],axis[1]))
 forces=[];positions=[];side=[];loads=R.contact_loads(m,d)
 for i in range(d.ncon):
  c=d.contact[i];names=[m.geom(int(g)).name for g in (c.geom1,c.geom2)]
  if 'ground' not in names or not any(n.startswith('foot_') for n in names):continue
  w=np.zeros(6);mujoco.mj_contactForce(m,d,i,w);F=c.frame.reshape(3,3).T@w[:3]
  if names[1]=='ground':F=-F
  if F[1]<=1e-7:continue
  forces.append(F);positions.append(c.pos.copy());side.append(0 if any(n.startswith('foot_r') for n in names) else 1)
 F=np.asarray(forces).reshape(-1,3);P=np.asarray(positions).reshape(-1,3);total=F.sum(axis=0);moment=np.cross(P,F).sum(axis=0)
 cop=np.array([moment[2],-moment[0]])/total[1] if total[1]>1 else np.full(2,np.nan);centroid=(P[:,[0,2]]*F[:,1,None]).sum(axis=0)/total[1] if total[1]>1 else np.full(2,np.nan);footcop=np.full((2,2),np.nan)
 for si in range(2):
  take=np.array(side)==si
  if take.any() and F[take,1].sum()>1:footcop[si]=(P[take][:,[0,2]]*F[take,1,None]).sum(axis=0)/F[take,1].sum()
 polygon=P[:,[0,2]];margin=inside_margin(polygon,com[[0,2]]) if len(P)>=3 else np.nan;xmargin=inside_margin(polygon,xcom) if len(P)>=3 else np.nan
 p=d.xipos[pelvis];rig=d.xfrc_applied[pelvis];gravity_mz=-S.BW*(com[0]-cop[0]) if np.isfinite(cop[0]) else np.nan;vertical_mz=(p[0]-cop[0])*rig[1];ap_mz=-p[1]*rig[0];rot_mz=rig[5]
 return dict(COM=com,COM_velocity=vel,XCoM_xz=xcom,CoP_xz=cop,normal_centroid_xz=centroid,foot_normal_centroid_xz=footcop,ground_force=total,contact_AP_bounds=np.array([polygon[:,0].min(),polygon[:,0].max()]) if len(P) else np.full(2,np.nan),COM_polygon_margin=np.array(margin),XCoM_polygon_margin=np.array(xmargin),pelvis_forelean_deg=np.array(lean(pelvis)),trunk_forelean_deg=np.array(lean(trunk)),pelvis_rotation_deg=np.rad2deg(R.rotation_error(d.qpos[3:7])),lumbar_rotation_deg=np.rad2deg(R.rotation_error(d.qpos[m.jnt_qposadr[m.joint('lumbar_ball').id]:m.jnt_qposadr[m.joint('lumbar_ball').id]+4])),hip_rotation_deg=np.array([np.rad2deg(R.rotation_error(d.qpos[m.jnt_qposadr[m.joint('thigh_'+s+'_joint').id]:m.jnt_qposadr[m.joint('thigh_'+s+'_joint').id]+4])) for s in ('r','l')]),knee_deg=np.rad2deg(list(dyn.coordinates(d).values())),rig=rig.copy(),rig_application_point=p.copy(),sagittal_external_moments_Nm=np.array([gravity_mz,vertical_mz,ap_mz,rot_mz]),contact_loads=loads)

def audit():
 OUT.mkdir(exist_ok=True);m,dyn=model();results={};initial=None
 for label,base in [('parent',S.BASE),('support085',S.ROOT)]:
  for mode in ('reference','mixed'):
   folder=base/mode/mode;a,scope=archived(folder,'ACTUATOR_BRIDGE.npz');b,bscope=archived(folder,'cad.npz');byindex={round(float(t)/.002):i for i,t in enumerate(a['t'])};bindex={round(float(t)/.002):i for i,t in enumerate(b['time'])};rows=[];parity=[]
   wanted=sorted(k for k in byindex if k%5==0 and k>0 and k-1 in byindex)
   for k in wanted:
    i=byindex[k];prior=byindex[k-1];d=mujoco.MjData(m);d.qpos[:]=a['qpos_full'][i];d.qvel[:]=a['qvel_full'][i];d.qfrc_applied[:]=a['qfrc'][prior];d.xfrc_applied[:]=a['rig'][prior];dyn.update(d);row=geometry(m,d,dyn);row['t']=np.array(a['t'][i]);rows.append(row)
    if k in bindex:parity.append(float(np.max(np.abs(row['contact_loads']-b['contact_loads'][bindex[k]]))))
   if initial is None:
    d=mujoco.MjData(m);d.qpos[:]=a['qpos_full'][0];d.qvel[:]=0;dyn.update(d);initial={key:np.asarray(value).tolist() for key,value in geometry(m,d,dyn).items()};initial['qpos']=d.qpos.tolist();initial['posture_target']='Standingtarget=q0 fromsourceframe1069; no laterjointreference changes. Root translation/orientation PD terms are discarded byinternalallocator; horizontalrig anchorsinitialrootx, angularrig targetsworldidentity.'
   data={key:np.asarray([r[key] for r in rows]) for key in rows[0]};np.savez_compressed(OUT/(label+'_'+mode+'_GEOMETRY.npz'),**data);tail=data['t']>=6;good=tail&np.isfinite(data['CoP_xz'][:,0]);summary=dict(samples=len(rows),scope=scope or 'Fullbridgehistory sampled100Hz',contact_channel_parity_samples=len(parity),contact_channel_max_error_N=max(parity) if parity else None,window_t_s=[float(data['t'][good].min()),float(data['t'][good].max())],tail_samples=int(good.sum()))
   for key in ('COM','COM_velocity','CoP_xz','XCoM_xz','contact_AP_bounds','pelvis_forelean_deg','trunk_forelean_deg','pelvis_rotation_deg','lumbar_rotation_deg','hip_rotation_deg','knee_deg','rig','rig_application_point','ground_force','contact_loads','sagittal_external_moments_Nm','COM_polygon_margin','XCoM_polygon_margin'):
    summary['tail_mean_'+key]=np.nanmean(data[key][good],axis=0).tolist()
   summary['tail_mean_COM_minus_CoP_x_m']=float(np.mean(data['COM'][good,0]-data['CoP_xz'][good,0]));summary['tail_mean_XCoM_minus_CoP_x_m']=float(np.mean(data['XCoM_xz'][good,0]-data['CoP_xz'][good,0]));results[label+'_'+mode]=summary
 result=dict(machine='easteregg2',scope='Existing-state mj_forward only, no timeintegration/newbehavioraltrials. CoP/contactpolygon are reconstructionvalidated against savedheel/toechannels; sparsefailedcase has shorterwindow.',XCoM_caution='ClassicalCOM+COMvelocity/sqrt(g/COMheight), kinematicindicator only; strongvertical/AP/rotationalharness invalidates unsupportedLIPbalanceinterpretation.',initial_posture=initial,cases=results)
 (OUT/'POSTURE_GEOMETRY_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(results,indent=2))

def report():
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 import xml.etree.ElementTree as ET
 result=json.loads((OUT/'POSTURE_GEOMETRY_AUDIT.json').read_text());source=C.ROBOT.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim'
 ranges={}
 if source.exists():
  for q in ET.parse(source).iter('Coordinate'):
   values=q.findtext('range')
   if values:ranges[q.get('name')]=np.rad2deg([float(v) for v in values.split()]).tolist()
 m,_=model();limits={m.joint(i).name:dict(limited=bool(m.jnt_limited[i]),range=m.jnt_range[i].tolist()) for i in range(m.njnt)}
 (OUT/'SOURCE_ROM_AUDIT.json').write_text(json.dumps(dict(CAD_joint_limits=limits,OpenSim_source=str(source),OpenSim_ranges_deg=ranges,warning='CAD hips/knees/lumbar do not carry finite anatomicalstops. OpenSimEulercoordinates are not CADball rotationvector components; humanranges alone do not certify hardwareROM.'),indent=2))
 plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False});fig,axes=plt.subplots(4,2,figsize=(12,10),sharex=True,constrained_layout=True)
 for col,label in enumerate(('parent','support085')):
  for mode,style,color in [('reference','--','#EE7733'),('mixed','-','#0077BB')]:
   d=load(OUT/(label+'_'+mode+'_GEOMETRY.npz'));t=d['t'];axes[0,col].plot(t,1000*(d['COM'][:,0]-d['CoP_xz'][:,0]),style,color=color,label=mode);axes[1,col].plot(t,1000*d['COM_polygon_margin'],style,color=color,label=mode);axes[2,col].plot(t,d['trunk_forelean_deg'],style,color=color,label=mode);axes[3,col].plot(t,d['rig'][:,0],style,color=color,label=mode)
  axes[0,col].set_title(('90%','85%')[col]+' feedforward support')
  for i,text in enumerate(('COM anterior − CoP (mm)','Loaded support-polygon margin (mm)','Trunk forelean (+anterior, °)','Actual AP harness force (+anterior, N)')):axes[i,col].set_ylabel(text);axes[i,col].grid(alpha=.2);axes[i,col].legend(frameon=False,fontsize=8);axes[i,col].set_xlim(0,8);axes[i,col].axvspan(6,8,color='.85',alpha=.4)
  axes[1,col].axhline(0,color='.3',lw=.8);axes[3,col].set_xlabel('Time (s)')
 fig.suptitle('Existing-trace posture audit: parity-validated contact reconstruction\nCOM position and trunk lean are different quantities; 85% reference geometry available only in retained windows',fontsize=12)
 for ext in ('png','pdf'):fig.savefig(OUT/('POSTURE_CONTACT_AUDIT.'+ext),dpi=200)
 plt.close(fig)
 lines=['# Supported CAD posture and foot-contact audit','','Machine: easteregg2. Existing saved states only: no timeintegration, neural/controller fit, posture change, or additional support reduction. Contactwrenches are reconstructed using the previous-step actual muscle generalizedforce and apparatuswrench, matching the original sensor timing. All four saved heel/toe channels agree within4.1e−12N. The failed85%ideal case has only a short final7.75–7.99s retained geometrywindow; the other cases cover6.01–7.99s.','','| Case | COM−CoP AP (mm) | COM polygonmargin (mm) | Trunk forelean ° | AP harness N | Sagittal moments gravity / verticalrig / APrig / angularrig (Nm) |','| --- | ---: | ---: | ---: | ---: | --- |']
 for name,v in result['cases'].items():
  moments=np.array(v['tail_mean_sagittal_external_moments_Nm']);lines.append(f"| {name} | {1000*v['tail_mean_COM_minus_CoP_x_m']:.2f} | {1000*v['tail_mean_COM_polygon_margin']:.2f} | {v['tail_mean_trunk_forelean_deg']:.3f} | {v['tail_mean_rig'][0]:.3f} | {np.round(moments,3).tolist()} |")
  print(name,dict(moments=moments.tolist(),sum=float(moments.sum()),rigpoint=v['tail_mean_rig_application_point'],heeltoe=v['tail_mean_contact_loads']))
 lines+=['','The85%pressurecase COM lies near the posterior contact boundary and approximately167mm behind the weightedCoP. Its loaded2D polygonmargin is−7.22mm. However, the trunk is only0.443° posterior, compared with0.607° posterior initially. Thus the posterior COM excursion is mainly translation and legconfiguration, not a large backward trunkangle. CoP proximity is strongly altered by the harness unloading geometry. Classical XCoM is reported only as COM+velocity/sqrt(g/height); its unsupported pendulum balance interpretation does not apply under these large external forces/torques.','','The Cartesian harness force acts at the pelvis bodyCOM (data.xipos), not its bodyframeorigin (data.xpos). This distinction materially changes the sagittal moment arm. The first derived table used the frameorigin and is preserved as POSTURE_GEOMETRY_AUDIT_before_forcepoint_fix.json; contact/CoP results are unchanged. The correctedmoment components use actualbodyCOM and are diagnostic, not an unsupported balance claim.','','The standing internaljoint reference is the initial sourceframe1069 pose. Free-root PD terms are discarded by the internal allocator; a separate horizontalrig spring anchors initialrootx, while angularrig targets worldidentity. The initialreference was selected from a route-envelope recording, not established as a ground-standingIK target. The CAD model has no explicit finite hip/knee/lumbar stops; the separate knee safety guard rejects hyperextension. SourceROM_AUDIT.json lists this limitation and the current OpenSimranges without equating rotationvectors to Eulercoordinates.','','A more anterior reference therefore requires a distinct offline kinematic/force-capacity feasibility check before any new closed-loop condition. Existing settling passes are not revised or relabeled.']
 (OUT/'EXPERIMENT_RECORD.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def candidate():
 from scipy.optimize import least_squares,linprog
 from scipy.spatial.transform import Rotation
 OUT.mkdir(exist_ok=True);registration=dict(type='Offlinegeometryandforce-capacity proposal only, not a newclosed-loopcondition',requested_anterior_COM_shift_m=.010,preserve_both_foot_body_poses=True,trunk_target='Remove fore/aftlean; retain remaining worldorientation',local_rotation_component_bounds_deg=5.,root_translation_component_bounds_m=.04,knee_search_interval_rad=[-.65,-.02],pair_weights_unchanged=S.N.U.weights(),pair_force_assumption='Freeze observed final neuraloutput referenceT=4530*u; no newneuralprediction or currentfit',other52='Knownidealreferencecapacity bounds .05*MIF..MIF, not calibratedBPAs',apparatus='Existing85%feedforward, originalheight/x/zanchors and allguides unchanged',no_timeintegration=True,no_hardware_ROM_certification=True,max_solver_nfev=40,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());path=OUT/'OFFLINE_REFERENCE_PREREGISTRATION.json'
 registration['reference_target']='Actual saved footposes and COM; exact kneeprojection used only as IK initialization. Prior projected-target draft preserved.'
 if path.exists():assert json.loads(path.read_text())==registration
 else:path.write_text(json.dumps(registration,indent=2))
 m,dyn=model();a=load(S.ROOT/'mixed/mixed/ACTUATOR_BRIDGE.npz');raw=a['qpos_full'][-1];d=mujoco.MjData(m);d.qpos[:]=raw;dyn.update(d);feet=[m.body('foot_'+s).id for s in ('r','l')];rawfeet=np.array([d.xpos[i].copy() for i in feet]);rawfootrot=[Rotation.from_matrix(d.xmat[i].reshape(3,3)) for i in feet];rawcom=d.subtree_com[m.body('pelvis').id].copy();rawtrunk=Rotation.from_matrix(d.xmat[m.body('lumbar_torso').id].reshape(3,3));rawaxis=d.xmat[m.body('lumbar_torso').id].reshape(3,3)[:,1];rawlean=np.arctan2(rawaxis[0],rawaxis[1]);originalk=np.array(list(dyn.coordinates(d).values()))
 for side,angle in zip(('r','l'),originalk):R.knee_pose(m,d,side,float(angle),d.qpos[[m.jnt_qposadr[m.joint(p+'_'+side+'_joint').id] for p in ('back','shank','front')]])
 dyn.update(d);base=d.qpos.copy();kbase=np.array(list(dyn.coordinates(d).values()));projectedfeet=np.array([d.xpos[i].copy() for i in feet]);footpos=rawfeet;footrot=rawfootrot;pelvis=m.body('pelvis').id;torso=m.body('lumbar_torso').id;com0=rawcom;trunkgoal=Rotation.from_rotvec([0,0,rawlean])*rawtrunk;goal=com0.copy();goal[0]+=.010
 hipdofs=[int(m.jnt_dofadr[m.joint('thigh_'+s+'_joint').id]) for s in ('r','l')];lumbar=int(m.jnt_dofadr[m.joint('lumbar_ball').id]);seeds=[base[[m.jnt_qposadr[m.joint(p+'_'+s+'_joint').id] for p in ('back','shank','front')]] for s in ('r','l')]
 def pose(x):
  d.qpos[:]=base;v=np.zeros(m.nv);v[:6]=x[:6]
  for j,di in enumerate(hipdofs):v[di:di+3]=x[6+3*j:9+3*j]
  v[lumbar:lumbar+3]=x[14:17];mujoco.mj_integratePos(m,d.qpos,v,1.)
  for j,side in enumerate(('r','l')):R.knee_pose(m,d,side,float(kbase[j]+x[12+j]),seeds[j])
  d.qvel[:]=0;dyn.update(d)
 def residual(x):
  pose(x);e=[]
  for i,bid in enumerate(feet):e.extend((d.xpos[bid]-footpos[i])/.001);e.extend((footrot[i].inv()*Rotation.from_matrix(d.xmat[bid].reshape(3,3))).as_rotvec()/np.deg2rad(1.))
  e.extend((d.subtree_com[pelvis][[0,2]]-goal[[0,2]])/.001);e.extend((trunkgoal.inv()*Rotation.from_matrix(d.xmat[torso].reshape(3,3))).as_rotvec()/np.deg2rad(1.));return np.r_[e,1e-5*x]
 limit=np.deg2rad(5);lo=np.full(17,-limit);hi=np.full(17,limit);lo[:3]=-.04;hi[:3]=.04
 for j in range(2):lo[12+j]=max(-limit,-.65-kbase[j]);hi[12+j]=min(limit,-.02-kbase[j])
 sol=least_squares(residual,np.zeros(17),bounds=(lo,hi),max_nfev=40,ftol=1e-10,xtol=1e-10,gtol=1e-10);pose(sol.x);q=d.qpos.copy();geo=geometry(m,d,dyn);poserr=max(float(np.linalg.norm(d.xpos[bid]-footpos[i])) for i,bid in enumerate(feet));roterr=max(float(np.linalg.norm((footrot[i].inv()*Rotation.from_matrix(d.xmat[bid].reshape(3,3))).as_rotvec())) for i,bid in enumerate(feet));geom_ok=bool(poserr<1e-5 and roterr<1e-4 and abs(d.subtree_com[pelvis][0]-goal[0])<1e-5)
 # Capacity equilibrium, including root balance and conservativefriction pyramid.
 config=json.loads((S.N.BASE/'ASSUMED_PAIR_CONFIGURATION.json').read_text());native=mujoco.MjModel.from_xml_path(str(C.ROBOT/'ah001_neural.xml'));caps=-native.actuator_gainprm[:,0].copy();ids=[m.actuator(n).id for n in C.PAIR];tids=m.actuator_trnid[:,0];pairT=caps[ids]*a['neural_output'][-1,ids];Nint=R.anatomical_knee_tangent(m,d);T=np.zeros((m.nv,6+Nint.shape[1]));T[:6,:6]=np.eye(6);T[6:,6:]=Nint;A=T.T@dyn.moment_arms[tids].T;positions=[];mus=[];contact_bodies=[]
 for c in d.contact:
  names=[m.geom(int(g)).name for g in (c.geom1,c.geom2)]
  if 'ground' in names and any(n.startswith('foot_') for n in names):positions.append(c.pos.copy());mus.append(float(c.friction[0])/np.sqrt(2));contact_bodies.append(int(m.geom_bodyid[c.geom2 if names[0]=='ground' else c.geom1]))
 anchors=a['qpos_full'][0];weight=m.body_mass.sum()*9.81;F=np.array([-20*(q[0]-anchors[0]),np.clip(.85*weight+800*(anchors[1]-q[1]),0,1.2*weight),-250*(q[2]-anchors[2])]);torque=-300*R.rotation_error(q[3:7]);Qrig=np.zeros(m.nv);mujoco.mj_applyFT(m,d,F,torque,d.xipos[pelvis],pelvis,Qrig);columns=[A]
 for P,bid in zip(positions,contact_bodies):
  J=np.zeros((3,m.nv));Jr=np.zeros_like(J);mujoco.mj_jac(m,d,J,Jr,P,bid);columns.append(T.T@J.T)
 E=np.column_stack(columns);bounds=[(.05*cap,cap) for cap in caps]
 for j,idx in enumerate(ids):bounds[idx]=(float(pairT[j]),float(pairT[j]))
 bounds.extend([(None,None),(0,weight),(None,None)]*len(positions));ineq=[]
 for j,mu in enumerate(mus):
  for axis in (0,2):
   for sign in (-1,1):row=np.zeros(E.shape[1]);row[m.nu+3*j+axis]=sign;row[m.nu+3*j+1]=-mu;ineq.append(row)
 objective=np.r_[1/caps,np.zeros(3*len(positions))];lp=linprog(objective,A_ub=np.asarray(ineq),b_ub=np.zeros(len(ineq)),A_eq=E,b_eq=T.T@(d.qfrc_bias-Qrig),bounds=bounds,method='highs') if positions else None
 bank=C.RoutedBPABank(C.PAIR,config,m.opt.timestep);pressure=[]
 for j,name in enumerate(C.PAIR):
  p,flag=C.inverse_pressure(float(pairT[j]),float(d.ten_length[tids[ids[j]]]+config['actuators'][name]['constant_remote_path_offset_m']),bank.parameters[j]);pressure.append(dict(name=name,T_N=float(pairT[j]),static_pressure_kpa=float(p),flag=int(flag)))
 result=dict(scope=registration['type'],geometric_candidate=geom_ok,solver_success=bool(sol.success),solver_nfev=int(sol.nfev),foot_position_max_error_m=poserr,foot_orientation_max_error_rad=roterr,initial_knee_projection_foot_shift_m=np.linalg.norm(projectedfeet-rawfeet,axis=1).tolist(),achieved_COM_shift_m=(d.subtree_com[pelvis]-com0).tolist(),local_parameters=sol.x.tolist(),knee_deg=np.rad2deg(list(dyn.coordinates(d).values())).tolist(),trunk_forelean_deg=float(geo['trunk_forelean_deg']),harness_force_N=F.tolist(),harness_torque_Nm=torque.tolist(),capacity_feasible=bool(lp is not None and lp.success),capacity_message=lp.message if lp is not None else 'No geometricgroundcontacts',pair_static_pressure=pressure,hardware_ROM_certified=False,no_neural_or_joint_gain_retune=True,no_closed_loop_execution=True)
 effective=weight-F[1];point=d.xipos[pelvis];com=d.subtree_com[pelvis];requiredcop=np.array([weight*com[0]-point[0]*F[1]+point[1]*F[0]-torque[2],weight*com[2]+point[1]*F[2]-point[2]*F[1]+torque[0]])/effective
 result.update(required_static_CoP_xz=requiredcop.tolist(),required_CoP_polygon_margin_m=inside_margin(np.asarray(positions)[:,[0,2]],requiredcop) if len(positions)>=3 else None,actual_contact_bodies=contact_bodies)
 if lp is not None and lp.success:result.update(equilibrium_residual_max=float(np.max(np.abs(E@lp.x-T.T@(d.qfrc_bias-Qrig)))),other52_total_tension_N=float(lp.x[:m.nu].sum()-pairT.sum()),contact_forces=lp.x[m.nu:].reshape(-1,3).tolist())
 np.savez_compressed(OUT/'OFFLINE_ANTERIOR_REFERENCE.npz',source_saved_qpos=raw,projected_qpos=base,candidate_qpos=q,local_delta=sol.x,pair_reference_force_N=pairT,capacity_solution=np.array([]) if lp is None or not lp.success else lp.x);(OUT/'OFFLINE_ANTERIOR_REFERENCE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--probe',action='store_true');p.add_argument('--audit',action='store_true');p.add_argument('--report',action='store_true');p.add_argument('--candidate',action='store_true');a=p.parse_args()
 if a.probe:probe()
 if a.audit:audit()
 if a.report:report()
 if a.candidate:candidate()
