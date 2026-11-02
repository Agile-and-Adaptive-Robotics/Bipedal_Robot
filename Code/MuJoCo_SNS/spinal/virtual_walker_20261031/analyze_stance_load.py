"""Read saved force/Jacobian/MN budgets before changing stance recruitment."""
import argparse,json
from pathlib import Path
import numpy as np
from import_virtual_walker import mujoco
import native_kernel as K
from recruitment import apply_recruitment

p=argparse.ArgumentParser();p.add_argument('folders',nargs='+');args=p.parse_args()
for item in args.folders:
 folder=Path(item);z=np.load(folder/'trajectory.npz');r=json.loads((folder/'RESULT.json').read_text());m=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));d=mujoco.MjData(m)
 bank=json.loads((folder.parent/'IMPORT_AUDIT.json').read_text())['muscles'];c=K.NativeController(r['source'],bank) if 'source' in r else K.NativeController(json.loads((folder/'SPIKING_COMPARISON.json').read_text())['source_configuration']['source'],bank)
 params=K.bank_parameters(bank);names=z['muscle_names'].tolist();neurons=z['neuron_names'].tolist();act=[i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];muscles=[bank[i]['name'] for i in act];N=len(z['t']);J=len(act)
 if (folder/'POOLED_FORCE_PARAMETERS.npz').exists():
  pooled=np.load(folder/'POOLED_FORCE_PARAMETERS.npz');assert pooled['muscle_names'].tolist()==names;params=pooled['pooled'].copy()
  for row in r['pooled_inherited_hip_physics']['legacy_hip_extensor_Ib_ledger']:
   for ii,(target,index,kind,coeff,body) in enumerate(c.inputs):
    if c.names[target]==row['target'] and kind==row['source_data_type']:c.inputs[ii]=(target,index,kind,np.array(row['polynomial_after_A_per_input']),body)
 tid=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+name) for name in muscles];joints=['hip_L','knee_L','ankle_L','hip_R','knee_R','ankle_R'];ji=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,name) for name in joints];qa=m.jnt_qposadr[ji];da=m.jnt_dofadr[ji];length=np.empty((N+1,J));jac=np.empty((N+1,J,6));tilt=np.empty(N+1);rootid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root')
 for i,q in enumerate(np.vstack([r['initial_qpos'],z['q']])):
  d.qpos[:]=q;mujoco.mj_kinematics(m,d);mujoco.mj_comPos(m,d);mujoco.mj_tendon(m,d);length[i]=d.ten_length[tid];jac[i]=d.ten_J.reshape(m.ntendon,m.nv)[np.ix_(tid,da)];R=d.xmat[rootid].reshape(3,3);tilt[i]=np.degrees(np.arctan2(R[0,1],R[1,1]))
 force=z['forces'][:,act];torque=-force[:,:,None]*jac[:-1];work=-force*np.diff(length,axis=0);bw=r['robot_mass_kg']*9.81;window=z['t']>=1.;legs={};rows=[]
 for side,index in [('L',0),('R',1)]:
  stance=window&(z['foot_vertical'][:,index]>.05*bw);single=stance&(z['foot_vertical'][:,1-index]<.025*bw)
  legs[side]=dict(stance_fraction=float(stance[window].mean()),single_support_fraction=float(single[window].mean()),single_support_mean_foot_BW=float(z['foot_vertical'][single,index].mean()/bw),mean_foot_BW=float(z['foot_vertical'][window,index].mean()/bw),peak_foot_BW=float(z['foot_vertical'][window,index].max()/bw))
 for col,idx in enumerate(act):
  name=bank[idx]['name'];side='L' if '_L_' in name else 'R';s=0 if side=='L' else 1;stance=window&(z['foot_vertical'][:,s]>.05*bw);pre,_,coeff=next(row for row in c.outputs if row[1]==idx);vm=K.polynomial(coeff,z['neuro'][:,pre]*.001);l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=params[idx];fl=np.maximum(0.,1.-((length[:-1,col]-l0)/width)**2);active=np.maximum(0.,B/(1.+np.exp(np.minimum(700.,C*(A-vm))))+D)*fl;full=np.maximum(0.,B/(1.+np.exp(min(700.,C*A)))+D)*fl;steady=np.clip((kpe*(length[:-1,col]-l0)+full)/(1+kpe/kse),0,maximum)
  rows.append(dict(muscle=name,actual_motor_neuron=neurons[pre],source_motor_rest_mV=float(c.params[pre,0]),source_force_sigmoid=dict(A_V=float(A),B_N=float(B),C_per_V=float(C),D_N=float(D)),source_force_feedback_adapters=[dict(neuron=c.names[target],data_type=kind,polynomial_A_per_input=coeff.tolist()) for target,index,kind,coeff,body in c.inputs if index==idx and kind in ('Ib','Tension')],stance_MN_mean_mV=float(z['neuro'][stance,pre].mean()),stance_force_mean_N=float(force[stance,col].mean()),stance_force_p95_N=float(np.percentile(force[stance,col],95)),stance_active_fraction_of_full_mean=float(np.mean(np.divide(active[stance],full[stance],out=np.zeros_like(active[stance]),where=full[stance]>0))),stance_full_command_quasistatic_force_mean_N=float(steady[stance].mean()),length_factor_mean=float(fl[window].mean()),length_factor_zero_fraction=float((fl[window]==0).mean()),positive_work_J=float(np.maximum(work[window,col],0).sum()),negative_work_J=float(-np.minimum(work[window,col],0).sum()),joint_torque_rms_Nm=np.sqrt(np.mean(torque[window,col]**2,axis=0)).tolist()))
 low_load=np.all(z['foot_vertical']<.025*bw,axis=1);double=np.all(z['foot_vertical']>.05*bw,axis=1);clear=np.stack([np.min(z['foot_clearance'][:,:2],axis=1),np.min(z['foot_clearance'][:,2:],axis=1)],axis=1)
 summary=dict(folder=str(folder),source_sha256=r['source_sha256'],net_support_BW=r['mean_net_vertical_assistance_BW'],both_feet_below_0p025BW_fraction=float(low_load[window].mean()),both_feet_above_0p05BW_fraction=float(double[window].mean()),both_whole_feet_above_1mm_fraction=float(np.all(clear[window]>.001,axis=1).mean()),mean_pelvis_anterior_tilt_deg=float(tilt[1:][window].mean()),pelvis_anterior_tilt_range_deg=[float(tilt[1:][window].min()),float(tilt[1:][window].max())],joint_order=joints,legs=legs,muscles=rows,scope='After 1 s. Source output-adapter mappings; actual applied Hill tension and pre-step tendon Jacobian. Quasistatic full-command estimate uses the same source Hill law at the observed route length with zero velocity, not a target force or a feasibility proof. Actual force includes dynamics/passive components. No dynamics or neural parameters changed.')
 (folder/'STANCE_LOAD_DIAGNOSIS.json').write_text(json.dumps(summary,indent=2));np.savez_compressed(folder/'STANCE_MECHANICS.npz',t=z['t'],muscles=muscles,joints=joints,torque_Nm=torque.sum(1),muscle_torque_Nm=torque,muscle_work_J=work,pelvis_anterior_tilt_deg=tilt[1:]);print(json.dumps({k:v for k,v in summary.items() if k not in ('muscles','scope')},indent=2))
 # Frozen-snapshot response fit. No state integration, output override, or applied force.
 meta=folder/'SPIKING_COMPARISON.json'
 if meta.exists() and len(neurons)==108:
  rec=json.loads(meta.read_text());cfg=rec['source_configuration'];c=K.NativeController(cfg['source'],bank,ankle_pf_crossed=cfg['ankle_pf_crossed']);apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale'])
  for edit in rec['spiking']['output_gain_edits']:
   edge=[e for e in c.edges if c.names[int(e[0])]==edit['source'] and c.names[int(e[1])]==edit['target']];assert len(edge)==1;edge[0][2]=edit['effective_gmax_uS'];edge[0][3]=edit.get('effective_threshold_mV',edge[0][3])
  activity=np.clip((z['neuro'][:,c.edges[:,0].astype(int)]-c.edges[:,3])/(c.edges[:,4]-c.edges[:,3]),0,1);fits=[];chemical=[]
  for post,target,coeff in c.outputs:
   if target not in act:continue
   side=0 if c.names[post].startswith('L ') else 1
   for k,e in enumerate(c.edges):
    if int(e[1])!=post:continue
    current=e[2]*activity[:,k]*(e[5]-z['neuro'][:,post]);phases={}
    for phase,mask in [('stance',window&(z['foot_vertical'][:,side]>.05*bw)),('swing',window&(z['foot_vertical'][:,side]<.025*bw))]:
     if mask.any():phases[phase]=dict(mean_nA=float(current[mask].mean()),minimum_nA=float(current[mask].min()),maximum_nA=float(current[mask].max()),depolarizing_fraction=float(np.mean(current[mask]>1e-9)),hyperpolarizing_fraction=float(np.mean(current[mask]<-1e-9)))
    chemical.append(dict(source=c.names[int(e[0])],target=c.names[post],gmax_uS=float(e[2]),threshold_mV=float(e[3]),saturation_mV=float(e[4]),reversal_mV=float(e[5]),phases=phases))
  (folder/'ACTUAL_MOTOR_CHEMICAL_CURRENTS.json').write_text(json.dumps(dict(rows=chemical,scope='All twelve actual muscle MNs, saved2ms proxy-voltage samples after1s, actual fitted gains/thresholds, source saturation and reversal. Positive means depolarizing at the observed postsynaptic voltage; a named inhibitory edge can depolarize below its reversal. Excludes leak, intrinsic ionic currents, and externally applied adapter current; not a full substep peak trace.'),indent=2))
  for col,idx in enumerate(act):
   pre,_,coeff=next(row for row in c.outputs if row[1]==idx);name=c.names[pre];side=0 if name.startswith('L ') else 1;incoming=c.edges[:,1].astype(int)==pre;G=(activity[:,incoming]*c.edges[incoming,2]).sum(1);pf=[k for k,e in enumerate(c.edges) if int(e[1])==pre and ' PF ' in c.names[int(e[0])] and not c.names[int(e[0])].endswith(' IN')];assert len(pf)==1;k=pf[0];e=c.edges[k];g=e[2]*activity[:,k];v=z['neuro'][:,pre];l0,width,kse,kpe,damping,maximum,A,B,C,D,*_=params[idx];fl=np.maximum(0,1-((length[:-1,col]-l0)/width)**2)
   def tension(voltage_mV):
    vm=K.polynomial(coeff,voltage_mV*.001);active=np.maximum(0,B/(1+np.exp(np.clip(C*(A-vm),-700,700)))+D)*fl
    return np.clip((kpe*(length[:-1,col]-l0)+active)/(1+kpe/kse),0,maximum)
   for multiplier in ((1.5,2.,3.) if ' Ank ' in name else (1.25,)):
    dg=(multiplier-1)*g;extra=dg*(e[5]-v);dv=extra/(1+G+dg);newforce=tension(v+dv);oldforce=tension(v);delta_torque=-(newforce-oldforce)[:,None]*jac[:-1,col]
    for phase in ('stance','swing'):
     mask=window&((z['foot_vertical'][:,side]>.05*bw) if phase=='stance' else (z['foot_vertical'][:,side]<.025*bw))
     if not mask.any():continue
     fits.append(dict(motor_neuron=name,phase=phase,multiplier=multiplier,source_PF=c.names[int(e[0])],baseline_gmax_uS=float(e[2]),candidate_gmax_uS=float(e[2]*multiplier),baseline_PF_current_mean_nA=float((g*(e[5]-v))[mask].mean()),added_PF_current_p95_nA=float(np.percentile(extra[mask],95)),baseline_MN_mV=float(v[mask].mean()),estimated_MN_mV=float((v+dv)[mask].mean()),actual_tension_mean_N=float(force[mask,col].mean()),zero_velocity_force_before_N=float(oldforce[mask].mean()),estimated_zero_velocity_force_after_N=float(newforce[mask].mean()),estimated_joint_torque_increment_Nm=delta_torque[mask].mean(0).tolist()))
  (folder/'MOTOR_OUTPUT_CURRENT_FIT.json').write_text(json.dumps(dict(source_sha256=r['source_sha256'],rows=fits,scope='Instantaneous frozen-other-neuron conductance equilibrium increments, using actual fitted PF edges, actual spike-filtered proxy voltages, source output adapters, pooled source Hill law, pre-step geometry, and zero-velocity tension clipping. Includes source passive-force offset and force-length factor. Not a dynamic prediction, trained primitive, applied target force, or allocator. Separate coupled air/ground validation is required.'),indent=2))
