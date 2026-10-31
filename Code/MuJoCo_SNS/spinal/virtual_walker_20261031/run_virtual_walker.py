"""Finite source-circuit MuJoCo trials; no force allocator or MN-current override."""
import argparse,json,time,copy,hashlib,sys
from pathlib import Path
import xml.etree.ElementTree as E
import numpy as np
from import_virtual_walker import build,HERE,safe,mujoco,number
import native_kernel as K
from contact_guides import add_guides,Measurements

def run(newer=False,mode='supported',duration=4.,feedback=True,reference=False,tag='',dt=.002,release_at=1.5,contact_si=False,pf_knee_flex_scale=1.,source=None,ankle_pf_crossed=False,symmetric_rg_ext_tonic=None,sagittal_harness=0.,rail_height_offset=0.,neural_time_scale=1.,guide_clearance=None,yaw_perturb_nm=0.,belt_motor_end_s=10.,v3_output_scale=1.,neural_pulse=None):
 path=build(newer,source);folder=path.parent;meta=json.loads((folder/'IMPORT_AUDIT.json').read_text());bank=meta['muscles'];tree=E.parse(path);root=tree.getroot();root.find('option').set('timestep',str(dt))
 if mode=='air':
  b=root.find("worldbody/body[@name='Root']");b.remove(b.find('freejoint'));p=np.array(list(map(float,b.get('pos').split())));p[1]+=.3;b.set('pos',' '.join(map(str,p)))
  for g in root.findall('.//geom'):g.set('contype','0');g.set('conaffinity','0')
 if rail_height_offset:
  for side in ('Left','Right'):
   b=root.find("worldbody/body[@name='"+side+"_support_rail']");p=np.array(list(map(float,b.get('pos').split())));p[1]+=rail_height_offset;b.set('pos',' '.join(map(str,p)))
 guide_geometry=None if guide_clearance is None else add_guides(root,guide_clearance)
 trial=tag or f'{mode}_{"feedback" if feedback else "deafferented"}_{"reference" if reference else "numba"}_{duration:g}s';trialdir=folder/trial;trialdir.mkdir(exist_ok=True)
 xml=trialdir/'model.xml';tree.write(xml,encoding='utf-8',xml_declaration=True);m=mujoco.MjModel.from_xml_path(str(xml));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
 controller=K.NativeController(meta['source_snapshot'],bank,contact_si,pf_knee_flex_scale,ankle_pf_crossed,symmetric_rg_ext_tonic,neural_time_scale,v3_output_scale);params=K.bank_parameters(bank);nb=len(bank);n=round(duration/dt);nj=len(controller.v)
 source_profile_path=Path(meta['source']).parent/'CANONICAL_WORKING_PROFILE.json';source_profile=json.loads(source_profile_path.read_text()) if source_profile_path.exists() else {};source_guides=bool(root.findall(".//body[@name='Guide_carriage_cap']"))
 if source_profile.get('II_homonymous_relays')==12:assert nj==108 and len(controller.edges)==216,'Canonical108 census mismatch: refuse stale runtime export'
 if source_guides and guide_clearance is not None:raise ValueError('Canonical source already contains guides; do not inject a second apparatus')
 if neural_pulse is not None:
  for name in neural_pulse['targets']:
   assert name in controller.names;controller.test_pulses.append((controller.names.index(name),float(neural_pulse['start_s']),float(neural_pulse['start_s'])+float(neural_pulse['duration_s']),float(neural_pulse['amplitude_nA'])))
 tid=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+safe(b['name'])) for b in bank]);act=[(i,mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,safe(b['name']))) for i,b in enumerate(bank) if b['type']=='LinearHillMuscle']
 joints=['hip_L','knee_L','ankle_L','hip_R','knee_R','ankle_R'];jids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,j) for j in joints];qadr=m.jnt_qposadr[jids]
 initial=np.array([-10.,10.,0.,20.,55.,0.]);rootid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');floorid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Treadmill_Belt')
 measurements=Measurements(m,rootid)
 robot_mass=float(m.body_subtreemass[rootid])
 if guide_clearance is not None:assert not any(c.geom1 in measurements.guides or c.geom2 in measurements.guides for c in d.contact),'Initial guide penetration: reject before running'
 rails=[i for i in range(m.ngeom) if 'support_rail' in (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,i) or '')];motor=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,'motor_Treadmill_Slide')
 physical_names=['foot_L','toe_L','foot_R','toe_R'];footids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,x) for x in physical_names];geomid=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+x) for x in physical_names]
 vm=np.full(nb,-.15);internal=np.zeros(nb);force=np.zeros(nb);ia=np.zeros(nb);ii=np.zeros(nb);ib=np.zeros(nb);lastlength=d.ten_length[tid].copy();velocity=np.zeros(nb)
 arrays=dict(t=np.empty(n),q=np.empty((n,m.nq)),angles=np.empty((n,6)),neuro=np.empty((n,nj)),forces=np.empty((n,nb)),root=np.empty((n,3)),rail_vertical=np.empty(n),foot_vertical=np.empty((n,2)),foot_clearance=np.empty((n,4)),harness_force=np.zeros((n,3)),harness_torque=np.zeros((n,3)))
 arrays.update({k:np.empty((n,)+shape) for k,shape in measurements.shapes.items()});arrays['perturb_torque']=np.zeros((n,3));arrays['contact_sensor_binary']=np.zeros((n,4))
 arrays['neural_stimulus_nA']=np.zeros((n,nj))
 native_values={}
 if mode in ('replay','native_kinematic'):
  native=np.load(Path(meta['source']).parent/'verification.npz');z=native['data'];labels=native['labels'].tolist();nt=z[:,labels.index('Time')];native_values={name:np.interp(np.arange(n)*dt,nt,z[:,i]) for i,name in enumerate(labels) if 'angle' in name or 'MN' in name or 'contacts' in name}
 # Warm compiled functions before timing; reset all mutable state afterward.
 K.neural_steps(controller.v.copy(),controller.m.copy(),controller.h.copy(),controller.newv.copy(),controller.gs.copy(),controller.gi.copy(),controller.params,controller.edges,controller.tonic,controller.applied,.2,1)
 K.hill_step(lastlength,velocity,vm,params,internal.copy(),force.copy(),ia.copy(),ii.copy(),ib.copy(),dt)
 oldneural=K.neural_steps;hill=K.hill_step.py_func if reference else K.hill_step
 if reference:K.neural_steps=K.neural_steps.py_func
 started=time.perf_counter();completed=n;failure=None;contacts={};contact_force=np.empty(6);released=False;root_z0=float(d.xpos[rootid,2]);root_R0=d.xmat[rootid].reshape(3,3).copy();spatial_velocity=np.empty(6);rotation_q=np.empty(4);rotation_v=np.empty(3)
 perturb_start=None;phase_index=controller.names.index('L RG flx');previous_phase_voltage=float(controller.v[phase_index])
 try:
  for step in range(n):
   t=step*dt
   if mode in ('released','drop') and not released and (mode=='released' or t>=release_at):
    for g in rails:m.geom_contype[g]=0;m.geom_conaffinity[g]=0
    released=True
   if mode=='native_kinematic':
    for name in ['hip_L','knee_L','ankle_L','toe_L','hip_R','knee_R','ankle_R','toe_R']:
     j=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,name);d.qpos[m.jnt_qposadr[j]]=native_values[name+'_angle'][step]
    mujoco.mj_forward(m,d)
    contacts={name.replace('_contacts',''):values[step] for name,values in native_values.items() if name.endswith('_contacts')}
   else:mujoco.mj_step1(m,d)
   length=d.ten_length[tid];velocity[:]=(length-lastlength)/dt;lastlength[:]=length
   if mode=='replay':
    for pre,target,coeff in controller.outputs:vm[target]=K.polynomial(coeff,native_values[controller.names[pre]][step])
   else:controller.step(t,ia,ii,ib,contacts,dt,feedback);controller.motor_voltages(vm);arrays['neural_stimulus_nA'][step]=controller.last_test_impulse_nA_s/dt
   hill(length,velocity,vm,params,internal,force,ia,ii,ib,dt)
   for b,a in act:d.ctrl[a]=force[b]
   if motor>=0:d.ctrl[motor]=.1 if .5<=t<belt_motor_end_s else 0.
   d.xfrc_applied[rootid,:]=0.
   if sagittal_harness and mode not in ('air','native_kinematic'):
    mujoco.mj_objectVelocity(m,d,mujoco.mjtObj.mjOBJ_BODY,rootid,spatial_velocity,0);delta=d.xmat[rootid].reshape(3,3)@root_R0.T;mujoco.mju_mat2Quat(rotation_q,delta.ravel());mujoco.mju_quat2Vel(rotation_v,rotation_q,1.)
    d.xfrc_applied[rootid,:]=0.;d.xfrc_applied[rootid,2]=sagittal_harness*(-2000*(d.xpos[rootid,2]-root_z0)-150*spatial_velocity[5]);d.xfrc_applied[rootid,3:5]=sagittal_harness*(-200*rotation_v[:2]-20*spatial_velocity[:2]);arrays['harness_force'][step]=d.xfrc_applied[rootid,:3];arrays['harness_torque'][step]=d.xfrc_applied[rootid,3:]
   if yaw_perturb_nm and perturb_start is None and t>=2. and previous_phase_voltage<=-59<controller.v[phase_index]:perturb_start=t
   previous_phase_voltage=float(controller.v[phase_index])
   if perturb_start is not None and perturb_start<=t<perturb_start+.1:d.xfrc_applied[rootid,4]+=yaw_perturb_nm;arrays['perturb_torque'][step,1]=yaw_perturb_nm
   if mode!='native_kinematic':mujoco.mj_step2(m,d);mujoco.mj_forward(m,d)
   measured=measurements.sample(d)
   for key,value in measured.items():arrays[key][step]=value
   contacts={};railforce=0.;feet=np.zeros(2)
   for c in range(d.ncon):
    con=d.contact[c];g1,g2=int(con.geom1),int(con.geom2);b1,b2=int(m.geom_bodyid[g1]),int(m.geom_bodyid[g2]);mujoco.mj_contactForce(m,d,c,contact_force);worldforce=con.frame.reshape(3,3).T@contact_force[:3]
    if g1 in rails or g2 in rails:railforce+=abs(worldforce[1])
    for i,bid in enumerate(footids):
     if bid in (b1,b2):
      if floorid in (b1,b2):feet[i//2]+=abs(worldforce[1])
   clearance=np.empty(4);belt_top=d.xpos[floorid,1]+m.geom_size[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt'),1]
   for i,g in enumerate(geomid):clearance[i]=d.geom_xpos[g,1]-np.abs(d.geom_xmat[g].reshape(3,3)[1])@m.geom_size[g]-belt_top
   # ContactCount is per other BODY, not one count per manifold point. Sense
   # the original ghost boxes geometrically without adding artificial forces.
   for name in physical_names:
    g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+name+'_contact')
    if g>=0:
     bottom=d.geom_xpos[g,1]-np.abs(d.geom_xmat[g].reshape(3,3)[1])@m.geom_size[g]
     belt_geom=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt')
     inside=abs(d.geom_xpos[g,0]-d.geom_xpos[belt_geom,0])<=m.geom_size[belt_geom,0]+.1 and abs(d.geom_xpos[g,2]-d.geom_xpos[belt_geom,2])<=m.geom_size[belt_geom,2]+.1
     contacts[name+'_contact']=float(bottom<=belt_top and inside)
   arrays['contact_sensor_binary'][step]=[contacts.get(name+'_contact',0.) for name in physical_names]
   railforce=float(measured['rail_force_world'][1]);feet=measured['foot_force_world'][:,1]
   arrays['t'][step]=t;arrays['q'][step]=d.qpos;arrays['angles'][step]=initial+np.degrees(d.qpos[qadr]);arrays['neuro'][step]=controller.v;arrays['forces'][step]=force;arrays['root'][step]=d.xpos[rootid];arrays['rail_vertical'][step]=railforce;arrays['foot_vertical'][step]=feet;arrays['foot_clearance'][step]=clearance
   if not np.isfinite(d.qpos).all() or not np.isfinite(force).all():failure='nonfinite';completed=step+1;break
   if mode!='native_kinematic' and min(arrays['angles'][step,[1,4]])<0:failure='knee hyperextension';completed=step+1;break
   if mode not in ('air','native_kinematic') and d.xpos[rootid,1]<.6:failure='root below0.6m: lost upright support';completed=step+1;break
 finally:K.neural_steps=oldneural
 elapsed=time.perf_counter()-started;arrays={k:v[:completed] for k,v in arrays.items()};np.savez_compressed(trialdir/'trajectory.npz',**arrays,neuron_names=controller.names,muscle_names=[b['name'] for b in bank],joint_names=joints)
 summary=dict(machine='easteregg2',source_sha256=meta['source_sha256'],mode=mode,feedback=feedback,reference=reference,contact_si=contact_si,physics_dt=dt,neural_dt=controller.dt,neurons=nj,nonspiking_edges=len(controller.edges),bank_units=nb,muscle_actuators=len(act),duration_requested_s=duration,duration_completed_s=completed*dt,warmed_wall_s=elapsed,simulated_s_per_wall_s=completed*dt/elapsed,failure=failure,root_y_min_m=float(arrays['root'][:,1].min()),root_y_final_m=float(arrays['root'][-1,1]),root_x_change_m=float(arrays['root'][-1,0]-arrays['root'][0,0]),anatomical_knee_min_deg=arrays['angles'][:,[1,4]].min(0).tolist(),joint_range_deg=np.ptp(arrays['angles'],axis=0).tolist(),mean_rail_support_BW=float(np.mean(arrays['rail_vertical'])/((meta['total_dynamic_mass_kg']-4)*9.81)),mean_foot_vertical_N=arrays['foot_vertical'].mean(0).tolist(),scope='Source-specific graded controller and Hill law. Vortex/MuJoCo contact dynamics not equivalent. No training or walking pass inferred from rhythm/translation alone.')
 if mode=='native_kinematic':
  comparisons={}
  for pre,target,coeff in controller.outputs:
   name=controller.names[pre];measured=native_values[name][:completed]*1000;predicted=arrays['neuro'][:,pre];comparisons[name]=dict(rmse_mV=float(np.sqrt(np.mean((predicted-measured)**2))),range_measured_mV=float(np.ptp(measured)),range_port_mV=float(np.ptp(predicted)))
  summary['native_MN_comparison']=comparisons
 summary['shared_pf_knee_flex_scale']=pf_knee_flex_scale;summary['tuned_edges']=controller.tuned_edges
 summary.update(ankle_pf_crossed=ankle_pf_crossed,ankle_phase_edits=controller.ankle_phase_edits,symmetric_rg_ext_tonic_nA=symmetric_rg_ext_tonic,sagittal_harness_scale=sagittal_harness,harness_force_rms_N=np.sqrt(np.mean(arrays['harness_force']**2,axis=0)).tolist(),harness_torque_rms_Nm=np.sqrt(np.mean(arrays['harness_torque']**2,axis=0)).tolist())
 summary['neural_time_scale']=neural_time_scale
 summary['robot_mass_kg']=robot_mass
 summary['belt_motor_end_s']=belt_motor_end_s
 summary.update(v3_output_scale=v3_output_scale,v3_output_edits=controller.v3_output_edits,rg_final_current_terms=controller.rg_current_terms(),rg_current_term_scope='Native-equivalent nA for Size1 source cells evaluated at final state, not the exact previous substep current; includes individually signed synaptic driving terms. Startup pulses have ended by final state.')
 summary.update(guide_clearance_m=guide_clearance,guide_geometry=guide_geometry,guide_force_rms_N=np.sqrt(np.mean(arrays['guide_force_world']**2,axis=0)).tolist(),guide_torque_rms_Nm=np.sqrt(np.mean(arrays['guide_torque_about_pelvis']**2,axis=0)).tolist(),guide_contact_fraction=float(np.mean(arrays['guide_contact_count']>0)),pelvis_euler_XYZ_range_deg=np.ptp(arrays['pelvis_euler_XYZ_deg'],axis=0).tolist(),pelvis_rotational_energy_max_J=float(np.max(arrays['pelvis_rotational_energy_J'])),yaw_perturb_nm=yaw_perturb_nm,yaw_perturb_start_s=perturb_start,yaw_perturb_duration_s=.1,mean_net_vertical_assistance_BW=float(np.mean(arrays['rail_force_world'][:,1]+arrays['guide_force_world'][:,1]+arrays['harness_force'][:,1])/((meta['total_dynamic_mass_kg']-4)*9.81)),measurement_axes='Euler XYZ decomposition Rz(pitchZ) Ry(yawY) Rx(rollX), world angular velocity; pelvis-only rotational energy using its world-oriented principal inertia. Foot forces include all environmental ground bodies, not just moving belt.')
 belt_q=m.jnt_qposadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'Treadmill_Slide')];belt_dx=float(-(arrays['q'][-1,belt_q]-arrays['q'][0,belt_q]));summary.update(rail_height_offset_m=rail_height_offset,belt_world_x_change_m=belt_dx,root_forward_relative_to_belt_m=summary['root_x_change_m']-belt_dx)
 summary['implementation_sha256']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['native_kernel.py','run_virtual_walker.py','import_virtual_walker.py','contact_guides.py']}
 summary['logging']='Post-step mj_forward updates kinematics and contacts before logging; timestamp t denotes command interval start.'
 summary.update(source_neural_clock_scale=source_profile.get('neural_clock_scale',1.),source_profile=source_profile,source_embedded_guides=source_guides,neural_pulse=neural_pulse,neural_pulse_logged_integral_nA_s=np.sum(arrays['neural_stimulus_nA'],axis=0).tolist())
 summary['neural_pulse_logged_integral_nA_s']=(np.sum(arrays['neural_stimulus_nA'],axis=0)*dt).tolist()
 summary['mean_net_vertical_assistance_BW']=float(np.mean(arrays['rail_force_world'][:,1]+arrays['guide_force_world'][:,1]+arrays['harness_force'][:,1])/(robot_mass*9.81))
 if source_guides:summary['guide_clearance_m']=source_profile.get('guide_clearance_m');summary['guide_geometry']=dict(source_embedded=True,side_foreaft_clearance_at_neutral_m=source_profile.get('guide_clearance_m'),scope='Actual canonical source apparatus, not additional runtime guidance')
 (trialdir/'RESULT.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));return summary
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--newer',action='store_true');p.add_argument('--source');p.add_argument('--mode',choices=['supported','released','drop','air','replay','native_kinematic'],default='supported');p.add_argument('--duration',type=float,default=4.);p.add_argument('--feedback-off',action='store_true');p.add_argument('--reference',action='store_true');p.add_argument('--tag',default='');p.add_argument('--dt',type=float,default=.002);p.add_argument('--release-at',type=float,default=1.5);p.add_argument('--contact-si',action='store_true');p.add_argument('--pf-knee-flex-scale',type=float,default=1.);p.add_argument('--ankle-pf-crossed',action='store_true');p.add_argument('--symmetric-rg-ext-tonic',type=float);p.add_argument('--sagittal-harness',type=float,default=0.);p.add_argument('--rail-height-offset',type=float,default=0.);p.add_argument('--neural-time-scale',type=float,default=1.);p.add_argument('--guide-clearance',type=float);p.add_argument('--yaw-perturb-nm',type=float,default=0.);p.add_argument('--belt-motor-end-s',type=float,default=10.);p.add_argument('--v3-output-scale',type=float,default=1.);p.add_argument('--config');a=p.parse_args()
 cfg={} if a.config is None else json.loads(Path(a.config).read_text());allowed={'contact_si','pf_knee_flex_scale','ankle_pf_crossed','symmetric_rg_ext_tonic','sagittal_harness','rail_height_offset','neural_time_scale','guide_clearance','yaw_perturb_nm','belt_motor_end_s','v3_output_scale'}
 if set(cfg)-allowed:raise ValueError('Unconsumed configuration keys: '+str(set(cfg)-allowed))
 explicit_flags={s.split('=')[0] for s in sys.argv[1:] if s.startswith('--')}
 for key,value in cfg.items():
  if '--'+key.replace('_','-') not in explicit_flags:setattr(a,key,value)
 run(a.newer,a.mode,a.duration,not a.feedback_off,a.reference,a.tag,a.dt,a.release_at,a.contact_si,a.pf_knee_flex_scale,a.source,a.ankle_pf_crossed,a.symmetric_rg_ext_tonic,a.sagittal_harness,a.rail_height_offset,a.neural_time_scale,a.guide_clearance,a.yaw_perturb_nm,a.belt_motor_end_s,a.v3_output_scale)
