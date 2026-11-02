"""Frozen8s matched CAD extension with predeclared stationarity gates."""
import json,hashlib,sys,argparse,shutil
from pathlib import Path
import numpy as np
import cad_pair_unclamped as U
from compact_trial_traces import CompactTraces
BASE=U.OUT;OUT=BASE/'nominal8';C=U.C;L=U.L;DURATION=8.

class RecoverableTraces(CompactTraces):
 def __exit__(self,kind,value,traceback):
  if kind is not None:
   self.materialize(full=True)
   (self.folder/'READOUT_RECOVERY.json').write_text(json.dumps(dict(error_type=kind.__name__,error=str(value),policy='Full RAM archives persisted on readout error; recover assessment without rerunning physics. Scientificparameters unchanged.'),indent=2))
  return super().__exit__(kind,value,traceback)

def read(p):return json.loads(p.read_text())
def protocol():
 OUT.mkdir(exist_ok=True);expected=read(BASE/'RECRUITMENT_PREREGISTRATION.json')
 for name,digest in expected['source_sha256'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest,'Frozensourcechanged:'+name
 core=read(BASE/'mixed/cad.json')
 for rel,digest in core['source_hashes'].items():assert hashlib.sha256((C.ROBOT.parent/rel).read_bytes()).hexdigest()==digest,'Frozenmodel/controllerchanged:'+rel
 p=dict(machine='easteregg2',duration_s=8.,stationary_window_s=[6.,8.],six_absolute_weights_uS=U.weights(),pair_clamp_current_nA=0.,pair_voltage_target_or_floor=False,other52='Originalhybridreferencecontext;unchanged5%floor',apparatus='Unchangednominal0.9BW/vertical800/lateral250/angular300;nominalcap1.2BW',geometry_and_pressure_assumptions_unchanged=True,baseline_gate=str(BASE/'FRESH_REPEAT_AUDIT.json'),both8s_cases_required_for_joint_promotion=True,pressure_case_diagnostic_if_reference8_fails=True,gates=dict(complete_duration=True,no_runtime_failure=True,knee_max_rad=0.,pressure_max_kpa=620.,minimum_visual_foot_mesh_y_m=-.008,mean_upward_rig_bw_max=1.05,mean_ground_load_bw_min=.05,tail_contact_fraction_above1N_min=.95,tail_COM_abs_slope_xyz_m_s=[.005,.002,.005],tail_knee_abs_slope_deg_s=.5,tail_knee_peak_to_peak_deg=2.,tail_pelvis_rotation_abs_slope_deg_s=.5,tail_pelvis_rotation_peak_to_peak_deg=2.,tail_internal_velocity_rms_rad_s=.05,tail_pressure_abs_slope_kpa_s=10.,tail_pressure_inverse_saturation_fraction_max=.01,tail_pair_mean_force_N_min=50.,pair_POST_abs_nA_max=0.,pair_equation_error_mV_max=1e-12),no_retuning_after_failure=True,no_harness_reduction=True,no_endurance_claim=True,retention='RAMfullgates;selectedcurves+compactwindowsforallcases;fullarchivesonlysteadygatepassers',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(U.__file__),Path(C.__file__),Path(L.__file__),BASE/'ASSUMED_PAIR_CONFIGURATION.json')})
 dest=OUT/'PREREGISTRATION.json'
 if dest.exists():assert read(dest)==p,'Frozen8sprotocolchanged'
 else:dest.write_text(json.dumps(p,indent=2))
 return p

def slope(t,x):return np.linalg.lstsq(np.column_stack((t-t.mean(),np.ones(len(t)))),x,rcond=None)[0][0]

def run(mode):
 p=protocol();case=OUT/mode
 if (case/'ASSESSMENT.json').exists():return read(case/'ASSESSMENT.json')
 case.mkdir(exist_ok=True);gate=U.checks();assert gate['passed'];U.checks=lambda:gate;U.OUT=case;C.OUT=case;L.OUT=case
 # This fixture points to the accepted2s prerequisite, not an additionalrun.
 (case/'reference').mkdir(exist_ok=True);r=read(BASE/'reference/BRIDGE_SUMMARY.json');r.update(reference_fixture_only_not_new_run=True,original_reference_summary=str(BASE/'reference/BRIDGE_SUMMARY.json'),role='Known2sprerequisite;matched8sresultisoutsidethisfixture');(case/'reference/BRIDGE_SUMMARY.json').write_text(json.dumps(r,indent=2))
 oldcompile=L.compile_opt_in
 def compiled():
  fn,parser,code,edits=oldcompile();parse=parser.parse_args
  def longer(args=None,*a,**kw):
   args=list(args);i=args.index('--duration');assert args[i+1]=='2';args[i+1]='8';return parse(args,*a,**kw)
  parser.parse_args=longer;return fn,parser,code,edits
 L.compile_opt_in=compiled;oldapply=L.Mixed.apply;contacts=[];oldcontact=C.R.contact_loads
 def contact_monitor(model,data):
  loads=oldcontact(model,data);contacts.append((float(data.time),loads.copy()));return loads
 C.R.contact_loads=contact_monitor
 def observe(self,data):
  previous=len(self.logs);failure=None
  try:oldapply(self,data)
  except ValueError as error:failure=error
  if len(self.logs)>previous:
   if not hasattr(self,'foot_mesh_vertices'):
    self.foot_mesh_vertices=[]
    for g in range(self.m.ngeom):
     mid=self.m.geom_dataid[g]
     if (self.m.geom(g).name or '').startswith('foot_') and self.m.geom_type[g]==C.mujoco.mjtGeom.mjGEOM_MESH and mid>=0:self.foot_mesh_vertices.append((g,self.m.mesh_vert[self.m.mesh_vertadr[mid]:self.m.mesh_vertadr[mid]+self.m.mesh_vertnum[mid]].copy()))
   row=self.logs[-1];pairQ=-self.dynamic.length_jacobian[self.tids[self.ids]].T@row['force'][self.ids];otherQ=row['qfrc']-pairQ;mesh=min(float((v@data.geom_xmat[g].reshape(3,3)[1]+data.geom_xpos[g,1]).min()) for g,v in self.foot_mesh_vertices);row.update(qpos_full=data.qpos.copy(),qvel_full=data.qvel.copy(),whole_COM=data.subtree_com[self.m.body('pelvis').id].copy(),foot_mesh_min_y_m=np.array(mesh),pair_Q=pairQ,other52_Q=otherQ,pair_mechanical_power=np.array(pairQ@data.qvel),other52_mechanical_power=np.array(otherQ@data.qvel))
  if failure is not None:raise failure
 L.Mixed.apply=observe
 with RecoverableTraces(case) as compact:
  summary=U.run(mode);folder=case/mode
  # Retain duration-template metadata explicitly rather than pretend2swas8s.
  for file in (case/'RECRUITMENT_PREREGISTRATION.json',folder/'LIVE_PREREGISTRATION.json'):
   d=read(file);d['inherited_template_duration_s']=d['duration_s'];d['duration_s']=8.;d['authoritative_extension_protocol']=str(OUT/'PREREGISTRATION.json');file.write_text(json.dumps(d,indent=2))
  with np.load(folder/'ACTUATOR_BRIDGE.npz') as z:a={k:z[k] for k in z.files}
  if contacts:b=dict(time=np.array([r[0] for r in contacts]),contact_loads=np.array([r[1] for r in contacts]))
  else:
   with np.load(folder/'cad.npz') as z:b={k:z[k] for k in z.files}
  with np.load(folder/'PAIR_NEURAL_LOG.npz') as z:n=z['samples'][~z['sample_is_warm']]
  t=a['t'];tail=t>=6.;tail_ok=bool(tail.sum()>2 and summary['duration_s']>=8.-1e-9);pelvis=np.array([C.R.rotation_error(q[3:7]) for q in a['qpos_full']]);ground=b['contact_loads'].sum(axis=(1,2));weight=L.ACTIVE.m.body_mass.sum()*9.81;gt=b['time']>=6.;ids=L.ACTIVE.ids;force=a['force'][:,ids];rig=a['rig'][:,L.ACTIVE.m.body('pelvis').id];knees=np.rad2deg(a['knee']);g=p['gates'];m=dict(duration_s=summary['duration_s'],failure=summary['failure'],mesh_min_m=float(a['foot_mesh_min_y_m'].min()),mean_upward_bw=float(np.maximum(rig[:,1],0).mean()/weight),mean_ground_load_bw=float(ground.mean()/weight),knee_max_rad=a['knee'].max(axis=0).tolist(),pressure_peak_kpa=float(a['pressure'].max()),pair_POST_current_max_nA=summary['actual_POST_current_max_nA'],MN_equation_error_mV=summary['equation_error_max_mV'],final_knee_deg=knees[-1].tolist(),final_pressure_kpa=a['pressure'][-1].tolist(),peak_pair_force_N=force.max(axis=0).tolist(),pair_work_J=float(np.trapz(a['pair_mechanical_power'],t)),other52_work_J=float(np.trapz(a['other52_mechanical_power'],t)),tail_window_available=tail_ok)
  gates=dict(complete=summary['duration_s']>=8.-1e-9,no_failure=summary['failure'] is None,knees=max(m['knee_max_rad'])<=0,pressure=m['pressure_peak_kpa']<=620,mesh=m['mesh_min_m']>=-.008,mean_upward=m['mean_upward_bw']<=1.05,mean_ground_load=m['mean_ground_load_bw']>=.05,pair_POST=m['pair_POST_current_max_nA']==0,MN_equation=m['MN_equation_error_mV']<=1e-12,tail_available=tail_ok)
  if tail_ok:
   m.update(tail_COM_slope_xyz_m_s=np.asarray(slope(t[tail],a['whole_COM'][tail])).tolist(),tail_knee_slope_deg_s=np.asarray(slope(t[tail],knees[tail])).tolist(),tail_knee_ptp_deg=np.ptp(knees[tail],axis=0).tolist(),tail_pelvis_rotation_slope_deg_s=np.rad2deg(slope(t[tail],pelvis[tail])).tolist(),tail_pelvis_rotation_ptp_deg=np.rad2deg(np.ptp(pelvis[tail],axis=0)).tolist(),tail_max_internal_velocity_rms_rad_s=float(np.sqrt(np.mean(a['qvel_full'][tail,6:]**2,axis=0)).max()),tail_pressure_slope_kpa_s=np.asarray(slope(t[tail],a['pressure'][tail])).tolist(),tail_saturation_fraction=float(np.mean(a['flags'][tail]!=0)),tail_contact_fraction=float(np.mean(ground[gt]>1)),tail_pair_mean_force_N=force[tail].mean(axis=0).tolist(),tail_pair_Qnorm_mean_Nm=float(np.linalg.norm(a['pair_Q'][tail,6:],axis=1).mean()),tail_other52_Qnorm_mean_Nm=float(np.linalg.norm(a['other52_Q'][tail,6:],axis=1).mean()))
   gates.update(COM_slope=bool(np.all(np.abs(m['tail_COM_slope_xyz_m_s'])<=np.array([.005,.002,.005]))),knee_slope=max(np.abs(m['tail_knee_slope_deg_s']))<=.5,knee_ptp=max(m['tail_knee_ptp_deg'])<=2,pelvis_slope=max(np.abs(m['tail_pelvis_rotation_slope_deg_s']))<=.5,pelvis_ptp=max(m['tail_pelvis_rotation_ptp_deg'])<=2,internal_velocity=m['tail_max_internal_velocity_rms_rad_s']<=.05,pressure_slope=max(np.abs(m['tail_pressure_slope_kpa_s']))<=10,inverse_saturation=m['tail_saturation_fraction']<=.01,contact=m['tail_contact_fraction']>=.95,pair_load=min(m['tail_pair_mean_force_N'])>=50)
  gates={key:bool(value) for key,value in gates.items()}
  assessment=dict(mode=mode,accepted=all(gates.values()),gates=gates,metrics=m,scope='Raw8snominalextension, steadygateonly;noenduranceorharnessweaningclaim');(case/'ASSESSMENT.json').write_text(json.dumps(assessment,indent=2))
  # Selected full-time curves are small; fullhighdimensionalarchives kept onlyonpass.
  other=np.ones(L.ACTIVE.m.nu,dtype=bool);other[ids]=False
  compact.original_save(case/'REPORT_CURVES.npz',t=t,knee_deg=knees,COM=a['whole_COM'],root_position=a['qpos_full'][:,:3],pair_force_N=force,pressure_kpa=a['pressure'],pressure_command_kpa=a['command'],pair_command=a['neural_output'][:,ids],rig_force_N=rig[:,:3],rig_torque_Nm=rig[:,3:],pair_power_W=a['pair_mechanical_power'],other52_power_W=a['other52_mechanical_power'],pair_Qnorm_Nm=np.linalg.norm(a['pair_Q'][:,6:],axis=1),other52_Qnorm_Nm=np.linalg.norm(a['other52_Q'][:,6:],axis=1),other52_activation_sum=a['neural_output'][:,other].sum(axis=1),other52_floor_fraction=(a['neural_output'][:,other]<=.0500001).mean(axis=1),other52_tension_sum_N=a['force'][:,other].sum(axis=1),mesh_y_m=a['foot_mesh_min_y_m'],ground_t=b['time'],ground_N=ground)
 compact.materialize(full=assessment['accepted']);print(json.dumps(assessment,indent=2));return assessment

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',choices=['reference','mixed']);a=p.parse_args()
 if a.case:run(a.case)
 else:print(json.dumps(protocol(),indent=2))
