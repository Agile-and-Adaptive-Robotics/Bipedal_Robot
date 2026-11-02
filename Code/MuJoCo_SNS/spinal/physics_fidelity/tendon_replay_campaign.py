"""Restartable Gait2392 force-backend replay; preserves one frozen control input.

This is a matched excitation-replay mechanics experiment, not a live SNS fit.
The external rigid/compliant comparators share all fiber laws and physical
parameters; only biological tendon compliance differs. Native-surrogate
replay is retained separately. No production files are modified.
"""
import os,sys,json,time,hashlib,argparse,xml.etree.ElementTree as ET
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent;S=HERE.parent
sys.path[:0]=[str(S),str(S/'functional_gait2392')]
import numpy as np,mujoco
import runner as R
from mechanics import CoupledArms,JOINTS
from series_elastic import SeriesElasticBank,rigid_bank

SOURCE=S/'hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/contact_feasible_start_v1/ground_guide'
OSIM=S.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'
DEFAULT_OUT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\physics')
INTEGRATOR='implicitfast'
PHYSICS_DT=.002
DYNAMIC_SOURCE=False

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source_parameters(names,beta):
 tree=ET.parse(OSIM);lookup={e.attrib.get('name'):e for e in tree.iter() if e.find('max_isometric_force') is not None};p=[]
 for name in names:
  e=lookup[name];l0=float(e.findtext('optimal_fiber_length'));p.append([float(e.findtext('max_isometric_force')),l0,float(e.findtext('tendon_slack_length')),l0*np.sin(float(e.findtext('pennation_angle_at_optimal'))),float(e.findtext('FmaxTendonStrain','.049')),beta])
 p=np.array(p)
 for i,n in enumerate(names):
  if n.endswith('_r'):
   j=names.index(n[:-2]+'_l');p[i]=p[j]=.5*(p[i]+p[j])
 return p

def build(q0,config,damping):
 raw=mujoco.MjModel.from_xml_path(str(R.MODEL));d=mujoco.MjData(raw);d.qpos[:]=q0;mujoco.mj_forward(raw,d)
 m=R.apply_harness(raw,d,kxy=1500.)
 m.opt.timestep=PHYSICS_DT
 if INTEGRATOR=='euler':m.opt.integrator=mujoco.mjtIntegrator.mjINT_EULER
 profile=json.loads(Path(config['force_law_parameters_file']).read_text())
 for row in profile['muscles'].values():
  for name in row['actuator_names']:
   a=m.actuator(name).id
   for k in ('gainprm','biasprm','dynprm','lengthrange'):getattr(m,'actuator_'+k)[a]=row[k]
 d=mujoco.MjData(m);d.qpos[:]=q0
 for name in ('pelvis_list','pelvis_rotation','hip_adduction_r','hip_adduction_l','hip_rotation_r','hip_rotation_l','subtalar_angle_r','subtalar_angle_l','lumbar_bending','lumbar_rotation'):
  j=m.joint(name);m.jnt_stiffness[j.id]+=800.;m.dof_damping[j.dofadr[0]]+=40.
 for side in ('r','l'):
  j=m.joint('knee_angle_'+side);m.jnt_limited[j.id]=True;m.jnt_range[j.id]=np.radians([-90.,0.]);m.jnt_margin[j.id]=np.radians(2.);m.jnt_solref[j.id]=[.004,1.];m.jnt_solimp[j.id]=[.999,.9999,.001,.5,2.]
 for name in ('hip_flexion','knee_angle','ankle_angle'):
  a,b=(m.joint(name+'_'+s).dofadr[0] for s in ('r','l'));v=.5*(m.dof_damping[a]+m.dof_damping[b])*damping;m.dof_damping[a]=m.dof_damping[b]=v
 mujoco.mj_forward(m,d)
 for axis in ('tx','ty','tz'):
  j=m.joint('pelvis_'+axis)
  if abs(d.xaxis[j.id,2])>.8:m.jnt_stiffness[j.id]=0.;m.dof_damping[j.dofadr[0]]=0.
 CoupledArms(m).project(d);mujoco.mj_forward(m,d)
 d.xfrc_applied[m.body('pelvis').id,2]=.7*m.body_mass.sum()*9.81
 return m,d

def run(case,duration,damping,beta,source,out):
 cfg=json.loads(source.with_suffix('.json').read_text());z=np.load(source.with_suffix('.npz'),allow_pickle=True)
 q0=z['qfull'][0].copy();m,d=build(q0,cfg,damping);names=[m.actuator(i).name for i in range(m.nu)];p=source_parameters(names,beta)
 dt=float(m.opt.timestep);n=round(duration/dt);t=z['t'];acts=z['ctrl'];activation=np.ascontiguousarray(np.column_stack([np.interp(np.arange(n)*dt,t,acts[:,i]) for i in range(m.nu)]))
 if np.any(activation<0) or np.any(activation>1):raise ValueError('Saved activation outside [0,1]')
 native_reference=None
 if case in ('native_external','native_affine'):native_reference,native_data=build(q0,cfg,damping)
 if case!='native':
  m.actuator_gaintype[:]=mujoco.mjtGain.mjGAIN_FIXED;m.actuator_gainprm[:]=0.;m.actuator_biastype[:]=mujoco.mjtBias.mjBIAS_NONE;m.actuator_biasprm[:]=0.
 if case=='native_affine':m.actuator_biastype[:]=mujoco.mjtBias.mjBIAS_AFFINE
 d.act[:]=activation[0];d.ctrl[:]=activation[0];mujoco.mj_forward(m,d)
 bank=SeriesElasticBank(p,dt);length=np.ascontiguousarray(d.actuator_length.copy());bank.initialize(length,activation[0]);rigidstate=np.maximum(length-p[:,2],.1*p[:,1]);forces=np.zeros(m.nu)
 # Warm compilation outside measured loop and discard these mutable states.
 dummy=rigidstate.copy();rigid_bank(dummy,length,activation[0],p,dt,forces)
 knee=np.array([m.joint('knee_angle_'+s).qposadr[0] for s in ('r','l')]);qlog=np.empty((n,m.nq));flog=np.empty((n,m.nu));alog=activation.copy();el=np.empty(n);bound=np.empty(n,dtype=int);kmax=d.qpos[knee].copy();count=0;reason=None
 feetlog=np.zeros((n,2));pitchlog=np.empty(n);velog=np.empty((n,6));worklog=np.empty(n);cf=np.zeros(6);feet={m.body(b+'_'+s).id:k for k,s in enumerate(('r','l')) for b in ('calcn','toes')};legvel=[m.joint(j+'_'+s).dofadr[0] for s in ('r','l') for j in ('hip_flexion','knee_angle','ankle_angle')]
 initial_force=None;instant_mapping_error=0.;initial_lengths=length.copy();tic=time.perf_counter()
 for k in range(n):
  d.ctrl[:]=activation[k];d.qfrc_applied[:]=0.;mujoco.mj_step1(m,d);alog[k]=d.act
  if case in ('rigid','compliant'):
   length[:]=d.actuator_length
   if case=='compliant':bank.step(length,d.act);forces[:]=bank.force;el[k]=np.max(abs(bank.error));bound[k]=np.count_nonzero(bank.bounds)
   else:rigid_bank(rigidstate,length,d.act,p,dt,forces);el[k]=0.;bound[k]=0
   if not np.all(np.isfinite(forces)) or np.any(forces>20*p[:,0]):reason='nonfinite_or_over_20F0_force_guard';break
   d.qfrc_applied[:]=d.actuator_moment.T@(-forces)
  elif case in ('native_external','native_affine'):
   native_data.qpos[:]=d.qpos;native_data.qvel[:]=d.qvel;native_data.act[:]=d.act;native_data.ctrl[:]=activation[k];mujoco.mj_forward(native_reference,native_data)
   signed=native_data.actuator_force.copy();forces[:]=-signed;injected=d.actuator_moment.T@signed;instant_mapping_error=max(instant_mapping_error,float(np.max(abs(injected-native_data.qfrc_actuator))));el[k]=0.;bound[k]=0
   if case=='native_external':d.qfrc_applied[:]=injected
   else:
    # A local affine actuator reproduces both force and the velocity tangent
    # used by implicitfast. This is an isolated native-law control, not yet
    # a compliant muscle tangent or a production backend.
    eps=1e-6;v=native_data.actuator_velocity.copy();native_data.actuator_velocity[:]=v+eps;mujoco.mj_fwdActuation(native_reference,native_data);plus=native_data.actuator_force.copy();native_data.actuator_velocity[:]=v-eps;mujoco.mj_fwdActuation(native_reference,native_data);slope=(plus-native_data.actuator_force)/(2*eps);native_data.actuator_velocity[:]=v
    m.actuator_biasprm[:,0]=signed-slope*d.actuator_velocity;m.actuator_biasprm[:,2]=slope
  else:mujoco.mj_fwdActuation(m,d);forces[:]=-d.actuator_force;el[k]=0.;bound[k]=0
  worklog[k]=float(np.sum(-forces*d.actuator_velocity))*dt
  mujoco.mj_step2(m,d);mujoco.mj_forward(m,d)
  if case=='native':forces[:]=-d.actuator_force
  if initial_force is None:initial_force=forces.copy()
  qlog[k]=d.qpos;flog[k]=forces;count=k+1;kmax=np.maximum(kmax,d.qpos[knee])
  pitchlog[k]=d.qpos[m.joint('pelvis_tilt').qposadr[0]];velog[k]=d.qvel[legvel]
  for ci,c in enumerate(d.contact):
   g1,g2=int(c.geom1),int(c.geom2);b1,b2=int(m.geom_bodyid[g1]),int(m.geom_bodyid[g2])
   if m.geom_type[g1]!=mujoco.mjtGeom.mjGEOM_PLANE and m.geom_type[g2]!=mujoco.mjtGeom.mjGEOM_PLANE:continue
   fb=b1 if b1 in feet else b2 if b2 in feet else None
   if fb is not None:mujoco.mj_contactForce(m,d,ci,cf);feetlog[k,feet[fb]]+=max(cf[0],0.)
  if np.any(d.qpos[knee]>0):reason='strict_knee_hyperextension';break
  if not np.all(np.isfinite(d.qpos)) or np.max(abs(d.qvel))>300:reason='nonfinite_or_velocity_guard';break
  if d.time<(k+1)*dt-1e-6:reason='physics_time_reset';break
  if bound[k]:reason='fiber_equilibrium_bound_hit';break
  if el[k]>1e-5:reason='fiber_equilibrium_residual';break
 elapsed=time.perf_counter()-tic
 a=dict(case=case,damping_multiplier=damping,fiber_beta=beta,duration_requested_s=duration,duration_completed_s=count*dt,warmed_wall_s=elapsed,simulated_seconds_per_wall_second=count*dt/elapsed,complete=bool(count==n and reason is None),guard=reason,strict_no_knee_hyperextension=bool(np.all(kmax<=0)),maximum_knee_extension_deg=np.degrees(kmax).tolist(),max_force_N=float(np.max(flog[:count])) if count else None,max_equilibrium_residual_N=float(np.max(el[:count])) if count else None,bound_hit_steps=int(np.sum(bound[:count])) if count else None,initial_force_N=None if initial_force is None else initial_force.tolist(),initial_length_m=initial_lengths.tolist(),instant_force_mapping_error=instant_mapping_error,muscle_names=names,physical_parameters=p.tolist(),scope='Frozen saved excitation replay, not live SNS. Same excitation array/pose/rig/timestep across all backends; native activation dynamics retained, initialized from first excitation; 70% constant upward support, original ground contacts, strict knee gate. Source-informed rigid/compliant approximate fiber laws share tied bilateral parameters. No controller refit or walking inference.')
 a.update(integrator=INTEGRATOR,physics_dt_s=dt,body_weight_N=float(m.body_mass.sum()*9.81),peak_abs_pelvis_tilt_deg=float(np.degrees(np.max(abs(pitchlog[:count])))) if count else None,leg_angular_velocity_rms_rad_s=float(np.sqrt(np.mean(velog[:count]**2))) if count else None,peak_foot_force_BW=float(np.max(feetlog[:count])/ (m.body_mass.sum()*9.81)) if count else None,mean_total_foot_force_BW=float(np.mean(feetlog[:count].sum(1))/(m.body_mass.sum()*9.81)) if count else None,muscle_net_mechanical_work_J=float(np.sum(worklog[:count])) if count else None)
 tag=f'{case}_damp{damping:g}_beta{beta:g}';np.savez_compressed(out/(tag+'.npz'),t=np.arange(count)*dt,qpos=qlog[:count],force=flog[:count],activation=alog[:count],equilibrium_residual=el[:count],bound_hits=bound[:count],foot_force_N=feetlog[:count],pelvis_tilt=pitchlog[:count],leg_velocity=velog[:count],muscle_work_J=worklog[:count]);(out/(tag+'.json')).write_text(json.dumps(a,indent=2));return a

def run_dynamic(case,duration,damping,beta,source,out):
 from thelen_dynamic import parameters_from_osim,ConstitutiveBackend
 if case not in ('rigid','compliant'):return run(case,duration,damping,beta,source,out)
 cfg=json.loads(source.with_suffix('.json').read_text());z=np.load(source.with_suffix('.npz'),allow_pickle=True);m,d=build(z['qfull'][0].copy(),cfg,damping)
 names=[m.actuator(i).name for i in range(m.nu)];p=parameters_from_osim(OSIM,names,beta);dt=float(m.opt.timestep);n=round(duration/dt)
 commands=np.column_stack([np.interp(np.arange(n)*dt,z['t'],z['ctrl'][:,i]) for i in range(m.nu)]);d.act[:]=commands[0];d.ctrl[:]=commands[0];mujoco.mj_forward(m,d)
 backend=ConstitutiveBackend(m,d,p,case);backend.bank.evaluate(d.actuator_length,d.act)
 qlog=np.empty((n,m.nq));flog=np.empty((n,m.nu));alog=np.empty((n,m.nu));residual=np.empty(n);count=0;reason=None;tic=time.perf_counter()
 for k in range(n):
  d.ctrl[:]=commands[k];alog[k]=d.act
  try:backend.step()
  except RuntimeError as exc:reason=str(exc);qlog[k]=d.qpos;flog[k]=backend.bank.force;residual[k]=max(abs(backend.bank.error));count=k+1;break
  qlog[k]=d.qpos;flog[k]=backend.bank.force;residual[k]=max(abs(backend.bank.error));count=k+1
 elapsed=time.perf_counter()-tic;diag=backend.diagnostics();r=dict(case=case,damping_multiplier=damping,fiber_beta=beta,duration_requested_s=duration,duration_completed_s=float(d.time),warmed_wall_s=elapsed,complete=bool(count==n and reason is None),guard=reason,maximum_knee_extension_deg=diag['maximum_knee_extension_deg'],diagnostics=diag,muscle_names=names,physical_parameters=p.tolist(),physics_dt_s=dt,scope='Source-matched Thelen FV with pooled bilateral physical parameters, predicted-length force/velocity tangent and actual-length fiber-state commit; fixed saved neural excitation, not live SNS. Original 70% support mechanical fixture.')
 tag=f'{case}_damp{damping:g}_beta{beta:g}';np.savez_compressed(out/(tag+'.npz'),t=(np.arange(count)+1)*dt,qpos=qlog[:count],force=flog[:count],activation=alog[:count],equilibrium_residual=residual[:count]);(out/(tag+'.json')).write_text(json.dumps(r,indent=2));return r

def analyze(folder,repeat):
 summary=json.loads((folder/'SUMMARY.json').read_text());rows=summary['cases'];answer={}
 n=np.load(folder/'native_damp1_beta0.05.npz');control=np.load(folder/'native_damp1_beta0.1.npz')
 answer['native_beta_negative_control_max_difference']={k:float(np.max(abs(n[k]-control[k]))) for k in ('qpos','force','activation')}
 if repeat:
  a=np.load(folder/'compliant_damp1_beta0.05.npz');b=np.load(repeat/'compliant_damp1_beta0.05.npz');answer['fresh_compliant_repeat_max_difference']={k:float(np.max(abs(a[k]-b[k]))) for k in ('qpos','force','activation','equilibrium_residual')}
 answer['bilateral_parameter_max_difference']=0.
 for r in rows:
  p=np.array(r['physical_parameters']);names=r['muscle_names']
  for i,name in enumerate(names):
   if name.endswith('_r'):answer['bilateral_parameter_max_difference']=max(answer['bilateral_parameter_max_difference'],float(np.max(abs(p[i]-p[names.index(name[:-2]+'_l')]))))
 answer['all_cases_complete']=all(r['complete'] for r in rows);answer['all_cases_no_hyperextension']=all(r['strict_no_knee_hyperextension'] for r in rows);answer['max_fiber_residual_N']=max(r['max_equilibrium_residual_N'] for r in rows)
 selected=[r for r in rows if r['damping_multiplier']==1 and r['fiber_beta']==.05];rigid=next(r for r in selected if r['case']=='rigid');comp=next(r for r in selected if r['case']=='compliant')
 answer['compliant_to_rigid_whole_replay_wall_ratio']=comp['warmed_wall_s']/rigid['warmed_wall_s'];answer['runtime_scope']='Full mechanical replay loop, includes force backend/contacts/logging, excludes neural computation, loading, compilation, and NPZ saving.'
 refinement=[folder.parent/'refinement_1ms',folder.parent/'refinement_0p5ms']
 if all((p/'compliant_damp1_beta0.05.npz').exists() for p in refinement):
  model=mujoco.MjModel.from_xml_path(str(R.MODEL));indices=[model.joint(j).qposadr[0] for j in JOINTS];coarse=np.load(folder/'compliant_damp1_beta0.05.npz');medium=np.load(refinement[0]/'compliant_damp1_beta0.05.npz');fine=np.load(refinement[1]/'compliant_damp1_beta0.05.npz');conv={}
  for label,z,stride in [('2ms',coarse,4),('1ms',medium,2)]:
   n=len(fine['qpos'])//stride;dq=z['qpos'][:n,indices]-fine['qpos'][stride-1::stride,indices];df=z['force'][:n]-fine['force'][::stride]
   conv[label]=dict(joint_RMS_difference_deg=float(np.degrees(np.sqrt(np.mean(dq**2)))),max_joint_difference_deg=float(np.degrees(np.max(abs(dq)))),all_muscle_force_RMS_difference_N=float(np.sqrt(np.mean(df**2))))
  answer['full_model_2s_refinement_vs_0p5ms']=conv
 (folder/'VALIDATION.json').write_text(json.dumps(answer,indent=2))
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'Arial','font.size':11})
 fig,axes=plt.subplots(3,1,figsize=(8,7),sharex=True)
 colors=['#4477AA','#EE6677','#228833'];styles=['-','--','-.']
 for r,color,style in zip(selected,colors,styles):
  z=np.load(folder/f"{r['case']}_damp1_beta0.05.npz");t=z['t']+r.get('physics_dt_s',.002) # states follow each command by one physics step
  bw=float(np.mean(z['foot_force_N'].sum(1))/r['mean_total_foot_force_BW']);index=r['muscle_names'].index('soleus_r');axes[0].plot(t,z['force'][:,index],color=color,ls=style,label=r['case']);axes[1].plot(t,z['foot_force_N'].sum(1)/bw,color=color,ls=style);axes[2].plot(t,np.degrees(z['pelvis_tilt']),color=color,ls=style)
 axes[0].set_ylabel('Right soleus force (N)');axes[1].set_ylabel('Foot normal load (BW)');axes[2].set_ylabel('Pelvis tilt (deg)');axes[2].set_xlabel('Physical time (s)');axes[0].legend(ncol=3,frameon=False)
 fig.suptitle('Gait2392: matched frozen neural-command replay\n70% prescribed upward support; external rigid/compliant laws share parameters',fontsize=12)
 for ax in axes:ax.grid(alpha=.2)
 fig.tight_layout();fig.savefig(folder/'PHYSICS_COMPARISON.png',dpi=180);plt.close(fig);print(json.dumps(answer,indent=2))

def main():
 global INTEGRATOR,PHYSICS_DT,DYNAMIC_SOURCE
 parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=DEFAULT_OUT/'tendon_replay_v1');parser.add_argument('--duration',type=float,default=2.);parser.add_argument('--suite',choices=['preflight','screen'],default='preflight');parser.add_argument('--integrator',choices=['implicitfast','euler'],default='implicitfast');parser.add_argument('--case',choices=['native','rigid','compliant']);parser.add_argument('--damping',type=float,default=1.);parser.add_argument('--beta',type=float,default=.05);parser.add_argument('--dt',type=float,default=.002);parser.add_argument('--source-dynamic',action='store_true');parser.add_argument('--analyze',action='store_true');parser.add_argument('--repeat',type=Path);args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True);INTEGRATOR=args.integrator;PHYSICS_DT=args.dt;DYNAMIC_SOURCE=args.source_dynamic
 if args.analyze:return analyze(args.out,args.repeat)
 manifest=dict(source_activation_trace=str(SOURCE.with_suffix('.npz')),source_sha256=sha(SOURCE.with_suffix('.npz')),source_config_sha256=sha(SOURCE.with_suffix('.json')),source_osim_sha256=sha(OSIM),script_sha256=sha(__file__),model_source=str(R.MODEL),source_model_sha256=sha(R.MODEL),axes='Legacy converted Gait2392: x anterior,y left,z superior; native knee flexion negative; explicit adapter retained.',suite=args.suite,duration=args.duration)
 manifest.update(integrator=INTEGRATOR,physics_dt_s=PHYSICS_DT,single_case=args.case,damping=args.damping,beta=args.beta,runner_sha256=sha(S/'runner.py'),series_elastic_sha256=sha(HERE/'series_elastic.py'),source_dynamic=DYNAMIC_SOURCE,thelen_dynamic_sha256=sha(HERE/'thelen_dynamic.py'));mf=args.out/'MANIFEST.json'
 if mf.exists() and json.loads(mf.read_text())!=manifest:raise RuntimeError('Output manifest differs: use a new named run directory')
 mf.write_text(json.dumps(manifest,indent=2));jobs=[('native',1.,.05),('native_external',1.,.05),('native_affine',1.,.05)] if args.suite=='preflight' else [(c,d,b) for c in ['native','rigid','compliant'] for d,b in [(1.,.05),(2.,.05),(1.,.1)]]
 if args.case:jobs=[(args.case,args.damping,args.beta)]
 rows=[]
 for c,d,b in jobs:
  f=args.out/f'{c}_damp{d:g}_beta{b:g}.json'
  r=json.loads(f.read_text()) if f.exists() else (run_dynamic if DYNAMIC_SOURCE else run)(c,args.duration,d,b,SOURCE,args.out);rows.append(r);print(json.dumps({k:r[k] for k in ['case','damping_multiplier','fiber_beta','complete','guard','duration_completed_s','warmed_wall_s','maximum_knee_extension_deg']}),flush=True)
 summary=dict(manifest=manifest,cases=rows)
 if args.suite=='preflight' and not args.case:
  a=np.load(args.out/'native_damp1_beta0.05.npz');b=np.load(args.out/'native_external_damp1_beta0.05.npz');n=min(len(a['qpos']),len(b['qpos']));summary['external_force_injection_max_qpos_error']=float(np.max(abs(a['qpos'][:n]-b['qpos'][:n]))) if n else None
  c=np.load(args.out/'native_affine_damp1_beta0.05.npz');n=min(len(a['qpos']),len(c['qpos']));summary['affine_velocity_tangent_max_qpos_error']=float(np.max(abs(a['qpos'][:n]-c['qpos'][:n]))) if n else None
 (args.out/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
if __name__=='__main__':main()
