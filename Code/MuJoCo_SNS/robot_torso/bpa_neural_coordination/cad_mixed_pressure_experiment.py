"""Opt-in hybrid CAD controller: two pressure BPAs,52 reference-tension routes."""
import os,sys,json,hashlib,argparse,time,ast
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=r'D:\Anaconda\envs\myo'
import numpy as np,mujoco
HERE=Path(__file__).resolve().parent;ROBOT=HERE.parent;sys.path.insert(0,str(ROBOT))
import run_neural_robot as R
from cad_interface import RoutedBPABank
from actuator_interface import force_at_pressure,inverse_pressure
from dynamic_routes import DynamicRoutes
OUT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\mixed_pressure_pair');PAIR=['vas_med_r','vas_med_l']

def save(name,d):OUT.mkdir(exist_ok=True);(OUT/name).write_text(json.dumps(d,indent=2),encoding='utf-8')
def pose(model,data,dyn):
 # Actual pose already represented in the source401-pose route envelope.
 data.qpos[:]=np.load(ROBOT/'runs/verified65_exact_fast.npz')['qpos'][1069];data.qvel[:]=0;dyn.update(data)

def configuration():
 proposal=json.loads((OUT.parent/'mount_subset/PAIRED_REMOTE_SLEEVE_MOUNT_PROPOSAL.json').read_text());bounds={r['name']:r for r in json.loads((OUT.parent/'route_envelopes/CAD_ROUTE_ENVELOPE_REPORT.json').read_text())['routes']};m=mujoco.MjModel.from_xml_path(str(ROBOT/'ah001_neural.xml'));rows={r['name']:r for r in proposal['routes']};act={}
 for name in PAIR:
  row=rows[name];rest=row['derived_sleeve_rest_m'];offset=row['straight_mount_span_m'];fixed=bounds[name]['route_max_m']+offset-rest
  act[name]=dict(diameter_mm=20,resting_length_m=rest,maximum_contraction_fraction=.25,nonmuscle_length_m=fixed,series_stiffness_N_m=2e5,parallel_count=1,pressure_max_kpa=620.,tau_fill_s=.04,tau_vent_s=.06,delay_s=.004,force_command_scale_N=float(-m.actuator_gainprm[m.actuator(name).id,0]),constant_remote_path_offset_m=offset,provenance='Explicit virtual engineering assumption: source-derived sleeve rest/remote span; fixed length chosen for0contraction at sampled original route maximum; single20mm,200kN/m,40/60ms,4ms delay not hardware identified.')
 assert all(act[PAIR[0]][k]==act[PAIR[1]][k] for k in ('resting_length_m','maximum_contraction_fraction','series_stiffness_N_m','parallel_count','pressure_max_kpa','tau_fill_s','tau_vent_s','delay_s','constant_remote_path_offset_m'))
 return dict(command_mode='inverse_force',actuators=act,scope='Mixed2assumed pressure-BPAs+52ideal reference-tension actuators; onlineallocator and MNsynaptic cancellation preserved/explicit, not autonomous SNS coordination or calibratedhardware.')

def audit():
 OUT.mkdir(exist_ok=True);cfg=configuration();m=mujoco.MjModel.from_xml_path(str(ROBOT/'ah001_neural.xml'));d=mujoco.MjData(m);dyn=DynamicRoutes(m,json.loads((ROBOT/'robotbody_routes.json').read_text()));pose(m,d,dyn);names=[m.actuator(i).name for i in range(m.nu)];ids=np.array([m.actuator(n).id for n in PAIR]);tids=m.actuator_trnid[:,0];caps=-m.actuator_gainprm[:,0].copy();null=R.anatomical_knee_tangent(m,d);matrix=null.T@dyn.moment_arms[tids,6:].T*caps;weight=m.body_mass.sum()*9.81;target=d.qfrc_bias.copy()
 for side in ('r','l'):
  jac=np.zeros((3,m.nv));ang=np.zeros_like(jac);mujoco.mj_jacBodyCom(m,d,jac,ang,m.body('foot_'+side).id);target-=.05*weight*jac[1]
 R.params.G.update(c1_gain=.15,v3_gain=0.,heel_rge=.15,toe_rge=.05);R.params.TAU['rg_nap_h']=.35
 desired=R.allocate(matrix,null.T@target[6:],np.zeros(m.nu));net=R.bn.build(names,dt=m.opt.timestep,interleg=True);net.compiled=R.ExactFastForward(net.compiled);mn=np.array([net.idx['MN_'+n] for n in names]);ports={n:i for i,n in enumerate(net.inputs)};inputs=net.make_inputs();net.compiled.V[mn]=R.params.E_HI*desired
 for _ in range(100):
  post=R.motor_target_current(net.compiled,mn,desired);inputs[:]=0
  for i,n in enumerate(names):inputs[ports['POST_'+n]]=post[i]
  inputs[ports['POSTURE']]=.15;net.step(inputs)
 u=np.clip(net.compiled.V[mn]/R.params.E_HI,0,1);bank=RoutedBPABank(PAIR,cfg,m.opt.timestep);offset=np.array([cfg['actuators'][n]['constant_remote_path_offset_m'] for n in PAIR]);lengths=d.ten_length[tids[ids]]+offset;rows=[];traces=[]
 for k in range(1000):
  force=bank.step(u[ids],lengths);traces.append(np.r_[k*m.opt.timestep,u[ids],lengths,bank.command,bank.pressure,force,bank.flags])
 for j,n in enumerate(PAIR):
  maximum=force_at_pressure(lengths[j],620.,bank.parameters[j]);requested=u[ids[j]]*caps[ids[j]];rows.append(dict(name=n,neural_command=float(u[ids[j]]),route_length_m=float(lengths[j]-offset[j]),full_ideal_remote_route_length_m=float(lengths[j]),reference_force_N=float(requested),static_capacity_N=float(maximum),actual_final_pressure_kpa=float(bank.pressure[j]),actual_final_tension_N=float(bank.force[j]),inverse_flag=int(bank.flags[j]),target_feasible=bool(bank.flags[j]==0 and np.isfinite(maximum) and maximum>=requested),relative_tracking_error=float(abs(bank.force[j]-requested)/max(requested,1.))))
 base=caps*u;mixed=base.copy();mixed[ids]=bank.force;Q0=-dyn.length_jacobian[tids].T@base;Q1=-dyn.length_jacobian[tids].T@mixed;state=dict(qpos=d.qpos.copy(),qvel=d.qvel.copy());np.savez_compressed(OUT/'STATIC_AUDIT.npz',**state,traces=np.asarray(traces),baseline_tension=base,mixed_tension=mixed,moment_arm=dyn.moment_arms[tids],Q_reference=Q0,Q_mixed=Q1)
 src=[Path(__file__),ROBOT/'run_neural_robot.py',ROBOT/'dynamic_routes.py',ROBOT/'robotbody_routes.json',ROBOT/'ah001_neural.xml',HERE/'cad_interface.py',HERE/'actuator_interface.py'];registration=dict(machine='easteregg2',scope=cfg['scope'],pair=PAIR,controller='Existing hybrid torque allocation+SNS MN-current-clamp; background synaptic cancellation retained. No autonomous neuralcoordination claim.',initial_knee_flexion_rad=.1,nominal_support_bw=.9,initial_pose_only_fixed_during_static_force_response=True,actual_ground_step=False,duration_live_s=2.,comparators=['default_native','external_reference_all54','mixed_pressure_pair'],gates=dict(all_finite=True,no_knee_hyperextension=True,no_domain_abort=True,pressure_max_kpa=620.,native_disabled_for_external=True,no_double_moving_correction=True),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in src},static_before_live=True)
 save('ASSUMED_PAIR_CONFIGURATION.json',cfg);save('PREREGISTRATION.json',registration);result=dict(ready_for_live=all(r['target_feasible'] and r['relative_tracking_error']<.05 for r in rows),routes=rows,neurons=len(net.idx),neural_inputs=len(net.inputs),weight_N=weight,knee_coordinates=dyn.coordinates(d),reference_reduced_residual_Nm=float(np.linalg.norm(null.T@(Q0-target)[6:])),mixed_reduced_residual_Nm=float(np.linalg.norm(null.T@(Q1-target)[6:])),incremental_joint_torque_norm_Nm=float(np.linalg.norm((Q1-Q0)[6:])),initial_qpos_preserved=True,geometry_frame='x anterior,y superior,z right',source_pressure_knee_plot='STATIC_AUDIT.npz',note='Fixed-pose sourceNN mapping and staticload/torque audit; support90%BW and ground10%BW are static assumed load split, not measured live forces.');save('STATIC_AUDIT.json',result);print(json.dumps(result,indent=2));return result

if __name__=='__main__':audit()
