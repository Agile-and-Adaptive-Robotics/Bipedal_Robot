"""Explicit in-memory opt-in hooks; original CAD hybrid runner stays unchanged."""
import ast,json,hashlib,argparse,time,sys
from pathlib import Path
import numpy as np
import cad_mixed_pressure_experiment as C
R=C.R;mujoco=C.mujoco;OUT=C.OUT;PAIR=C.PAIR
ACTIVE=None
class Mixed:
 def __init__(self,m,d,dynamic,capacities,mode):
  self.m=m;self.d=d;self.dynamic=dynamic;self.mode=mode;self.caps=capacities.copy();self.names=[m.actuator(i).name for i in range(m.nu)];self.ids=np.array([m.actuator(n).id for n in PAIR]);self.tids=m.actuator_trnid[:,0].copy();self.cfg=json.loads((OUT/'ASSUMED_PAIR_CONFIGURATION.json').read_text());self.bank=C.RoutedBPABank(PAIR,self.cfg,m.opt.timestep);self.offset=np.array([self.cfg['actuators'][n]['constant_remote_path_offset_m'] for n in PAIR]);self.logs=[];self.native_gain=m.actuator_gainprm.copy();self.native_bias=m.actuator_biasprm.copy();self.native_gaintype=m.actuator_gaintype.copy();self.native_biastype=m.actuator_biastype.copy()
  if mode!='disabled':m.actuator_gainprm[:]=0;m.actuator_biasprm[:]=0;m.actuator_gaintype[:]=0;m.actuator_biastype[:]=0
 def reported(self,data):
  if self.mode=='disabled':return data.actuator_force
  tensions=self.caps*data.act
  if self.mode=='mixed':tensions[self.ids]=self.bank.force
  return -tensions
 def apply(self,data):
  if self.mode=='disabled':self.dynamic.apply_force_correction(data);return
  self.dynamic.update(data);tensions=self.caps*data.act
  if self.mode=='mixed':
   L=data.ten_length[self.tids[self.ids]]+self.offset;tensions[self.ids]=self.bank.step(data.ctrl[self.ids],L)
   if np.any(self.bank.flags==3):raise ValueError('Requested neural tension below feasible sleeve domain; no implicit vent or force floor')
  Q=-self.dynamic.length_jacobian[self.tids].T@tensions;data.qfrc_applied[:]+=Q
  self.logs.append(dict(t=float(data.time),pressure=self.bank.pressure.copy(),command=self.bank.command.copy(),force=tensions.copy(),neural_output=data.ctrl.copy(),length=data.ten_length[self.tids].copy(),velocity=self.dynamic.velocities[self.tids].copy(),qfrc=Q.copy(),rig=data.xfrc_applied.copy(),knee=np.array(list(self.dynamic.coordinates(data).values())),flags=self.bank.flags.copy()))
  if max(self.logs[-1]['knee'])>0.:raise ValueError('Knee hyperextension safety stop')

def initialize(model,data,dynamic,capacities,args):
 global ACTIVE
 ACTIVE=Mixed(model,data,dynamic,capacities,MODE)
def source_pose(data):data.qpos[:]=np.load(C.ROBOT/'runs/verified65_exact_fast.npz')['qpos'][1069];data.qvel[:]=0
def force(data):return ACTIVE.reported(data)
def apply(dynamic,data):ACTIVE.apply(data)

def compile_opt_in():
 source=Path(R.__file__).read_text();tree=ast.parse(source);run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run');edits=dict(pose=0,init=0,forces=0,apply=0)
 class Transform(ast.NodeTransformer):
  def visit_Assign(self,node):
   node=self.generic_visit(node)
   if any(isinstance(t,ast.Name) and t.id=='data' for t in node.targets) and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Attribute) and node.value.func.attr=='MjData':edits['pose']+=1;return [node,ast.parse('_source_pose(data)').body[0]]
   if any(isinstance(t,ast.Name) and t.id=='capacities' for t in node.targets):edits['init']+=1;return [node,ast.parse('_initialize(model,data,dynamic,capacities,args)').body[0]]
   return node
  def visit_Attribute(self,node):
   if isinstance(node.value,ast.Name) and node.value.id=='data' and node.attr=='actuator_force':edits['forces']+=1;return ast.copy_location(ast.parse('_force(data)',mode='eval').body,node)
   return self.generic_visit(node)
  def visit_Expr(self,node):
   if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Attribute) and isinstance(node.value.func.value,ast.Name) and node.value.func.value.id=='dynamic' and node.value.func.attr=='apply_force_correction':edits['apply']+=1;return ast.copy_location(ast.parse('_apply(dynamic,data)').body[0],node)
   return self.generic_visit(node)
 modified=Transform().visit(run);assert edits==dict(pose=1,init=1,forces=2,apply=1),edits
 namespace=dict(R.__dict__);namespace.update(_source_pose=source_pose,_initialize=initialize,_force=force,_apply=apply);ast.fix_missing_locations(modified);exec(compile(ast.Module(body=[modified],type_ignores=[]),str(Path(__file__)),'exec'),namespace)
 main=next(n for n in tree.body if isinstance(n,ast.If) and isinstance(n.test,ast.Compare) and '__name__' in ast.unparse(n.test));exec(compile(ast.Module(body=main.body[:-1],type_ignores=[]),'<original CAD CLI>','exec'),namespace)
 return namespace['run'],namespace['parser'],ast.unparse(modified),edits

def parity():
 assert json.loads((OUT/'STATIC_AUDIT.json').read_text())['ready_for_live'];cfg=json.loads((OUT/'ASSUMED_PAIR_CONFIGURATION.json').read_text());m=mujoco.MjModel.from_xml_path(str(C.ROBOT/'ah001_neural.xml'));d=mujoco.MjData(m);dyn=C.DynamicRoutes(m,json.loads((C.ROBOT/'robotbody_routes.json').read_text()));C.pose(m,d,dyn);caps=-m.actuator_gainprm[:,0].copy();d.act[:]=np.linspace(.01,.1,m.nu);dyn.update(d);native=d.qfrc_actuator.copy();d.qfrc_applied[:]=0;dyn.apply_force_correction(d);complete=native+d.qfrc_applied;J=dyn.length_jacobian[m.actuator_trnid[:,0]];external=-J.T@(caps*d.act);Qerror=float(np.max(np.abs(complete-external)));velocity=np.linspace(-.01,.01,m.nv);powererror=float(abs(external@velocity-(-caps*d.act)@(J@velocity)));disabled=Mixed(m,d,dyn,caps,'disabled');native_before=d.actuator_force.copy();same=bool(np.array_equal(disabled.reported(d),native_before) and np.array_equal(m.actuator_gainprm,disabled.native_gain) and np.array_equal(m.actuator_biasprm,disabled.native_bias));out=dict(default_force_and_model_exact=same,full_moving_generalized_force_max_error=Qerror,virtual_power_error_W=powererror,passed=same and Qerror<1e-9 and powererror<1e-10,scope='Actual54 CAD routes at source pose1069; fixed-site nativeforce+complete moving correction equals all-external fullJacobian. Disabled interface preserves model/force arrays. Not a fulltrajectory default repeat.');C.save('MOVING_AND_DEFAULT_PARITY.json',out);return out

def run(mode):
 global MODE
 MODE=mode;assert parity()['passed'];fn,parser,code,edits=compile_opt_in();folder=OUT/mode;folder.mkdir(exist_ok=True);(OUT/'source_snapshot/opt_in_run.py').write_text(code,encoding='utf-8');args=parser.parse_args(['--mode','stand','--duration','2','--support','.9','--support-cap','1.2','--anatomical-knee-tangent','--load-compensation','1','--output',str(folder/'cad')]);failure=None
 if mode=='mixed':
  baseline=json.loads((OUT/'reference/BRIDGE_SUMMARY.json').read_text());assert not baseline['failure'] and baseline['duration_s']>=2.-1e-10 and max(baseline['knee_max_rad'])<=0.,'External reference baseline must pass first'
 sources=[Path(__file__),Path(C.__file__),Path(R.__file__),OUT/'ASSUMED_PAIR_CONFIGURATION.json'];protocol=dict(machine='easteregg2',initial_pose_source=str(C.ROBOT/'runs/verified65_exact_fast.npz'),source_pose_index=1069,initial_velocity_zero=True,initial_pose_choice='Existing geometry-envelope sample with bothknees negative and nearlevelpelvis. Original0.1rad newpose cold-domain failure preserved.',duration_s=2.,same_controller='Original262cell SNS with torqueallocator and explicit MNsynaptic cancellation retained',pair_pressure_initial_kpa=[0,0],mode=mode,other52='Reference tension via30ms activation filter',arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},stop_knee_hyperextension=True,strict_invalid_inverse=True,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},no_physical_gain_optimization=True,static_metadata_correction='Earlier PREREGISTRATION.json initial_knee_flexion_rad0.1 text belonged to failedfirstsetup; successful STATIC_AUDIT.npz and this liveprotocol use sourceframe1069 withoutgeometryresizing.')
 C.save(mode+'/LIVE_PREREGISTRATION.json',protocol)
 try:fn(args)
 except ValueError as e:failure=str(e)
 if (folder/'cad.json').exists():
  core=json.loads((folder/'cad.json').read_text());failure=failure or core.get('failed');core['base_runner_force_model_description']=core.get('force_model');core['force_model']=('2assumed single20mm pressure-BPA routes+52ideal reference-tension routes; complete externalmovingJacobian' if mode=='mixed' else '54reference-tension actuators; complete externalmovingJacobian' if mode=='reference' else core.get('force_model'));core['controller_scope']='Hybrid online torqueallocator and MNbackground-synaptic-current cancellation; not autonomous SNScoordination';(folder/'cad.json').write_text(json.dumps(core,indent=2))
 a=ACTIVE;logs=a.logs
 if logs:
  arrays={k:np.asarray([r[k] for r in logs]) for k in logs[0]};np.savez_compressed(folder/'ACTUATOR_BRIDGE.npz',**arrays);rig=arrays['rig'][:,a.m.body('pelvis').id];weight=a.m.body_mass.sum()*9.81;summary=dict(mode=mode,scope='Mixed2virtualpressureBPAs+52reference tensions; original online allocator/MN synaptic cancellation retained, not autonomous neural coordination',failure=failure,steps=len(logs),duration_s=float(a.d.time),knee_max_rad=arrays['knee'].max(axis=0).tolist(),pressure_peak_kpa=float(arrays['pressure'].max()),mean_upward_rig_bw=float(np.maximum(rig[:,1],0).mean()/weight),mean_downward_rig_bw=float(np.maximum(-rig[:,1],0).mean()/weight),lateral_rig_rms_bw=float(np.sqrt(np.mean(rig[:,2]**2))/weight),mean_abs_rig_torque_Nm=np.abs(rig[:,3:]).mean(axis=0).tolist(),pair_force_range_N=np.column_stack((arrays['force'][:,a.ids].min(axis=0),arrays['force'][:,a.ids].max(axis=0))).tolist(),native_force_disabled=bool(mode=='disabled' or(not np.any(a.m.actuator_gainprm) and not np.any(a.m.actuator_biasprm))),source_hook_edits=edits)
 else:summary=dict(mode=mode,failure=failure,steps=0,scope='No live actuated sample completed')
 C.save(mode+'/BRIDGE_SUMMARY.json',summary);print(json.dumps(summary,indent=2));return summary
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['disabled','reference','mixed']);args=p.parse_args()
 if args.mode:run(args.mode)
 else:print(json.dumps(parity(),indent=2))
