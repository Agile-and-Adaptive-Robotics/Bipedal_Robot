"""Independent installed-OpenSim FV/transient and constitutive tangent checks."""
import os,sys,json,time,math,argparse
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
import opensim as osim
from thelen_dynamic import force_velocity,evaluate_one,DynamicThelenBank,RigidThelenBank

def muscle(name='test',angle=0.):
 m=osim.Thelen2003Muscle(name,1000.,.1,.2,angle);m.set_FmaxTendonStrain(.033);m.set_KshapeActive(.5);m.set_KshapePassive(4.);m.set_Af(.3);m.set_Flen(1.8);m.setMaxContractionVelocity(10.);m.set_fv_linear_extrap_threshold(.95)
 return m

def validate(out):
 out.mkdir(parents=True,exist_ok=True);report={'opensim_version':osim.GetVersion()};cases=[]
 for angle in (0.,.3):
  model=osim.Model();model.setGravity(osim.Vec3(0));m=muscle(angle=angle);m.addNewPathPoint('o',model.getGround(),osim.Vec3(0));m.addNewPathPoint('i',model.getGround(),osim.Vec3(.302,0,0));model.addForce(m);model.finalizeConnections();s=model.initSystem()
  for a in (.01,.05,.3,.8,1.):
   for nl in (.8,1.,1.2):
    for v in (-15.,-8.,-3.,-.1,0.,.1,1.,5.,15.):
     length=.1*nl;co=np.sqrt(1-(.1*np.sin(angle)/length)**2);fl=np.exp(-(nl-1)**2/.5);ref=float(m.calcActiveFiberForceAlongTendon(a,length,v*.1));fv,_=force_velocity(v,a,.3,1.8,10.);ours=1000*a*fl*fv*co
     cases.append([angle,a,nl,v,ref,ours,ours-ref])
 arr=np.array(cases);report['force_velocity_grid_cases']=len(cases);report['max_FV_active_force_error_N']=float(np.max(abs(arr[:,-1])));np.savez_compressed(out/'force_velocity_grid.npz',columns=np.array(['pennation','activation','normalized_length','fiber_lengths_per_s','opensim_force_N','force_N','difference_N']),data=arr)
 # Analytic dF/dL must agree with independent root finite differences away
 # from source curve kinks and unilateral tendon slack boundaries.
 tests=[];rng=np.random.default_rng(301031)
 p=np.array([1000.,.1,.2,.025,.033,.05,.3,1.8,10.])
 for _ in range(120):
  a=float(rng.uniform(.03,.9));L=float(rng.uniform(.29,.32));old=float(rng.uniform(.085,.11));s,f,e,b,k,v=evaluate_one(old,L,a,p,.002)
  if b or abs(v)<1e-3 or f<1e-3:continue
  ep=2e-8;plus=evaluate_one(old,L+ep,a,p,.002)[1];minus=evaluate_one(old,L-ep,a,p,.002)[1];fd=(plus-minus)/(2*ep)
  tests.append(dict(length=L,old=old,activation=a,force_N=f,velocity=v,analytic_N_per_m=k,finite_difference_N_per_m=fd,relative_error=abs(k-fd)/max(1,abs(fd)),residual_N=e))
 report['tangent_test_cases']=len(tests);report['max_tangent_relative_error']=max(t['relative_error'] for t in tests);report['max_tangent_residual_N']=max(abs(t['residual_N']) for t in tests);(out/'TANGENT_CHECKS.json').write_text(json.dumps(tests,indent=2))
 velocity_tests=[]
 for cls in (RigidThelenBank,DynamicThelenBank):
  for beta in (0.,.05):
   for _ in range(35):
    p[5]=beta;a=float(rng.uniform(.03,.9));v=float(rng.uniform(-.04,.04));L=float(rng.uniform(.299,.309));bank=cls(p[None,:],.002);bank.initialize(np.array([L]),np.array([a]));old=bank.state.copy();pred=np.array([L+.002*v]);bank.evaluate(pred,np.array([a]));force=bank.force[0];analytic=-.002*bank.tangent[0]
    if force<1e-3 or abs(v)<.001:continue
    eps=1e-5;bank.evaluate(pred+.002*eps,np.array([a]));plus=-bank.force[0];bank.evaluate(pred-.002*eps,np.array([a]));minus=-bank.force[0];fd=(plus-minus)/(2*eps)
    velocity_tests.append(dict(bank=cls.__name__,beta=beta,analytic_Ns_per_m=float(analytic),finite_difference_Ns_per_m=float(fd),relative_error=float(abs(analytic-fd)/max(1,abs(fd))),state_unmodified=bool(np.array_equal(old,bank.state))))
 report['predictive_velocity_tangent_cases']=len(velocity_tests);report['max_predictive_velocity_tangent_relative_error']=max(t['relative_error'] for t in velocity_tests);report['evaluations_never_commit_state']=all(t['state_unmodified'] for t in velocity_tests);(out/'VELOCITY_TANGENT_CHECKS.json').write_text(json.dumps(velocity_tests,indent=2))
 # Prescribed sinusoidal MT length in actual OpenSim, adaptive reference
 # integration. Activation is controlled; prototype sees measured reference
 # activation and MT length at each sample, without fitting to its force.
 transients=[]
 for a,angle in ((.05,0.),(.3,0.),(.3,.3)):
  model=osim.Model();model.setGravity(osim.Vec3(0));body=osim.Body('moving',1.,osim.Vec3(0),osim.Inertia(1.))
  joint=osim.SliderJoint('slide',model.getGround(),osim.Vec3(.302,0,0),osim.Vec3(0),body,osim.Vec3(0),osim.Vec3(0));coordinate=joint.updCoordinate();coordinate.setName('x');coordinate.setPrescribedFunction(osim.Sine(.002,2*math.pi,math.pi/2));coordinate.setDefaultIsPrescribed(True);model.addBody(body);model.addJoint(joint)
  m=muscle(angle=angle);m.setDefaultActivation(a);m.addNewPathPoint('o',model.getGround(),osim.Vec3(0));m.addNewPathPoint('i',body,osim.Vec3(0));model.addForce(m);controller=osim.PrescribedController();controller.addActuator(m);controller.prescribeControlForActuator('test',osim.Constant(a));model.addController(controller);model.finalizeConnections();state=model.initSystem();model.equilibrateMuscles(state)
  length0=float(m.getLength(state));fiber0=float(m.getFiberLength(state));params=np.array([[1000.,.1,.2,.1*np.sin(angle),.033,0.,.3,1.8,10.]])
  manager=osim.Manager(model);manager.setIntegratorAccuracy(1e-9);manager.initialize(state);reference=[];tic=time.perf_counter()
  for step in range(1,501):
   s=manager.integrate(step*.002);model.realizeDynamics(s);reference.append([s.getTime(),m.getLength(s),m.getActivation(s),m.getTendonForce(s),m.getFiberLength(s),m.getFiberVelocity(s)])
  z=np.array(reference);results={}
  for dt,stride in ((.002,1),(.001,None)):
   bank=DynamicThelenBank(params,dt);bank.state[0]=math.sqrt(fiber0**2-params[0,3]**2);trace=[];hits=0
   for step in range(1,round(1/dt)+1):
    t=step*dt;L=.302+.002*math.sin(2*math.pi*t+math.pi/2);bank.commit(np.array([L]),np.array([a]));hits+=int(bank.bounds[0]);trace.append([t,bank.force[0],math.hypot(bank.state[0],params[0,3]),bank.error[0]])
   trace=np.array(trace);aligned=trace if dt==.002 else trace[1::2];error=aligned[:,1]-z[:,3];tag=f'a{a:g}_pen{angle:g}_dt{dt:g}';np.savez_compressed(out/(tag+'.npz'),reference=z,prototype=trace)
   results[str(dt)]=dict(force_RMS_error_N=float(np.sqrt(np.mean(error**2))),force_max_error_N=float(np.max(abs(error))),fiber_RMS_error_m=float(np.sqrt(np.mean((aligned[:,2]-z[:,4])**2))),bounds=hits,max_equilibrium_residual_N=float(np.max(abs(trace[:,3]))))
  transients.append(dict(activation=a,pennation=angle,initial_MT_length=length0,initial_fiber_length=fiber0,reference_wall_s=time.perf_counter()-tic,results=results));print('TRANSIENT',transients[-1],flush=True)
 report['prescribed_motion_transients']=transients;report['scope']='Thelen2003 source-matched force-velocity and passive/tendon/FL parameters; zero extra fiber damping for OpenSim comparisons. Fixed-width pennation near ordinary lengths; not a complete OpenSim activation/constraint replica.'
 (out/'VALIDATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
 if report['max_FV_active_force_error_N']>.02:raise RuntimeError('FV source-equivalence gate')
 if report['max_tangent_relative_error']>1e-4:raise RuntimeError('Constitutive tangent finite-difference gate')
 if report['max_predictive_velocity_tangent_relative_error']>1e-4 or not report['evaluations_never_commit_state']:raise RuntimeError('Predictive velocity tangent or state gate')
 if any(v['bounds'] for r in transients for v in r['results'].values()):raise RuntimeError('Unexpected fiber bound')

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\dynamic_source_v1'));args=parser.parse_args();validate(args.out)
