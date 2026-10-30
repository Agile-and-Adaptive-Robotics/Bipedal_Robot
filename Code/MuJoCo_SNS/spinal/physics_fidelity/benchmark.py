"""Run isolated equilibrium, runtime, MuJoCo coupling, and damping checks."""
import os,sys,json,time,math,xml.etree.ElementTree as ET
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import expm
import mujoco
from series_elastic import SeriesElasticBank,step_bank,residual,rigid_bank

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent/'runs/audit_20261030'
OUT.mkdir(parents=True,exist_ok=True)
SOURCE=ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'
tree=ET.parse(SOURCE)
muscles=[e for e in tree.iter() if e.find('max_isometric_force') is not None]
params=[];names=[]
for e in muscles:
    f0=float(e.findtext('max_isometric_force'));l0=float(e.findtext('optimal_fiber_length'));slack=float(e.findtext('tendon_slack_length'));angle=float(e.findtext('pennation_angle_at_optimal'))
    assert float(e.findtext('KshapeActive','.5'))==.5
    assert float(e.findtext('KshapePassive','4'))==4.
    assert float(e.findtext('FmaxMuscleStrain','.6'))==.6
    strain=float(e.findtext('FmaxTendonStrain','.049'))
    params.append([f0,l0,slack,l0*math.sin(angle),strain,.05]);names.append(e.attrib['name'])
params=np.asarray(params);assert len(params)==92,len(params)
report=dict(source_parameters=str(SOURCE),muscles=len(params),mujoco_version=mujoco.__version__,dt_s=.002,
    scope='Thelen-informed tendon/FL curves, approximate FV and added fiber damping; no walker integration or Millard/Hyfydy equivalence.')

# Prescribed length/activation, 92 units, 16 s. JIT warmup is separate.
n=len(params);base=params[:,2]+np.sqrt(params[:,1]**2-params[:,3]**2)+params[:,2]*.049*math.sqrt(.15)
act=np.full(n,.15);length=base.copy();bank=SeriesElasticBank(params)
start=time.perf_counter();bank.initialize(length,act);bank.step(length,act);report['jit_and_initialization_s']=time.perf_counter()-start
history=np.zeros((8000,6));max_error=0.;bound_count=0
start=time.perf_counter()
for k in range(8000):
    t=k*.002;length[:]=base+.04*params[:,1]*math.sin(2*math.pi*t);act[:]=.15+.1*math.sin(2*math.pi*t+.7)
    bank.step(length,act);max_error=max(max_error,float(np.max(abs(bank.error))));bound_count+=int(np.count_nonzero(bank.bounds))
    history[k]=[t,bank.state[0],bank.force[0],length[0],act[0],max_error]
elapsed=time.perf_counter()-start
report['bank_16s_wall_s']=elapsed;report['bank_real_time_factor']=16/elapsed
report['maximum_equilibrium_residual_N']=max_error;report['bound_hits']=bound_count
assert max_error<1e-5 and bound_count==0
rigid_state=base-params[:,2];rigid_force=np.zeros(n)
rigid_bank(rigid_state,base,act,params,.002,rigid_force)
start=time.perf_counter()
for k in range(8000):
    t=k*.002;length[:]=base+.04*params[:,1]*math.sin(2*math.pi*t);act[:]=.15+.1*math.sin(2*math.pi*t+.7)
    rigid_bank(rigid_state,length,act,params,.002,rigid_force)
report['rigid_bank_16s_wall_s']=time.perf_counter()-start
report['muscle_only_compliant_to_rigid_runtime_ratio']=elapsed/report['rigid_bank_16s_wall_s']

# Independent SciPy reference trajectory for one muscle, not only root residual.
p=params[names.index('soleus_r')];b=SeriesElasticBank(p[None,:]);L=np.array([p[2]+math.sqrt(p[1]**2-p[3]**2)+p[2]*.049*math.sqrt(.2)]);a=np.array([.2]);b.initialize(L,a);ref=b.state[0];errors=[]
for k in range(500):
    L[0]=p[2]+math.sqrt(p[1]**2-p[3]**2)+p[2]*.049*math.sqrt(.2)+.002*math.sin(k*.002*2*math.pi)
    a[0]=.2+.1*math.sin(k*.002*2*math.pi)
    old=ref;ref=brentq(lambda s:residual.py_func(s,old,L[0],a[0],*p[1:],.002,True),.1*p[1],L[0],xtol=1e-14)
    b.step(L,a);errors.append(abs(ref-b.state[0]))
report['scipy_reference_max_fiber_projection_error_m']=max(errors);assert max(errors)<1e-10

# Backward-Euler convergence against a refined prescribed-length trajectory.
refined={}
for dt in (.002,.001,.0005):
    bb=SeriesElasticBank(p[None,:],dt);aa=np.array([.2]);ll=np.array([p[2]+math.sqrt(p[1]**2-p[3]**2)+.003])
    bb.initialize(ll,aa);trace=np.zeros(round(2./dt))
    for k in range(len(trace)):
        t=(k+1)*dt;ll[0]=p[2]+math.sqrt(p[1]**2-p[3]**2)+.003+.002*math.sin(2*math.pi*t);aa[0]=.2+.1*math.sin(2*math.pi*t)
        trace[k]=bb.step(ll,aa)[0]
    refined[dt]=trace
coarse_error=float(np.sqrt(np.mean((refined[.002]-refined[.0005][3::4])**2)))
fine_error=float(np.sqrt(np.mean((refined[.001]-refined[.0005][1::2])**2)))
report['time_refinement_force_RMS_error_N']=dict(dt_2ms_vs_halfms=coarse_error,dt_1ms_vs_halfms=fine_error)
assert fine_error<coarse_error

# Static equilibrium tends to rigid-tendon force as compliance tends to zero.
rigid=[]
for strain in (.049,.01,.001,.0001):
    p0=np.array([[1000.,.1,.2,0.,strain,.05]]);bb=SeriesElasticBank(p0);bb.initialize(np.array([.3]),np.array([.3]))
    rigid.append(dict(strain_at_F0=strain,force_N=float(bb.force[0]),rigid_reference_N=300.))
assert abs(rigid[-1]['force_N']-300)<.01
report['rigid_limit']=rigid

# Interpreted outer loop vs JIT outer loop, both use identical compiled scalar solver.
# This is deliberately labeled; it is not a full Python-vs-Numba claim.
timings={}
for label,fn in [('python_outer_compiled_scalar',step_bank.py_func),('numba_outer',step_bank)]:
    bb=SeriesElasticBank(params);bb.initialize(base,np.full(n,.2));aa=np.full(n,.2)
    start=time.perf_counter()
    for k in range(200):fn(bb.state,base,aa,params,.002,bb.force,bb.error,bb.bounds)
    timings[label]=time.perf_counter()-start
report['200_step_outer_loop_benchmark_s']=timings

# Damping verification: same spring/mass at 2 ms; compare exact linear solution.
damping=[]
for c in (0.,2.,10.):
    xml=f'<mujoco><option timestep="0.002" integrator="implicitfast" gravity="0 0 0"/><worldbody><body><joint name="slide" type="slide" axis="1 0 0" stiffness="100" damping="{c}"/><geom type="sphere" size=".05" mass="1" contype="0" conaffinity="0"/></body></worldbody></mujoco>'
    m=mujoco.MjModel.from_xml_string(xml);d=mujoco.MjData(m);d.qpos[0]=.02;E=[];error=[]
    A=np.array([[0.,1.],[-100.,-c]])
    for k in range(1000):
        mujoco.mj_step(m,d);E.append(.5*d.qvel[0]**2+50*d.qpos[0]**2)
        exact=expm(A*((k+1)*.002))@np.array([.02,0.]);error.append(abs(d.qpos[0]-exact[0]))
    damping.append(dict(damping_Ns_per_m=c,max_position_error_m=max(error),final_to_initial_energy=E[-1]/.02,max_energy_increase=float(np.max(np.diff(E)))))
assert damping[-1]['final_to_initial_energy']<1e-7
assert damping[-1]['max_position_error_m']<.0002
report['joint_damping']=damping

# Actual MuJoCo spatial-tendon coupling (external muscle state, explicit split).
xml='<mujoco><option timestep="0.002" integrator="implicitfast" gravity="0 0 0"/><worldbody><site name="anchor" pos="0 0 0"/><body pos=".3 0 0"><joint type="slide" axis="1 0 0" stiffness="2000" damping="10"/><geom type="sphere" size=".02" mass="5" contype="0" conaffinity="0"/><site name="tip"/></body></worldbody><tendon><spatial name="mt"><site site="anchor"/><site site="tip"/></spatial></tendon><actuator><general tendon="mt" gainprm="1"/></actuator></mujoco>'
m=mujoco.MjModel.from_xml_string(xml);d=mujoco.MjData(m);mujoco.mj_forward(m,d)
bb=SeriesElasticBank(np.array([[200.,.1,.2,0.,.049,.05]]));aa=np.array([.15]);ll=np.array([d.ten_length[0]]);bb.initialize(ll,aa)
coupled=np.zeros((2000,5));coupled_bounds=0;start=time.perf_counter()
for k in range(2000):
    mujoco.mj_step1(m,d);ll[0]=d.ten_length[0];aa[0]=.15+.1*math.sin(2*math.pi*k*.002)
    bb.step(ll,aa);coupled_bounds+=int(np.count_nonzero(bb.bounds));d.ctrl[0]=-bb.force[0];mujoco.mj_step2(m,d)
    coupled[k]=[d.time,d.qpos[0],d.qvel[0],bb.force[0],bb.error[0]]
report['mujoco_coupling']=dict(simulated_s=4.,wall_s=time.perf_counter()-start,finite=bool(np.all(np.isfinite(coupled))),max_equilibrium_residual_N=float(np.max(abs(coupled[:,4]))),bound_hits=coupled_bounds,max_displacement_m=float(np.max(abs(coupled[:,1]))))
assert report['mujoco_coupling']['finite'] and report['mujoco_coupling']['max_equilibrium_residual_N']<1e-5
assert coupled_bounds==0
np.savez_compressed(OUT/'prototype_traces.npz',bank=history,coupled=coupled)
(OUT/'RESULTS.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
