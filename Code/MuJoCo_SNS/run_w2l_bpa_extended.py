"""One W2L BPA trial with graded or hybrid spiking neurons; no Optuna study."""
import os,sys,argparse,json,socket,time,hashlib
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent;S=HERE/'spinal';BODY=S/'w2l_mujoco'
sys.path[:0]=[str(BODY),str(S),str(S/'w2l_cpg')]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--neurons',choices=['graded','spiking'],default='graded')
p.add_argument('--backend',choices=['native','preallocated'],default='preallocated')
p.add_argument('--output',type=Path,required=True);p.add_argument('--duration',type=float,default=20.)
p.add_argument('--support',type=float,default=.7);p.add_argument('--drive-ext',type=float,default=3.)
p.add_argument('--drive-flex',type=float,default=4.);p.add_argument('--tau',type=float,default=.25)
p.add_argument('--comm',type=float,default=1.);p.add_argument('--diameter',type=float,default=20.)
p.add_argument('--count',type=int,default=2);a=p.parse_args()
if a.duration<=8.5:p.error('This drop protocol scores air and retained-support windows; use --duration 12 or longer.')
if a.neurons=='spiking' and a.backend!='native':p.error('Use --backend native for the spiking model; the constructed fast kernel is graded only.')
output=a.output.resolve()
if output.exists() or output.with_suffix('.result.json').exists():p.error('Use a new output path')
output.parent.mkdir(parents=True,exist_ok=True)
started=time.perf_counter()
import _curriculum_bpa as protocol
if a.neurons=='spiking':import build_w2l_split_net_spiking as builder
else:import build_w2l_split_net as builder
knobs=dict(protocol.SEED,dur=a.duration,support=a.support,te=a.drive_ext,tf=a.drive_flex,tau=a.tau,comm=a.comm,dia=a.diameter,count=a.count)
timing={'tau':a.tau} if a.neurons=='spiking' else {'tau_rg_nap_h':a.tau}
net=builder.build(comm=a.comm,template_path=str(S/'connectome_templates.json'),**timing)
if a.backend=='preallocated':
    from preallocated_exact_eval import PreallocatedExactForward
    net.compiled=PreallocatedExactForward(net.compiled)
trial_started=time.perf_counter()
metrics=protocol.run_bpa(net,knobs,save_npz=str(output),verbose=False)
trial_seconds=time.perf_counter()-trial_started
body=protocol.mujoco.MjModel.from_xml_path(protocol.MJCF)
pressure_channels=sorted(body.actuator(i).name for i in range(body.nu))
sources=[Path(__file__),Path(protocol.__file__),Path(builder.__file__),Path(protocol.MJCF),
    S/'connectome_templates.json',S/'params.py',S/'build_network.py',
    BODY/'build_w2l_orig_net.py',S/'w2l_cpg'/'build_w2l_net.py',HERE/'bpa_muscle.py',HERE/'bpa_mujoco.py']
report=dict(machine=socket.gethostname(),neurons=a.neurons,backend=a.backend,parameters=knobs,metrics=metrics,
    scope='Existing virtual-walker drop protocol; not the corrected Hunt five-stage curriculum',physics_timestep_s=.001,controller_call_interval_s=.002,
    trial_wall_seconds=trial_seconds,setup_wall_seconds=trial_started-started,pressure_channels=pressure_channels,
    timing_scope='Trial includes first-call JIT if needed; setup is reported separately',
    source_hashes={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources if path.is_file()})
output.with_suffix('.result.json').write_text(json.dumps(report,indent=2,default=float));print(json.dumps(report,default=float))
