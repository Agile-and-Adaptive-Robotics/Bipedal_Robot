"""Opt-in launcher for the existing spinal families; production sources stay unchanged."""
import os,sys,argparse,ast,re,json,time,hashlib,socket
from pathlib import Path

HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1];SPINAL=HERE/'spinal'
parser=argparse.ArgumentParser(description=__doc__,epilog='Put existing runner.py options after --. Use --eval to avoid its legacy fixed-name PNG.')
parser.add_argument('--family',choices=['s3k','w2lvar','syn6','spiking'],default='s3k')
parser.add_argument('--backend',choices=['native','preallocated','sparse'],default='preallocated')
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--model',type=Path,default=REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml')
args,runner_args=parser.parse_known_args()
if runner_args[:1]==['--']:runner_args=runner_args[1:]
if args.family=='spiking' and args.backend!='native':parser.error('The constructed fast kernels support graded networks only; use --backend native for the hybrid spiking mirror.')
known={n.value for n in ast.walk(ast.parse((SPINAL/'runner.py').read_text(encoding='utf-8'))) if isinstance(n,ast.Constant) and isinstance(n.value,str) and re.fullmatch(r'--[\w-]+',n.value)}
unknown=[x for x in runner_args if x.startswith('--') and x not in known]
if unknown:parser.error('Unknown runner options: '+', '.join(unknown))
output=args.output.resolve();model=args.model.resolve()
if output.suffix!='.npz':parser.error('--output must end in .npz')
if output.exists() or output.with_suffix('.launch.json').exists():parser.error('Choose a new output name; this launcher does not overwrite runs.')
if not model.is_file():parser.error('Model does not exist: '+str(model))
output.parent.mkdir(parents=True,exist_ok=True)
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
os.environ['AARL_NET']='' if args.family=='s3k' else args.family
os.environ['AARL_MODEL']=str(model);os.environ['AARL_NPZ']=str(output)
sys.path[:0]=[str(SPINAL),str(SPINAL/'w2l_cpg')];os.chdir(SPINAL)
import numpy as np
import build_network as BN,runner as R
original=BN.build;capture={}

def build(*a,**kw):
    net=original(*a,**kw)
    capture.update(neurons=int(net.compiled.V.size),inputs=len(net.inputs),actuators=list(a[0]))
    if args.backend!='native':
        from preallocated_exact_eval import PreallocatedExactForward
        from sparse_pairwise_eval import SparsePairwiseForward
        raw=net.compiled;fast=(PreallocatedExactForward if args.backend=='preallocated' else SparsePairwiseForward)(raw)
        forward=fast.forward
        def advance(x=None):
            # Existing family graphs are static. Preserve changes to existing
            # conductances; newly added runtime edges require another adapter.
            fast.gmax[:]=raw.g_max_non[fast.post,fast.pre]
            return forward(x)
        fast.forward=advance;net.compiled=fast
    return net
BN.build=build
sources=[Path(__file__),model,SPINAL/'runner.py',SPINAL/'params.py',*SPINAL.glob('build_network*.py'),
    SPINAL/'w2l_cpg'/'preallocated_exact_eval.py',SPINAL/'w2l_cpg'/'sparse_pairwise_eval.py',*SPINAL.glob('*.json')]
meta=dict(machine=socket.gethostname(),family=args.family,backend=args.backend,model=str(model),output=str(output),runner_arguments=runner_args,
    timestep_s=.002,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources if p.is_file()},
    hash_scope='Source/model plus candidate root-level JSON files; inclusion alone does not prove a JSON file was loaded')
output.with_suffix('.launch.json').write_text(json.dumps(meta,indent=2))
start=time.perf_counter();result=R.main(runner_args)
meta.update(wall_s=time.perf_counter()-start,build=capture,metrics=result)
output.with_suffix('.launch.json').write_text(json.dumps(meta,indent=2,default=float))
print(json.dumps(dict(output=str(output),backend=args.backend,wall_s=meta['wall_s'],build=capture),default=float),flush=True)
