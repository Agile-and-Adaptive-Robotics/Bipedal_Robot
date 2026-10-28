"""Parity and warmed serial timings; no physics or production source edits."""
import os,sys,json,copy,time,argparse
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE/'w2l_cpg')]
import numpy as np
import params as P
import build_network as B
from exact_fast_eval import ExactFastForward
from preallocated_exact_eval import PreallocatedExactForward
p=argparse.ArgumentParser();p.add_argument('family');p.add_argument('body');a=p.parse_args()
audit=json.loads((HERE/'reports_20261028_family_audit'/(a.family+'_'+a.body+'.json')).read_text())
P.G.update(audit['gains']);P.TAU.update(audit['time_constants'])
if a.family=='s3k':net=B.build(audit['actuators'])
elif a.family=='w2lvar':
    import build_network_w2lvar as V
    net=V.build(audit['actuators'])
else:
    import build_network_syn6 as V
    net=V.build(audit['actuators'])
ref=copy.deepcopy(net.compiled);pre=PreallocatedExactForward(copy.deepcopy(net.compiled))
rng=np.random.default_rng(20261028);inputs=rng.uniform(0,.3,(400,len(net.inputs)))
inputs[:,net.inputs.index('DRIVE')]=3.;inputs[:,net.inputs.index('POSTURE')]=.15
maximum=0.;output_max=0.;gate_max=0.
for x in inputs:
    ref.forward(x.tolist());pre.forward(x)
    maximum=max(maximum,float(np.max(abs(ref.V-pre.V))))
    if ref.outputs.size:output_max=max(output_max,float(np.max(abs(ref.outputs-pre.outputs))))
    if ref.gated:gate_max=max(gate_max,float(np.max(abs(ref.b_gate-pre.b_gate))))
assert maximum==output_max==gate_max==0.,(maximum,output_max,gate_max)
times={}
for label,cls in [('reference',lambda c:c),('exact_fast',ExactFastForward),('preallocated',PreallocatedExactForward)]:
    obj=cls(copy.deepcopy(net.compiled))
    for x in inputs[:10]:obj.forward(x)
    samples=[]
    for _ in range(3):
        start=time.perf_counter()
        for x in inputs:obj.forward(x)
        samples.append((time.perf_counter()-start)/len(inputs))
    times[label]=float(np.median(samples))
result=dict(family=a.family,body=a.body,neurons=ref.V.size,steps=400,
    max_voltage_difference=maximum,max_output_difference=output_max,max_h_gate_difference=gate_max,
    bit_exact=True,input_onehot=bool(pre.input_onehot),output_onehot=bool(pre.output_onehot),
    median_step_s=times,reference_speedup=times['reference']/times['preallocated'],
    previous_fast_speedup=times['exact_fast']/times['preallocated'])
(HERE/'reports_20261028_family_audit'/(a.family+'_'+a.body+'_preallocated.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
