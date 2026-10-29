"""Pairwise sparse reduction against the previously verified backend."""
import os,sys,json,copy,time
from pathlib import Path
H=Path(__file__).resolve().parent;S=H.parent;sys.path[:0]=[str(H),str(S),str(S/'w2l_cpg')]
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
import bsolve_ik as B
import network as F
from sparse_pairwise_eval import SparsePairwiseForward
from preallocated_exact_eval import PreallocatedExactForward
m,d=B.build_air_model(False);acts=[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(m.nu)]
n=F.build(acts,options=dict(axes=True,biarticular=True))
first=PreallocatedExactForward(copy.deepcopy(n.compiled));second=SparsePairwiseForward(copy.deepcopy(n.compiled))
rng=np.random.default_rng(2810);maximum=0.
currents=rng.uniform(0,1,(1000,len(n.inputs)))
for u in currents:
    first.forward(u);second.forward(u)
    maximum=max(maximum,float(np.max(abs(first.V-second.V))))
assert maximum==0.,maximum
timings={}
for name,c in [('preallocated',first),('sparse_pairwise',second)]:
    start=time.perf_counter()
    for u in currents:c.forward(u)
    timings[name]=time.perf_counter()-start
result=dict(maximum_voltage_difference=maximum,steps=1000,timings_s=timings,
            speed_ratio=timings['preallocated']/timings['sparse_pairwise'],neurons=len(first.V))
(S/'reports_20261028_functional_pf/sparse_pairwise_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
