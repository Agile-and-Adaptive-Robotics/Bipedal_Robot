import os
import sys
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import json
import time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'w2l_cpg'))
import build_w2l_net as B
from exact_fast_eval import ExactFastForward

def main():
    knobs=json.loads((HERE/'trial458_d30.json').read_text())['delivered_gains']
    a=B.build(template_key='bilateralrg2',ib_ports=True,gains=knobs).compiled
    b=B.build(template_key='bilateralrg2',ib_ports=True,gains=knobs).compiled
    fast=ExactFastForward(b)
    u=np.random.default_rng(42).uniform(0,4,a.input_connectivity.shape[1]).tolist()
    worst=0.
    t0=time.perf_counter()
    for _ in range(5000):
        a.forward(u); fast.forward(u)
        worst=max(worst,float(np.max(np.abs(a.V-fast.V))))
    assert worst==0.,worst
    t0=time.perf_counter()
    for _ in range(10000):
        a.forward(u)
    t1=time.perf_counter()
    for _ in range(10000):
        fast.forward(u)
    t2=time.perf_counter()
    result=dict(steps_compared=5000,max_voltage_difference=worst,
                dense_us=(t1-t0)*100,exact_fast_us=(t2-t1)*100,
                speedup=(t1-t0)/(t2-t1))
    (HERE/'benchmark_exact.json').write_text(json.dumps(result,indent=2))
    print(result)

if __name__=='__main__':
    main()
