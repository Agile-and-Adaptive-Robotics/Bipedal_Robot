import os
import sys
from pathlib import Path
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import time
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'dirusso_gait2392'))
import controller as C
from controller_numba import CompiledDiRusso

assert np.allclose(C.primitives(.1),[1.,0.,0.,0.,0.],atol=1e-14)
assert np.allclose(C.primitives(0.),[.5,0.,0.,0.,.5],atol=1e-14)
W = np.random.default_rng(42).uniform(-.2,.4,(18,5))
gains = {k:(v if np.isscalar(v) else v*.2) for k,v in C.default_gains().items()}
a = C.DiRussoController(W,.002,omega=2*np.pi,gamma=1.,gains=gains)
b = CompiledDiRusso(W,omega=2*np.pi,gamma=1.,gains=gains)
sens = np.zeros((2,4,9)); strikes = np.zeros(2,dtype=np.bool_)
maxdiff = 0.
for i in range(2000):
    t=i*.002
    sens[:,0] = .15*(1.+np.sin(t*6.28))
    sens[:,1] = .5
    sens[:,2] = .1
    sens[:,3] = .4
    pitch=.01*np.sin(t)
    strikes[:] = False
    if i==300:
        strikes[0]=True
    sa = {s:tuple(sens[j,:3])+ (sens[j,3,0],) for j,s in enumerate('LR')}
    aa = a.step(sa,dict(zip('LR',strikes)),pitch).copy()
    bb = b.step_array(sens,strikes,pitch).copy()
    maxdiff=max(maxdiff,float(np.max(np.abs(aa-bb))))
    assert np.allclose(a.phi,b.phi,rtol=0,atol=1e-13)
assert maxdiff < 1e-12, maxdiff
z = C.DiRussoController(np.zeros((18,5)),.002,omega=2*np.pi,gamma=1.)
zero={s:(np.zeros(9),np.zeros(9),np.zeros(9),0.) for s in 'LR'}
for i in range(500):
    z.step(zero,dict(L=False,R=False),0.)
assert np.allclose(z.phi%1.,[0.,.5],atol=1e-12),z.phi
t0=time.perf_counter()
for _ in range(2000):
    a.step(sa,dict(L=False,R=False),0.)
t1=time.perf_counter()
for _ in range(2000):
    b.step_array(sens,strikes,0.)
t2=time.perf_counter()
result=dict(max_output_difference=maxdiff, python_us=(t1-t0)*500,
            numba_us=(t2-t1)*500, speedup=(t1-t0)/(t2-t1))
import json
(Path(__file__).resolve().parent/'dirusso_neural_benchmark.json').write_text(json.dumps(result,indent=2))
print(result)
