"""Regression: shared learning must pool observations, not average fitted weights."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.optimize import lsq_linear
from bilateral import tied_solutions,assert_tied

def main():
    calls=[]
    def solve(*args,**kwargs):calls.append(1);return lsq_linear(*args,**kwargs)
    names=['PRE_r','POST_r','PRE_l','POST_l'];pre=np.array([0,2]);post=np.array([1,3]);lookup={1:np.array([0]),3:np.array([1])}
    weights,_=tied_solutions({1:np.ones((100,1)),3:np.ones((100,1))*2},
        {1:np.ones(100)*2,3:np.ones(100)*8},lookup,names,pre,post,solve,10.)
    assert len(calls)==1 and len(set(weights.values()))==1
    value=next(iter(weights.values()));assert abs(value-3.6)<.01 and abs(value-3.)>.5
    try:assert_tied({'PRE_r -> POST_r':1.,'PRE_l -> POST_l':2.})
    except ValueError:pass
    else:raise AssertionError('Asymmetric conductances accepted')
    print('Pooled weighted optimum, one shared solve, and asymmetric-file rejection passed:',value)
if __name__=='__main__':main()
