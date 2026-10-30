"""Meaningful force-pattern checks; synthetic cases are not walker results."""
import numpy as np
from ground_contact_check import check_samples
t=np.arange(0.,16.,.002);x=np.zeros((len(t),21));x[:,0]=t;bw=800.
phase=(t/1.2)%1.;right=phase<.62;left=((phase+.5)%1.)<.62
cases={'alternating':(right,left,True),'synchronous':(right,right,False),'permanent_double_stance':(np.ones(len(t),bool),np.ones(len(t),bool),False),'airborne':(np.zeros(len(t),bool),np.zeros(len(t),bool),False)}
for name,(r,l,expected) in cases.items():
    x[:,7]=r*.4*bw;x[:,8]=l*.4*bw;result=check_samples(x,bw)
    assert result['accepted']==expected,(name,result)
    print(name,'passed expected outcome',expected)
