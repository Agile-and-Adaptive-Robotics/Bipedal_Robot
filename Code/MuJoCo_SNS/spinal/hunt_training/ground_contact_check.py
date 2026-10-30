"""Independently verify load-bearing contact alternation from measured forces."""
import numpy as np
from scipy.ndimage import median_filter

def check(path,bw):
    return check_samples(np.load(path)['samples'],bw)

def check_samples(samples,bw):
    x=samples;keep=(x[:,0]>=5)&(x[:,0]<13);x=x[keep]
    if not len(x):return dict(accepted=False,reason='No scored force window')
    t=x[:,0];onsets=[];duties=[]
    for side in (0,1):
        f=median_filter(x[:,7+side]/bw,size=21);state=False;contact=np.zeros(len(t),bool)
        for i,value in enumerate(f):
            if value>.05:state=True
            elif value<.025:state=False
            contact[i]=state
        boundaries=np.r_[0,np.flatnonzero(contact[1:]!=contact[:-1])+1,len(t)];events=[]
        for k in range(1,len(boundaries)-1):
            begin,end=boundaries[k:k+2];before=boundaries[k-1]
            if contact[begin] and (begin-before)*.002>=.08 and (end-begin)*.002>=.1:events.append(float(t[begin]))
        onsets.append(np.asarray(events));duties.append(float(contact.mean()))
    phases=[];coverage=[]
    for a,b in ((0,1),(1,0)):
        valid=0
        for begin,end in zip(onsets[a][:-1],onsets[a][1:]):
            events=onsets[b][(onsets[b]>begin)&(onsets[b]<end)]
            if len(events)==1:phases.append(float((events[0]-begin)/(end-begin)));valid+=1
        coverage.append(valid/max(1,len(onsets[a])-1))
    error=float(np.mean(abs(np.asarray(phases)-.5))) if phases else .5
    load=np.mean(x[:,7:9],axis=0)/bw;loaded=bool(load.sum()>.2 and np.all(np.max(x[:,7:9],axis=0)/bw>.1))
    accepted=bool(loaded and min(map(len,onsets))>=3 and all(.35<d<.85 for d in duties) and min(coverage)>=.8 and error<.15)
    return dict(accepted=accepted,onsets_s=[a.tolist() for a in onsets],cycles=[len(a) for a in onsets],duty=duties,
        phase_error=error,coverage=coverage,mean_foot_bw=load.tolist(),loaded=loaded,
        mean_harness_vertical_bw=float(np.mean(x[:,1]+x[:,2])/bw),
        scope='Once-counted vertical GRF; 5-13 s window, 0.05/0.025 BW hysteresis, >=80 ms preceding swing and >=100 ms stance. Requires separate movement and frontal gates.')
