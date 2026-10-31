"""Isolated, approximate damped-equilibrium muscle/tendon prototype.

Not an implementation of Millard2012 or MuJoCo's fitted muscle curves.
All lengths are SI; tension is positive. Fiber projection is the state.
Thelen-informed exponential-to-linear tendon; source-specific strain at F0.
Constant fiber thickness gives pennation. Fiber force has Gaussian active
FL, exponential passive FL, bounded approximate FV, and linear viscous
damping. A backward-Euler equilibrium root is solved for every fiber state.
No production model hooks or changes to the pinned MuJoCo integrator.
"""
import math
import numpy as np
from numba import njit


@njit(cache=True)
def tendon_force(strain,strain_at_f0):
    toe_fraction=99.*math.exp(3.)/(166.*math.exp(3.)-67.)
    toe_end=toe_fraction*strain_at_f0
    if strain<=0.:return 0.
    if strain<toe_end:return .33*math.expm1(3.*strain/toe_end)/math.expm1(3.)
    return .33+.67*(strain-toe_end)/(strain_at_f0-toe_end)


@njit(cache=True)
def residual(s,old,L,a,l0,slack,height,strain,beta,dt,dynamic):
    fiber=math.sqrt(s*s+height*height)
    cosine=s/fiber
    normalized=fiber/l0
    velocity=(s-old)*cosine/(dt*l0) if dynamic else 0.
    active=math.exp(-(normalized-1.)**2/.5)
    fv=1.+math.tanh(velocity/10.)
    passive=math.expm1(4.*(normalized-1.)/.6)/math.expm1(4.) if normalized>1. else 0.
    tendon=tendon_force((L-s)/slack-1.,strain)
    return (a*active*fv+passive+beta*velocity)*cosine-tendon


@njit(cache=True)
def solve_one(old,L,a,l0,slack,height,strain,beta,dt,dynamic):
    lo=.1*l0;hi=max(lo+1e-9,L)
    flo=residual(lo,old,L,a,l0,slack,height,strain,beta,dt,dynamic)
    fhi=residual(hi,old,L,a,l0,slack,height,strain,beta,dt,dynamic)
    bound=0
    if flo>0.:s=lo;bound=1
    elif fhi<0.:s=hi;bound=2
    else:
        for k in range(44):
            s=.5*(lo+hi)
            f=residual(s,old,L,a,l0,slack,height,strain,beta,dt,dynamic)
            if f>0.:hi=s
            else:lo=s
        s=.5*(lo+hi)
    force=tendon_force((L-s)/slack-1.,strain)
    error=residual(s,old,L,a,l0,slack,height,strain,beta,dt,dynamic)
    return s,force,error,bound


@njit(cache=True)
def step_bank(state,length,activation,parameters,dt,force,error,bounds,dynamic=True):
    """Preallocated in-place bank; parameters columns F0,l0,slack,height,strain,beta."""
    for i in range(len(state)):
        f0,l0,slack,height,strain,beta=parameters[i]
        s,t,e,b=solve_one(state[i],length[i],activation[i],l0,slack,height,strain,beta,dt,dynamic)
        state[i]=s;force[i]=f0*t;error[i]=f0*e;bounds[i]=b


@njit(cache=True)
def rigid_bank(state,length,activation,parameters,dt,force):
    """Same approximate fiber law, rigid biological tendon; runtime comparator."""
    for i in range(len(state)):
        f0,l0,slack,height,strain,beta=parameters[i]
        s=max(length[i]-slack,.1*l0);fiber=math.sqrt(s*s+height*height);cosine=s/fiber
        velocity=(s-state[i])*cosine/(dt*l0);normalized=fiber/l0
        passive=math.expm1(4.*(normalized-1.)/.6)/math.expm1(4.) if normalized>1. else 0.
        force[i]=f0*max(0.,activation[i]*math.exp(-(normalized-1.)**2/.5)*(1.+math.tanh(velocity/10.))+passive+beta*velocity)*cosine
        state[i]=s


class SeriesElasticBank:
    def __init__(self,parameters,dt=.002):
        self.parameters=np.ascontiguousarray(parameters,dtype=np.float64)
        if self.parameters.ndim!=2 or self.parameters.shape[1]!=6:raise ValueError('Six parameter columns required')
        if np.any(self.parameters[:,:3]<=0) or np.any(self.parameters[:,4:]<=0):raise ValueError('Positive F0, lengths, strain, damping required')
        self.dt=float(dt);n=len(parameters)
        self.state=np.zeros(n);self.force=np.zeros(n);self.error=np.zeros(n);self.bounds=np.zeros(n,dtype=np.int64)
    def initialize(self,length,activation):
        step_bank(self.state,length,activation,self.parameters,self.dt,self.force,self.error,self.bounds,False)
    def step(self,length,activation):
        step_bank(self.state,length,activation,self.parameters,self.dt,self.force,self.error,self.bounds)
        return self.force
