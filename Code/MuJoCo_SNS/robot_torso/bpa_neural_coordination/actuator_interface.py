"""Isolated static-BPA inverse and pressure-state prototype; not identified dynamics.

The characterized surface and Fmax conventions are taken from bpa_muscle.py.
This fixture uses zero Xi0/Xi3 and a specified finite series stiffness. It does
not silently substitute these settings for a calibrated CAD route.
"""
from __future__ import annotations
import math
import numpy as np
from numba import njit

# parameters: rest, kmax_fraction, nonmuscle_length, Fmax_total, stiffness,
# a0, a1, a3, pressure_limit. Lengths m; forces N; pressure kPa.
@njit(cache=True)
def force_at_pressure(length, pressure, p):
    lm=length-p[2]
    rel0=(p[0]-lm)/(p[0]*p[1])
    if rel0 < 0.0:
        return np.nan  # outside the declared fixture's contraction domain
    if rel0 >= 1.0:
        return 0.0
    lo=0.0;hi=lm-p[0]*(1.0-p[1])
    for _ in range(46):
        r=0.5*(lo+hi)
        rel=(p[0]-lm+r)/(p[0]*p[1])
        f=max(0.0,p[3]*(p[5]*(math.exp(-p[6]*rel)-1.0)+(pressure/620.0)*math.exp(-p[7]*rel*rel)))
        if f-p[4]*r > 0.0:lo=r
        else:hi=r
    r=0.5*(lo+hi)
    rel=(p[0]-lm+r)/(p[0]*p[1])
    return max(0.0,p[3]*(p[5]*(math.exp(-p[6]*rel)-1.0)+(pressure/620.0)*math.exp(-p[7]*rel*rel)))

@njit(cache=True)
def inverse_pressure(target,length,p):
    """Analytic inverse of the fixture's static series-equilibrium law.

    At equilibrium r=target/K, so contraction is known at requested force.
    Pressure is affine at fixed contraction. Return command and saturation flag.
    """
    if target<=0.0:return 0.0,0
    lm=length-p[2]
    rel=(p[0]-lm+target/p[4])/(p[0]*p[1])
    if (p[0]-lm)<0.0:return 0.0,3
    if rel>=1.0:return p[8],2
    pressure=620.0*(target/p[3]-p[5]*(math.exp(-p[6]*rel)-1.0))*math.exp(p[7]*rel*rel)
    # Bisection roundoff in a forward-generated force can put the inverse a
    # few ulps beyond the endpoint; do not label that an infeasible request.
    if pressure>p[8]+1e-8:return p[8],2
    if pressure<0.0:return 0.0,1
    return min(pressure,p[8]),0

@njit(cache=True)
def pressure_step(actual,command,dt,tau_fill,tau_vent):
    tau=tau_fill if command>actual else tau_vent
    return command+(actual-command)*math.exp(-dt/tau)

@njit(cache=True)
def neural_step(v,scratch,inputs,gmax,elo,ehi,erev,tf,gm,vrest,bias,icon):
    """Preallocated graded-neuron update from SNS-Toolbox compiled arrays.

    No gated/spiking/electrical features; caller verifies those are absent.
    Old voltages are used throughout the simultaneous explicit update.
    """
    scratch[:]=v
    for i in range(v.size):
        sum_g=0.0;sum_ge=0.0;iapp=0.0
        for j in range(v.size):
            span=ehi[i,j]-elo[i,j]
            if span==0.0:span=1.0
            g=min(gmax[i,j],max(0.0,gmax[i,j]*(scratch[j]-elo[i,j])/span))
            sum_g+=g;sum_ge+=g*erev[i,j]
        for j in range(inputs.size):iapp+=icon[i,j]*inputs[j]
        v[i]=scratch[i]+tf[i]*(-gm[i]*(scratch[i]-vrest[i])+bias[i]+sum_ge-scratch[i]*sum_g+iapp)

def packed_parameters(diameter=20,rest=.3,kmax_fraction=.25,nonmuscle_length=.06,stiffness=2e5):
    # Import the existing model only at setup, not from inside compiled loops.
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from bpa_muscle import FESTO_COEFFS,maxBPAforce
    return np.array([rest,kmax_fraction,nonmuscle_length,maxBPAforce(rest,diameter),stiffness,*FESTO_COEFFS[diameter],620.],dtype=float)
