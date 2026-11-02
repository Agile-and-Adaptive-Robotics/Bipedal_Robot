"""Explicit alternatives for the rat's source continuous linear-Hill equation.

These are not default W2L/controller changes. Parameters use the same fifteen
columns as native_kernel.bank_parameters. Negative internal tension is retained;
only externally applied force is clipped, as in the native source.
"""
import math
import numpy as np
from numba import njit


@njit(cache=True)
def activation_force(length, voltage, p):
    l0, width = p[0], p[1]
    tl = max(0., 1. - ((length-l0)/width)**2)
    return max(0., p[7]/(1.+math.exp(min(700., p[8]*(p[6]-voltage))))+p[9])*tl


@njit(cache=True)
def linear_coefficients(kse, kpe, damping, dt):
    rate = (kse+kpe)/damping
    z = rate*dt
    decay = math.exp(-z)
    phi = -math.expm1(-z)/rate
    if z < 1e-5:
        psi = dt*dt*(.5-z/6.+z*z/24.)
    else:
        psi = (dt-phi)/rate
    return decay, phi, psi


@njit(cache=True)
def analytic_step(length, velocity, voltage, params, internal, force, dt):
    """Exact scalar tension update with length/velocity/activation held per step."""
    for i in range(len(internal)):
        p=params[i];kse=p[2];kpe=p[3];b=p[4]
        decay,phi,_=linear_coefficients(kse,kpe,b,dt)
        forcing=kse/b*(kpe*(length[i]-p[0])+b*velocity[i]+activation_force(length[i],voltage[i],p))
        internal[i]=decay*internal[i]+phi*forcing
        force[i]=min(p[5],max(0.,internal[i]))


@njit(cache=True)
def implicit_coefficients(length, voltage, params, internal, dt, constant, viscous):
    """T_next=constant+viscous*v_next for a linear length segment.

The native active force-length factor is frozen over this numerical step. Its
continuous equation/parameters are unchanged. Passive length and velocity terms
are integrated analytically; positive velocity response permits implicit body
coupling through an affine tendon-actuator velocity bias.
"""
    for i in range(len(internal)):
        p=params[i];kse=p[2];kpe=p[3];b=p[4]
        decay,phi,psi=linear_coefficients(kse,kpe,b,dt)
        constant[i]=decay*internal[i]+phi*kse/b*(kpe*(length[i]-p[0])+activation_force(length[i],voltage[i],p))
        viscous[i]=kse*phi+kse*kpe/b*psi


@njit(cache=True)
def explicit_euler_step(length,velocity,voltage,params,internal,force,dt):
    """Reference native discretization; retained for convergence comparisons."""
    for i in range(len(internal)):
        p=params[i];kse=p[2];kpe=p[3];b=p[4]
        internal[i]+=dt*kse/b*(kpe*(length[i]-p[0])+b*velocity[i]-(1+kpe/kse)*internal[i]+activation_force(length[i],voltage[i],p))
        force[i]=min(p[5],max(0.,internal[i]))


def source_bank(source):
    """Read actual source SI constitutive fields without reusing stale metadata."""
    import xml.etree.ElementTree as E
    root=E.parse(source).getroot();bank=[]
    for body in root.iter('RigidBody'):
        if body.findtext('Type')!='LinearHillMuscle':continue
        def val(path,default=None):return float(body.findtext(path,default))
        p=[val('LengthTension/RestingLength'),val('LengthTension/Lwidth'),val('Kse'),val('Kpe'),val('B'),val('MaximumTension')]
        p += [val('StimulusTension/'+q)for q in 'ABCD']
        p += [val('LengthTension/PeLength'),val('LengthTension/MinPeLength'),val('IaDischarge','100'),val('IIDischarge','100'),val('IbDischarge','100')]
        bank.append({'name':body.findtext('Name'),'id':body.findtext('ID'),'enabled':body.findtext('Enabled')!='False','params':p})
    result=np.array([b['params']for b in bank])
    assert np.all(result[:,2]>0)and np.all(result[:,3]>=0)and np.all(result[:,4]>0)
    return bank,result
