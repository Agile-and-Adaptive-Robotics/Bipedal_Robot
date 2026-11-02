"""Source-matched Thelen FV extension and analytic constitutive tangent.

Independent closed-form inversion of OpenSim 4.4 Thelen2003's piecewise
force/velocity equations; source parameters must be supplied explicitly.
https://github.com/opensim-org/opensim-core/blob/4.4/OpenSim/Actuators/Thelen2003Muscle.cpp
This module does not modify the earlier approximate series_elastic.py.
"""
import math
import numpy as np
from numba import njit
from series_elastic import tendon_force

def parameters_from_osim(path,names,beta=0.):
    import xml.etree.ElementTree as ET
    root=ET.parse(path);lookup={e.attrib.get('name'):e for e in root.iter() if e.find('max_isometric_force') is not None};rows=[]
    for name in names:
        e=lookup[name];l0=float(e.findtext('optimal_fiber_length'))
        if float(e.findtext('KshapeActive','.5'))!=.5 or float(e.findtext('KshapePassive','4'))!=4 or float(e.findtext('FmaxMuscleStrain','.6'))!=.6:raise ValueError('Static fiber curves differ from implemented source parameters')
        if float(e.findtext('fv_linear_extrap_threshold','.95'))!=.95 or float(e.findtext('minimum_activation','.01'))!=.01:raise ValueError('Source FV threshold/activation floor mismatch')
        rows.append([float(e.findtext('max_isometric_force')),l0,float(e.findtext('tendon_slack_length')),l0*np.sin(float(e.findtext('pennation_angle_at_optimal'))),float(e.findtext('FmaxTendonStrain')),beta,float(e.findtext('Af')),float(e.findtext('Flen')),float(e.findtext('max_contraction_velocity'))])
    p=np.array(rows)
    for i,name in enumerate(names):
        if name.endswith('_r'):
            j=names.index(name[:-2]+'_l');p[i]=p[j]=.5*(p[i]+p[j])
    return p

@njit(cache=True)
def force_velocity(velocity,activation,af,flen,vmax):
    """Return multiplier and derivative; velocity is fiber lengths/second."""
    a=max(.01,activation);c=.25+.75*a;z=velocity/vmax
    if z<=-c:return 0.,0.
    if z<=0.:
        den=c-z/af
        return (c+z)/den,c*(1.+1./af)/(den*den*vmax)
    k=(2.+2./af)/(flen-1.);edge=.95*flen
    threshold=c*(edge-1.)/(k*(flen-edge))
    if z<=threshold:
        den=c+z*k
        return (c+z*k*flen)/den,c*k*(flen-1.)/(den*den*vmax)
    slope=k*(flen-edge)**2/(c*(flen-1.))
    return edge+(z-threshold)*slope,slope/vmax

@njit(cache=True)
def tendon_slope(strain,strain0):
    if strain<=0:return 0.
    edge=99.*math.exp(3.)/(166.*math.exp(3.)-67.)*strain0
    if strain<edge:return .33*3./edge*math.exp(3.*strain/edge)/math.expm1(3.)
    return .67/(strain0-edge)

@njit(cache=True)
def balance_and_slope(s,old,length,activation,p,dt,dynamic=True):
    f0,l0,slack,height,strain0,beta,af,flen,vmax=p
    a=max(.01,activation);fiber=math.hypot(s,height);co=s/fiber;dc=height*height/(fiber**3);nl=fiber/l0
    vel=(s-old)*co/(dt*l0) if dynamic else 0.
    dv=(co+(s-old)*dc)/(dt*l0) if dynamic else 0.
    fl=math.exp(-(nl-1.)**2/.5);dfl=-4.*(nl-1.)*fl
    fv,dfv=force_velocity(vel,a,af,flen,vmax)
    fp=math.expm1(4.*(nl-1.)/.6)/math.expm1(4.) if nl>1. else 0.
    dfp=(4./.6)*math.exp(4.*(nl-1.)/.6)/math.expm1(4.) if nl>1. else 0.
    strain=(length-s)/slack-1.;ft=tendon_force(strain,strain0);kt=tendon_slope(strain,strain0)/slack
    fiberforce=a*fl*fv+fp+beta*vel
    derivative=co*((a*dfl*fv+dfp)*co/l0+(a*fl*dfv+beta)*dv)+fiberforce*dc+kt
    return fiberforce*co-ft,derivative,ft,kt,vel

@njit(cache=True)
def evaluate_one(old,length,activation,p,dt,dynamic=True):
    low=.1*p[1];high=max(low+1e-9,length)
    a=balance_and_slope(low,old,length,activation,p,dt,dynamic)[0]
    b=balance_and_slope(high,old,length,activation,p,dt,dynamic)[0]
    hit=0
    if a>0:s=low;hit=1
    elif b<0:s=high;hit=2
    else:
        for _ in range(48):
            s=.5*(low+high)
            if balance_and_slope(s,old,length,activation,p,dt,dynamic)[0]>0:high=s
            else:low=s
        s=.5*(low+high)
    residual,ds,ft,kt,velocity=balance_and_slope(s,old,length,activation,p,dt,dynamic)
    tangent=p[0]*kt*(1.-kt/ds) if ds>1e-16 else 0.
    return s,p[0]*ft,p[0]*residual,hit,tangent,velocity

@njit(cache=True)
def evaluate_bank(old,length,activation,parameters,dt,state,force,error,bounds,tangent,dynamic=True):
    for i in range(len(old)):
        s,f,e,b,k,v=evaluate_one(old[i],length[i],activation[i],parameters[i],dt,dynamic)
        state[i]=s;force[i]=f;error[i]=e;bounds[i]=b;tangent[i]=k

class DynamicThelenBank:
    """Nine columns: F0,l0,slack,height,strain,beta,Af,Flen,Vmax.

    Evaluation never advances committed state. A caller may evaluate predictors
    or finite differences repeatedly and then commit once at actual new length.
    """
    def __init__(self,parameters,dt=.002):
        self.parameters=np.ascontiguousarray(parameters,dtype=float);self.dt=float(dt)
        if self.parameters.ndim!=2 or self.parameters.shape[1]!=9:raise ValueError('Nine explicit parameter columns required')
        if np.any(self.parameters[:,:3]<=0) or np.any(self.parameters[:,4]<=0) or np.any(self.parameters[:,5]<0):raise ValueError('Invalid source parameters')
        n=len(parameters);self.state=np.zeros(n);self.proposed=np.zeros(n);self.force=np.zeros(n);self.error=np.zeros(n);self.bounds=np.zeros(n,dtype=np.int64);self.tangent=np.zeros(n)
    def initialize(self,length,activation):
        self.evaluate(length,activation,dynamic=False);self.state[:]=self.proposed
    def evaluate(self,length,activation,dynamic=True):
        evaluate_bank(self.state,length,activation,self.parameters,self.dt,self.proposed,self.force,self.error,self.bounds,self.tangent,dynamic)
        return self.force
    def commit(self,length,activation):
        self.evaluate(length,activation);self.state[:]=self.proposed
        return self.force

@njit(cache=True)
def evaluate_rigid_bank(old,length,activation,parameters,dt,state,force,error,bounds,tangent,dynamic=True):
    for i in range(len(old)):
        p=parameters[i];s=max(length[i]-p[2],.1*p[1]);r,k,ft,kt,v=balance_and_slope(s,old[i],length[i],activation[i],p,dt,dynamic)
        f=p[0]*(r+ft);state[i]=s;force[i]=max(0.,f);error[i]=0.;bounds[i]=int(s==.1*p[1]);tangent[i]=p[0]*(k-kt) if f>0 else 0.

@njit(cache=True)
def fiber_components(state,old,activation,p,dt,components):
    for i in range(len(state)):
        f0,l0,slack,height,strain,beta,af,flen,vmax=p[i];s=state[i];l=math.hypot(s,height);co=s/l;nl=l/l0;vel=(s-old[i])*co/(dt*l0)
        fv,_=force_velocity(vel,activation[i],af,flen,vmax);components[0,i]=f0*max(.01,activation[i])*math.exp(-(nl-1.)**2/.5)*fv*co
        components[1,i]=f0*math.expm1(4.*(nl-1.)/.6)/math.expm1(4.)*co if nl>1 else 0.
        components[2,i]=f0*beta*vel*co

class RigidThelenBank(DynamicThelenBank):
    def initialize(self,length,activation):
        self.state[:]=np.maximum(length-self.parameters[:,2],.1*self.parameters[:,1]);self.evaluate(length,activation,False)
    def evaluate(self,length,activation,dynamic=True):
        evaluate_rigid_bank(self.state,length,activation,self.parameters,self.dt,self.proposed,self.force,self.error,self.bounds,self.tangent,dynamic)
        return self.force

class ConstitutiveBackend:
    """Opt-in split-step force/velocity-tangent adapter for one model/data pair.

    A predicted end-step length supplies one Newton linearization. The fiber
    state is committed only once, using the actual resulting length. Source
    curve kinks/contacts make this an approximation, tested by refinement.
    """
    def __init__(self,model,data,parameters,mode='compliant',predictive=True):
        import mujoco
        self.mujoco=mujoco;self.model=model;self.data=data;self.parameters=np.asarray(parameters);self.dt=float(model.opt.timestep);self.predictive=predictive
        self.bank=(DynamicThelenBank if mode=='compliant' else RigidThelenBank)(parameters,self.dt);self.bank.initialize(data.actuator_length.copy(),data.act.copy())
        self.length=np.zeros(model.nu);self.activation=np.zeros(model.nu);self.predicted_force=np.zeros(model.nu);self.length_tangent=np.zeros(model.nu)
        self.old_state=np.zeros(model.nu);self.components=np.zeros((3,model.nu))
        self.max_residual=0.;self.bound_steps=0;self.max_prediction_error=0.;self.max_prediction_error_F0=0.;self.steps=0
        model.actuator_gaintype[:]=mujoco.mjtGain.mjGAIN_FIXED;model.actuator_gainprm[:]=0.;model.actuator_biastype[:]=mujoco.mjtBias.mjBIAS_AFFINE;model.actuator_biasprm[:]=0.
        self.knees=np.array([model.joint('knee_angle_'+side).qposadr[0] for side in ('r','l')]);self.knee_max=data.qpos[self.knees].copy()
    def step(self):
        m,d,j=self.model,self.data,self.mujoco;j.mj_step1(m,d);self.activation[:]=d.act;self.old_state[:]=self.bank.state
        self.length[:]=d.actuator_length+(self.dt*d.actuator_velocity if self.predictive else 0.)
        self.bank.evaluate(self.length,self.activation);self.predicted_force[:]=self.bank.force;self.length_tangent[:]=self.bank.tangent
        self._check()
        slope=-self.dt*self.length_tangent if self.predictive else np.zeros(m.nu)
        m.actuator_biasprm[:,0]=-self.predicted_force-slope*d.actuator_velocity;m.actuator_biasprm[:,2]=slope
        j.mj_step2(m,d);j.mj_forward(m,d)
        # The previous committed fiber state remains unchanged through both
        # predictor and correction, so a physical step is never counted twice.
        self.bank.commit(d.actuator_length,self.activation);self._check()
        fiber_components(self.bank.state,self.old_state,self.activation,self.parameters,self.dt,self.components)
        defect=abs(self.bank.force-(self.predicted_force+self.length_tangent*(d.actuator_length-self.length)))
        self.max_prediction_error=max(self.max_prediction_error,float(np.max(defect)));self.max_prediction_error_F0=max(self.max_prediction_error_F0,float(np.max(defect/self.parameters[:,0])))
        m.actuator_biasprm[:,0]=-self.bank.force-slope*d.actuator_velocity
        j.mj_forward(m,d);self.steps+=1;self.knee_max=np.maximum(self.knee_max,d.qpos[self.knees])
        if np.any(self.knee_max>0.):raise RuntimeError('Strict zero-hyperextension gate')
        if not np.all(np.isfinite(d.qpos)) or np.max(abs(d.qvel))>300:raise RuntimeError('Nonfinite/velocity guard')
    def _check(self):
        self.max_residual=max(self.max_residual,float(np.max(abs(self.bank.error))));self.bound_steps+=int(np.any(self.bank.bounds))
        if np.any(self.bank.bounds):raise RuntimeError('Fiber bound gate')
        if np.max(abs(self.bank.error))>1e-5:raise RuntimeError('Fiber equilibrium gate')
        if not np.all(np.isfinite(self.bank.force)) or np.any(self.bank.force>20*self.parameters[:,0]):raise RuntimeError('Fiber force guard')
    def diagnostics(self):
        return dict(steps=self.steps,max_equilibrium_residual_N=self.max_residual,bound_evaluations=self.bound_steps,max_linearization_error_N=self.max_prediction_error,max_linearization_error_fraction_F0=self.max_prediction_error_F0,maximum_knee_extension_deg=np.degrees(self.knee_max).tolist(),predictive_velocity_tangent=self.predictive,activation_floor=.01,force_velocity='Source-matched Thelen2003, explicit Af/Flen/Vmax')
