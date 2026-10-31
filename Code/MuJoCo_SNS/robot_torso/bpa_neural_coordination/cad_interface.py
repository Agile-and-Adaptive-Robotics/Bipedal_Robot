"""Opt-in pressure-state plumbing for routed MuJoCo BPA actuators.

This module has no torque allocator and never edits neural state. Route-specific
hardware/assumption parameters must be supplied; no human muscle length or MIF
is silently used as a BPA dimension or measured capacity.
"""
from __future__ import annotations
import numpy as np
from actuator_interface import packed_parameters,force_at_pressure,inverse_pressure,pressure_step

REQUIRED=('diameter_mm','resting_length_m','maximum_contraction_fraction',
          'nonmuscle_length_m','series_stiffness_N_m','parallel_count',
          'pressure_max_kpa','tau_fill_s','tau_vent_s','delay_s','provenance')

def validate_configuration(configuration,names):
    routes=configuration.get('actuators',{})
    problems={}
    for name in names:
        r=routes.get(name,{})
        missing=[k for k in REQUIRED if r.get(k) is None or r.get(k)=='']
        if configuration.get('command_mode','direct_pressure')=='inverse_force' and r.get('force_command_scale_N') is None:
            missing.append('force_command_scale_N')
        if missing:problems[name]=missing
    return problems

class RoutedBPABank:
    """Pressure and BPA force bank; caller owns SNS, support, and physics steps."""
    def __init__(self,names,configuration,dt):
        problems=validate_configuration(configuration,names)
        if problems:raise ValueError('Missing route-specific BPA parameters: '+str(problems))
        self.names=list(names);self.dt=float(dt);self.mode=configuration.get('command_mode','direct_pressure')
        if self.mode not in ('direct_pressure','inverse_force'):raise ValueError('Unknown command mode')
        params=[];fill=[];vent=[];delay=[];scales=[]
        for name in names:
            r=configuration['actuators'][name]
            assert r['diameter_mm'] in (10,20,40)
            assert r['resting_length_m']>0 and 0<r['maximum_contraction_fraction']<1
            assert r['nonmuscle_length_m']>=0 and r['series_stiffness_N_m']>0
            assert int(r['parallel_count'])==r['parallel_count'] and r['parallel_count']>=1
            assert 0<r['pressure_max_kpa']<=620 and r['tau_fill_s']>0 and r['tau_vent_s']>0 and r['delay_s']>=0
            p=packed_parameters(r['diameter_mm'],r['resting_length_m'],r['maximum_contraction_fraction'],r['nonmuscle_length_m'],r['series_stiffness_N_m'])
            p[3]*=r['parallel_count'];p[8]=r['pressure_max_kpa'];params.append(p)
            fill.append(r['tau_fill_s']);vent.append(r['tau_vent_s']);delay.append(round(r['delay_s']/dt));scales.append(r.get('force_command_scale_N',0.))
        self.parameters=np.array(params);self.fill=np.array(fill);self.vent=np.array(vent);self.delay=np.array(delay,dtype=int);self.scales=np.array(scales)
        self.pressure=np.zeros(len(names));self.command=np.zeros(len(names));self.force=np.zeros(len(names));self.flags=np.zeros(len(names),dtype=int)
        self.history=np.zeros((len(names),int(max(delay,default=0))+1));self.step_count=0

    def step(self,normalized_neural_output,lengths):
        u=np.asarray(normalized_neural_output)
        if u.shape!=self.pressure.shape or np.asarray(lengths).shape!=u.shape:raise ValueError('Bank input size mismatch')
        if not np.isfinite(u).all() or not np.isfinite(lengths).all():raise ValueError('Nonfinite bank input')
        for i in range(len(self.names)):
            activation=min(1.,max(0.,float(u[i])))
            if self.mode=='direct_pressure':self.command[i]=activation*self.parameters[i,8];self.flags[i]=0
            else:self.command[i],self.flags[i]=inverse_pressure(activation*self.scales[i],lengths[i],self.parameters[i])
            self.history[i,self.step_count%self.history.shape[1]]=self.command[i]
            delayed=self.history[i,(self.step_count-self.delay[i])%self.history.shape[1]] if self.step_count>=self.delay[i] else 0.
            self.pressure[i]=pressure_step(self.pressure[i],delayed,self.dt,self.fill[i],self.vent[i])
            self.force[i]=force_at_pressure(lengths[i],self.pressure[i],self.parameters[i])
        self.step_count+=1
        if not np.isfinite(self.force).all():raise ValueError('BPA length is outside the declared contraction domain')
        return self.force

    @staticmethod
    def generalized_force(length_jacobian,tensions):
        """Virtual work: tensile force acts to SHORTEN the path: Q=-J_L^T T.

        Use DynamicRoutes.length_jacobian for CAD moving-point terms, not the
        fixed-site MuJoCo derivative alone. Add the result to other applied
        forces; do not erase a harness/controller's existing force contribution.
        """
        return -np.asarray(length_jacobian).T@np.asarray(tensions)

    def assert_native_force_disabled(self,model):
        import mujoco
        for name in self.names:
            aid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_ACTUATOR,name)
            if aid<0:raise ValueError('Missing model actuator '+name)
            if model.actuator_gaintype[aid]!=0 or model.actuator_biastype[aid]!=0 or np.any(model.actuator_gainprm[aid]!=0):
                raise ValueError('Disable the native reference force explicitly before applying BPA forces: '+name)
