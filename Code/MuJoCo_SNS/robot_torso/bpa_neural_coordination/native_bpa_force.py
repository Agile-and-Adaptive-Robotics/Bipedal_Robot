"""Opt-in native-actuator pressure force adapter, without neural-state edits.

The caller must provide a complete RoutedBPABank configuration. No hardware
dimensions, pressure response, or actuator gains are inferred here. A user-bias
callback returns negative tensile force, preserving MuJoCo's native moment and
force-sensor conventions. The caller explicitly updates pressure once per step;
the callback itself has no state integration or optimization.

Only use this adapter when native MuJoCo path derivatives represent all site
motion. AH001's externally updated DynamicRoutes additionally need their
moving/conditional-point generalized-force correction and a separate parity
gate; this adapter alone is not a complete CAD pressure interface.
"""
from __future__ import annotations
import numpy as np
import mujoco
from cad_interface import RoutedBPABank
from actuator_interface import inverse_pressure,force_at_pressure


class NativeBPAForce:
    def __init__(self,model,configuration,dt,disabled_names=()):
        self.model=model
        self.names=list(configuration['actuators'])
        self.ids=np.array([model.actuator(n).id for n in self.names],dtype=int)
        self.disabled=np.array([model.actuator(n).id for n in disabled_names],dtype=int)
        if set(self.ids)&set(self.disabled):raise ValueError('Configured and disabled actuator sets overlap')
        if len(set(self.ids))!=len(self.ids):raise ValueError('Duplicate actuator mapping')
        if set(self.ids)|set(self.disabled)!=set(range(model.nu)):
            raise ValueError('Explicitly configure or disable every actuator; no silent biological/BPA mixture')
        # A pressure-force magnitude is physical path tension only for unit gear.
        if not np.allclose(model.actuator_gear[self.ids,0],1.,rtol=0,atol=1e-14):
            raise ValueError('This adapter requires unit positive tendon gear')
        if np.any(model.actuator_trntype[self.ids]!=mujoco.mjtTrn.mjTRN_TENDON):
            raise ValueError('Only explicitly routed tendon actuators are supported')
        self.bank=RoutedBPABank(self.names,configuration,dt)
        self.configuration=configuration
        uncertainty=configuration.get('uncertainty',{})
        allowed={'force_gain','directional_gain','direction_tau_s','velocity_scale_m_s'}
        if set(uncertainty)-allowed:raise ValueError('Unconsumed actuator-uncertainty parameters')
        self.force_gain=float(uncertainty.get('force_gain',1.))
        self.direction_gain=float(uncertainty.get('directional_gain',0.))
        self.direction_tau=float(uncertainty.get('direction_tau_s',.08))
        self.velocity_scale=float(uncertainty.get('velocity_scale_m_s',.03))
        if not(self.force_gain>0 and 0<=self.direction_gain<1 and self.direction_tau>0 and self.velocity_scale>0):raise ValueError('Invalid declared uncertainty parameters')
        self.direction_state=np.zeros(len(self.ids));self.effective_gain=np.ones(len(self.ids))*self.force_gain
        self.nominal_force=np.zeros(len(self.ids))
        self.uncertain=self.force_gain!=1. or self.direction_gain!=0.
        if self.uncertain:
            self.bank.force_parameters=self.bank.parameters.copy()
            self.bank.force_parameters[:,3]=self.bank.parameters[:,3]*self.force_gain
        self.lookup={int(a):i for i,a in enumerate(self.ids)}
        self.previous=mujoco.get_mjcb_act_bias()
        self.touched=np.r_[self.ids,self.disabled]
        self.original={k:getattr(model,'actuator_'+k)[self.touched].copy()
                       for k in ('gainprm','biasprm','gaintype','biastype','forcelimited')}
        model.actuator_gaintype[self.touched]=mujoco.mjtGain.mjGAIN_FIXED
        model.actuator_gainprm[self.touched]=0.
        model.actuator_biasprm[self.touched]=0.
        model.actuator_biastype[self.ids]=mujoco.mjtBias.mjBIAS_USER
        model.actuator_biastype[self.disabled]=mujoco.mjtBias.mjBIAS_NONE
        # The characterized pressure-force law supplies the capacity. Keeping
        # the old biological force limiter would silently introduce a second cap.
        model.actuator_forcelimited[self.touched]=False
        def callback(m,d,index):
            if m is self.model and index in self.lookup:
                return -float(self.bank.force[self.lookup[index]])
            return float(self.previous(m,d,index)) if self.previous is not None else 0.
        self.callback=callback
        mujoco.set_mjcb_act_bias(self.callback)
        self.closed=False

    def advance(self,data,normalized_commands):
        """Call once before physics, after refreshing current path geometry."""
        u=np.asarray(normalized_commands)
        if u.shape!=(self.model.nu,):raise ValueError('Expected one neural command per native actuator')
        lengths=data.actuator_length[self.ids]
        if self.direction_gain:
            target=np.tanh(data.actuator_velocity[self.ids]/self.velocity_scale)
            self.direction_state+=(1.-np.exp(-self.bank.dt/self.direction_tau))*(target-self.direction_state)
            self.effective_gain[:]=self.force_gain*(1.+self.direction_gain*self.direction_state)
        if self.uncertain:self.bank.force_parameters[:,3]=self.bank.parameters[:,3]*self.effective_gain
        try:self.bank.step(u[self.ids],lengths)
        except ValueError as error:
            invalid=np.flatnonzero(~np.isfinite(self.bank.force))
            details=[]
            for i in invalid:
                p=self.bank.force_parameters[i];maximum=p[0]+p[2]+p[3]*(self.bank.pressure[i]/620.)/p[4]
                details.append(dict(name=self.names[i],route_length_m=float(lengths[i]),zero_stretch_route_rest_length_m=float(p[0]+p[2]),pressure_supported_maximum_route_length_m=float(maximum),excess_over_pressure_supported_length_mm=float(1000*(lengths[i]-maximum)),actual_pressure_kpa=float(self.bank.pressure[i])))
            raise ValueError(str(error)+'; failing routes: '+str(details)) from error
        if self.uncertain:
            for i in range(len(self.ids)):self.nominal_force[i]=force_at_pressure(lengths[i],self.bank.pressure[i],self.bank.parameters[i])
        else:np.copyto(self.nominal_force,self.bank.force)
        return self.bank.force

    def precharge(self,data,normalized_commands):
        """Explicit steady pressure initial state; not an inferred measurement."""
        u=np.asarray(normalized_commands)[self.ids]
        for i,a in enumerate(self.ids):
            length=float(data.actuator_length[a]);p=self.bank.parameters[i]
            if self.bank.mode=='inverse_force':
                command,flag=inverse_pressure(float(np.clip(u[i],0,1))*self.bank.scales[i],length,p)
            else:command=float(np.clip(u[i],0,1))*p[8];flag=0
            self.bank.command[i]=command;self.bank.flags[i]=flag
            self.bank.pressure[i]=command;self.bank.history[i,:]=command
            self.nominal_force[i]=force_at_pressure(length,command,p)
            self.bank.force[i]=force_at_pressure(length,command,self.bank.force_parameters[i])
        self.bank.step_count=self.bank.history.shape[1]
        if not np.isfinite(self.bank.force).all():raise ValueError('Precharge length outside configured domain')

    def configure_response(self,configuration):
        """Change declared response uncertainty without changing stored state.

        Geometry, force-reference mapping, and stiffness must match. Pressure,
        command history, and direction memory are retained, not reinitialized.
        """
        if configuration.get('command_mode','direct_pressure')!=self.bank.mode:raise ValueError('Cannot switch pressure command interpretation')
        fixed=('diameter_mm','resting_length_m','maximum_contraction_fraction','nonmuscle_length_m','series_stiffness_N_m','parallel_count','pressure_max_kpa','force_command_scale_N')
        for i,name in enumerate(self.names):
            a=self.configuration['actuators'][name];b=configuration['actuators'][name]
            if any(a.get(k)!=b.get(k) for k in fixed):raise ValueError('Response switch cannot change geometry/energy/reference mapping: '+name)
            delay=round(b['delay_s']/self.bank.dt)
            if delay>=self.bank.history.shape[1]:raise ValueError('Delay history capacity insufficient')
            if b['tau_fill_s']<=0 or b['tau_vent_s']<=0 or delay<0:raise ValueError('Invalid pressure response')
            self.bank.fill[i]=b['tau_fill_s'];self.bank.vent[i]=b['tau_vent_s'];self.bank.delay[i]=delay
        u=configuration.get('uncertainty',{})
        if set(u)-{'force_gain','directional_gain','direction_tau_s','velocity_scale_m_s'}:raise ValueError('Unconsumed uncertainty setting')
        self.force_gain=float(u.get('force_gain',1.));self.direction_gain=float(u.get('directional_gain',0.));self.direction_tau=float(u.get('direction_tau_s',.08));self.velocity_scale=float(u.get('velocity_scale_m_s',.03))
        if not(self.force_gain>0 and 0<=self.direction_gain<1 and self.direction_tau>0 and self.velocity_scale>0):raise ValueError('Invalid uncertainty setting')
        self.uncertain=self.force_gain!=1. or self.direction_gain!=0.
        self.effective_gain[:]=self.force_gain*(1.+self.direction_gain*self.direction_state)
        self.bank.force_parameters=self.bank.parameters.copy() if self.uncertain else self.bank.parameters
        if self.uncertain:self.bank.force_parameters[:,3]=self.bank.parameters[:,3]*self.effective_gain
        self.configuration=configuration

    def close(self):
        if self.closed:return
        mujoco.set_mjcb_act_bias(self.previous)
        for k,v in self.original.items():getattr(self.model,'actuator_'+k)[self.touched]=v
        self.closed=True

    def __enter__(self):return self
    def __exit__(self,*exc):self.close()
