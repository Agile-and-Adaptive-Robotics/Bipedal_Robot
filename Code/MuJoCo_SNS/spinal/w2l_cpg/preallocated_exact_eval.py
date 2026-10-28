"""Exact-order SNS stepping with preallocated one-hot input/output maps.

Extends the validated sparse-edge Numba backend. Synaptic np.sum ordering,
NaP equations, and timestep remain unchanged. General dense maps fall back
to NumPy matmul; the family builds use one input per driven neuron and one
neuron per output, for which dense products are unnecessary.
"""
import numpy as np
from exact_fast_eval import ExactFastForward,_edges

class PreallocatedExactForward(ExactFastForward):
    def __init__(self,compiled):
        if compiled.gated and type(compiled).__name__!='SNS_NumpyFixedTau':
            raise ValueError('This backend implements fixed-tau NaP only; voltage-dependent paper kinetics require their own verified kernel.')
        super().__init__(compiled)
        c=self.c
        ip,iq=np.nonzero(c.input_connectivity)
        op,oq=np.nonzero(c.output_voltage_connectivity)
        self.input_onehot=(len(set(ip))==len(ip) and np.all(c.input_connectivity[ip,iq]==1.))
        self.output_onehot=(len(op)==c.output_voltage_connectivity.shape[0] and len(set(op))==len(op)
                            and np.all(c.output_voltage_connectivity[op,oq]==1.))
        self.ip,self.iq=ip,iq
        self.oq=oq
        self.i_app=np.zeros_like(c.V)
        self.input_values=np.empty(len(iq))
        self.syn_ge_sum=np.empty_like(c.V);self.syn_g_sum=np.empty_like(c.V)
        self.last_voltage=np.empty_like(c.V)
        self.output_values=np.empty(c.output_voltage_connectivity.shape[0])

    def forward(self,x=None):
        c=self.c
        np.copyto(self.last_voltage,c.V);c.V_last=self.last_voltage
        if x is None:i_app=0
        elif self.input_onehot:
            np.take(x,self.iq,out=self.input_values)
            self.i_app[self.ip]=self.input_values;i_app=self.i_app
        else:i_app=np.matmul(c.input_connectivity,x)
        _edges(c.V_last,self.post,self.pre,self.gmax,self.elo,self.erange,self.reversal,self.g,self.ge)
        np.sum(self.ge,axis=1,out=self.syn_ge_sum)
        np.sum(self.g,axis=1,out=self.syn_g_sum)
        i_syn=self.syn_ge_sum-c.V_last*self.syn_g_sum
        if c.gated:
            a_inf=1/(1+c.k_a*np.exp(c.slope_a*(c.e_a-c.V_last)))
            b_inf=1/(1+c.k_b*np.exp(c.slope_b*(c.e_b-c.V_last)))
            c_inf=1/(1+c.k_c*np.exp(c.slope_c*(c.e_c-c.V_last)))
            tau_c=c.tau_max_c*c_inf*np.sqrt(c.k_c*np.exp(c.slope_c*(c.e_c-c.V_last)))
            c.b_gate_last,c.c_gate_last=np.copy(c.b_gate),np.copy(c.c_gate)
            c.b_gate=c.b_gate_last+c.dt*((b_inf-c.b_gate_last)/c.tau_max_b)
            c.c_gate=c.c_gate_last+c.dt*((c_inf-c.c_gate_last)/tau_c)
            i_ion=c.g_ion*(a_inf**c.pow_a)*(c.b_gate**c.pow_b)*(c.c_gate**c.pow_c)*(c.e_ion-c.V_last)
            i_gated=np.sum(i_ion,axis=0)
            c.V=c.V_last+c.time_factor_membrane*(-c.g_m*(c.V_last-c.V_rest)+c.i_b+i_syn+i_app+i_gated)
        else:
            c.V=c.V_last+c.time_factor_membrane*(-c.g_m*(c.V_last-c.V_rest)+c.i_b+i_syn+i_app)
        if self.output_onehot:
            np.take(c.V,self.oq,out=self.output_values);c.outputs=self.output_values
        else:c.outputs=np.matmul(c.output_voltage_connectivity,c.V)
        return c.outputs
