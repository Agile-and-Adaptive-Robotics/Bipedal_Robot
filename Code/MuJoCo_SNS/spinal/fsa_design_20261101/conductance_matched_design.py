"""Explicit extension beyond literal Table29: match actual GLIF firing rate.

No numerical optimization or empirical gain multiplier. A monotone scalar
equation from Eq13 determines mean conductance, retaining voltage-dependent
synaptic current and conductance-dependent membrane time constant.
Applicable here only to fixed-threshold pure GLIF neurons (m=0).
"""
from dataclasses import replace
import math
from table2_glif import derive_synapse,constant_threshold_rate_hz

def rate_matched_synapse(pre_domain,post_design,gain,reversal_mV,delta):
    if post_design.m!=0:raise NotImplementedError('Adaptive theta needs its implicit Eq14 counterpart before this correction is applicable')
    if gain<=0:raise NotImplementedError('This transmission design is excitatory only; signed fan-in needs an explicit baseline/function equation')
    if reversal_mV<=post_design.activation.theta0_mV:raise ValueError('Reversal cannot drive threshold')
    target=gain*pre_domain.fmax_hz
    def rate(g):return constant_threshold_rate_hz(post_design,0.,conductance_uS=g,conductance_times_reversal_nA=g*reversal_mV)
    low=0.;high=post_design.leak_uS
    for _ in range(80):
        if rate(high)>=target:break
        high*=2
    else:raise ValueError('No finite bracket for specified rate')
    for _ in range(80):
        middle=(low+high)/2
        if rate(middle)<target:low=middle
        else:high=middle
    mean=(low+high)/2
    # Reuse only the documented tau/delta choice, not Table29's Gmax.
    template=derive_synapse(pre_domain,post_design.leak_uS,gain,reversal_mV,delta)
    peak=mean/(template.tau_s*template.fmax_hz*(1-template.actual_delta_at_fmax))
    return replace(template,target_mean_uS=mean,peak_table29_uS=template.peak_table29_uS,
        peak_finite_delta28_uS=peak,selected_peak_uS=peak,
        calibration_mode='full_conductance_Eq13_rate_match_at_Fmax_plus_Eq28_NOT_literal_Table29')

def high_rate_meanfield_conductance_uS(post_design,target_rate_hz,reversal_mV):
    theta=post_design.activation.theta0_mV
    if post_design.m!=0:raise NotImplementedError('This check is for constant threshold')
    if reversal_mV<=theta/2:raise ValueError('Invalid reversal for positive rate gain')
    return post_design.capacitance_nF/1000*theta*target_rate_hz/(reversal_mV-theta/2)
