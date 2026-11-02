"""Explicitly authorized electrical design probes only. Never starts a body."""
import argparse,json,time,math
from pathlib import Path
from dataclasses import asdict
import numpy as np
from table2_glif import ActivationDomain,derive_neuron,derive_synapse,constant_threshold_rate_hz,average_reset_synapse_uS,GLIFPopulationRuntime
from gait_design_profiles import write_profiles

def write(p,v):p.write_text(json.dumps(v,indent=2),encoding='utf-8')

def constant_current(d,current,duration=4.,pop=16,edge=None,dt=.00002):
    neurons=[d] if edge is None else [d,d];edges=[] if edge is None else [(0,1,edge)]
    a=GLIFPopulationRuntime(neurons,edges,pop,dt,allow_design_probe=True);before=np.zeros(len(neurons));times=[];threshold=[];trace=[];started=time.perf_counter()
    steps=round(duration/dt);settle=round((duration-1.)/dt)
    for k in range(steps):
        a.step([current] if edge is None else [current,0.])
        if k==settle-1:before=a.total_events.copy()
        if a.events[0,0]:times.append(a.time_s);threshold.append(a.theta[0,0])
        if k%100==0:trace.append([a.time_s,*a.U[:,0],*a.theta[:,0],*a.total_events])
    rates=(a.total_events-before)/pop
    return dict(rate_hz=rates.tolist(),wall_s=time.perf_counter()-started,finite=bool(np.isfinite(a.U).all() and np.isfinite(a.theta).all()),actual_events=a.total_events.tolist(),spike_times=np.array(times),spike_thresholds=np.array(threshold),trace=np.array(trace)),a

def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    protocol=dict(scope='Pure electrical GLIF/source-ledger design probes only; no optimizer,Hunt,curriculum,body,or existing-controller mutation.',source='Szczecinski2020 Table2 PDFp11/Eqs23,26–29',activation=dict(fmax_hz=100.,range_mV=20.,theta0_mV=1.),leak_uS=1.,delta=.01,step_s=.00002,populations=[1,16],gates=dict(constant_current_rate_error_hz_max=.5,paper_derivation_abs_error_max=1e-12,transmission_gain_relative_error_max=.03,transient_theta_RMSE_mV_max=.03,transient_tau_relative_error_max=.10),warning='Preregistered design expectations, not evidence. Fmax is not a refractory clamp; transient/loaded-domain validity must be measured. Literal Table29 and finite-delta correction kept distinct.',training_allowed=False)
    write(out/'PROTOCOL.json',protocol);domain=ActivationDomain(100.,20.,1.);m0=derive_neuron(domain,0.,1.,.5);ma=derive_neuron(domain,-5.,1.,.5);syn=derive_synapse(domain,1.,1.,160.,.01)
    expected=dict(m0_bias=.5,m0_tau=.2,adaptive_bias=1/7,adaptive_tau=.7,adaptive_theta_tau=1.75,syn_tau=-1/(100*math.log(.01)))
    values=dict(m0_bias=m0.bias_nA,m0_tau=m0.tau_mem_s,adaptive_bias=ma.bias_nA,adaptive_tau=ma.tau_mem_s,adaptive_theta_tau=ma.tau_theta_s,syn_tau=syn.tau_s)
    assert all(abs(values[k]-v)<1e-12 for k,v in expected.items())
    unit=dict(published_examples_exact=True,values=values,paper_synapse_peak_uS=syn.selected_peak_uS,paper_synapse_rounding='.658uS',capacitance_units='200nF and700nF for1uS leak',no_ionic_support=True)
    rejected=[]
    for label,fn in [('ionic',lambda:derive_neuron(domain,0.,1.,.5,ionic_conductance_uS=12.)),('bad_m',lambda:derive_neuron(domain,2.,1.,.5)),('bad_delta',lambda:derive_synapse(domain,1.,1.,160.,1.)),('negative_G',lambda:derive_synapse(domain,1.,1.,10.,.01)),('unregistered_probe',lambda:GLIFPopulationRuntime([m0],[],1,.00002))]:
        try:fn()
        except (ValueError,NotImplementedError,RuntimeError):rejected.append(label)
        else:raise AssertionError(label)
    unit['invalid_domains_rejected']=rejected;write(out/'UNIT_DERIVATION.json',unit)
    rows=[];data={}
    for fraction in (0.,.25,.5,1.):
        r,a=constant_current(m0,20*fraction);exact=constant_threshold_rate_hz(m0,20*fraction);row=dict(case='m0_'+str(fraction),current_nA=20*fraction,predicted_exact_rate_hz=exact,measured_rate_hz=r['rate_hz'][0],approximate_rate_hz=100*fraction,passed=abs(r['rate_hz'][0]-exact)<=.5,wall_s=r['wall_s'],actual_events=r['actual_events']);rows.append(row)
        if fraction in (0.,1.):data[row['case']]=r['trace']
    adaptive,a=constant_current(ma,20.,duration=6.);t=adaptive['spike_times'];v=adaptive['spike_thresholds'];pred=ma.theta_star_approx_mV+(1-ma.theta_star_approx_mV)*np.exp(-t/.5);keep=(t>.2)&(t<2.5);rmse=float(np.sqrt(np.mean((v[keep]-pred[keep])**2)));delta=v-ma.theta_star_approx_mV;fit=(t>.3)&(t<2.)&(delta>1e-4);fitted=-1/np.polyfit(t[fit],np.log(delta[fit]),1)[0]
    transient=dict(theta_RMSE_mV=rmse,source_target_tau_s=.5,fitted_tau_s=float(fitted),last_second_rate_hz=adaptive['rate_hz'][0],passed=bool(rmse<=.03 and abs(fitted/.5-1)<=.10));data['adaptive_trace']=adaptive['trace'];data['adaptive_spike_t']=t;data['adaptive_spike_theta']=v;write(out/'TRANSIENT.json',transient)
    transfer=[]
    for N in (1,16):
        r,a=constant_current(m0,20.,duration=4.,pop=N,edge=syn);pre,post=r['rate_hz'];row=dict(population_size=N,pre_rate_hz=pre,post_rate_hz=post,gain=post/pre if pre else None,target_gain=1.,passed=bool(pre and abs(post/pre-1)<=.03),wall_s=r['wall_s'],synapse_rule='Literal Table29/reset-to-peak, not finite-delta corrected');transfer.append(row);data['transfer_'+str(N)]=r['trace']
    write(out/'ELECTRICAL_RESULTS.json',dict(unit=unit,current_cases=rows,transient=transient,transmission=transfer,all_registered_design_gates_pass=bool(all(r['passed'] for r in rows) and transient['passed'] and all(r['passed'] for r in transfer)),training_allowed=False,scope='Passing or failing these electrical examples never releases the all-family training gate. No Gait pathway or body pass is implied.'))
    np.savez_compressed(out/'SELECTED_ELECTRICAL_TRACES.npz',**data)
    write_profiles(out/'source_profiles')
    print(json.dumps(dict(output=str(out),current_cases=rows,transient=transient,transmission=transfer,training_allowed=False)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run-design-probes',action='store_true',required=True);p.add_argument('--output',required=True);a=p.parse_args();run(a.output)
