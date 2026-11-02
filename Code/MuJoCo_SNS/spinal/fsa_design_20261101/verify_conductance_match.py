"""Design-only scalar/clock, range, and fan-in checks. No model integration."""
import argparse,json,time
from pathlib import Path
import numpy as np
from table2_glif import ActivationDomain,derive_neuron,derive_synapse,GLIFPopulationRuntime
from conductance_matched_design import rate_matched_synapse,high_rate_meanfield_conductance_uS
from verify_table2 import constant_current,write

def parity(d,s):
    a=GLIFPopulationRuntime([d,d],[(0,1,s)],16,.00002,allow_design_probe=True);b=GLIFPopulationRuntime([d,d],[(0,1,s)],16,.00002,allow_design_probe=True);maximum=0.;event_match=True
    for k in range(5000):
        u=[20.,0.];a.step(u);b.step_reference(u);maximum=max(maximum,float(abs(a.U-b.U).max()),float(abs(a.theta-b.theta).max()),float(abs(a.gs-b.gs).max()));event_match=event_match and np.array_equal(a.events,b.events)
    return dict(maximum_state_error=maximum,all_events_exact=bool(event_match),passed=bool(maximum<1e-10 and event_match))

def fanin(d,s):
    rows=[]
    for fractions in ((1.,1.),(.25,.75),(.5,0.)):
        a=GLIFPopulationRuntime([d,d,d],[(0,2,s),(1,2,s)],16,.00002,allow_design_probe=True);before=None
        for k in range(200000):
            a.step([20*fractions[0],20*fractions[1],0.])
            if k==149999:before=a.total_events.copy()
        rates=(a.total_events-before)/16;target=.5*(rates[0]+rates[1]);rows.append(dict(input_fractions=fractions,rates_hz=rates.tolist(),target_mean_hz=target,relative_error=float(abs(rates[2]/target-1)),passed=bool(abs(rates[2]/target-1)<=.03)))
    return rows

def run(out):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    if (out/'PROTOCOL.json').exists() and not (out/'RESULTS.json').exists() and not (out/'CLOCK_GUARD_REJECTION.json').exists():
        write(out/'CLOCK_GUARD_REJECTION.json',dict(prior_protocol=json.loads((out/'PROTOCOL.json').read_text()),failure='R5/theta1 has tau_mem50ms;20us exceeds paper policy min(tau_mem)*1e-4=5us. Runtime rejected before this contrast; no successful contrast result existed.',resolution='Use5us in the R5 electrical contrast; keep runtime guard unchanged. Other paper-domain examples retain20us.'))
    write(out/'PROTOCOL.json',dict(scope='Design-only correction beyond literal Table29; no training/body/controller mutation.',source_equations='Eq13 exact constant-threshold rate including Gtotal and G*E, then Eq28 reset-synapse average. Monotone analytical root, not empirical fit.',preserved_negative='../gait_table2/ELECTRICAL_RESULTS.json',registered_gates=dict(gain_relative_error_max=.03,scalar_state_error_max=1e-10,scalar_events_exact=True),cases=['paperR20/theta1 range25/50/100%,1cell/16cells','dt20us vs10us literal transmission','two independent0.5gain inputs','Gait-likeR5/theta1/E40 domain contrast; not an accepted family profile'],training_allowed=False))
    domain=ActivationDomain(100.,20.,1.);d=derive_neuron(domain,0.,1.,.5);literal=derive_synapse(domain,1.,1.,160.,.01);matched=rate_matched_synapse(domain,d,1.,160.,.01);checks=parity(d,literal);clock=[]
    for dt in (.00002,.00001):
        r,_=constant_current(d,20.,edge=literal,dt=dt);clock.append(dict(step_s=dt,rates_hz=r['rate_hz'],gain=r['rate_hz'][1]/r['rate_hz'][0]))
    rows=[]
    for N in (1,16):
        for frac in (.25,.5,1.):
            r,_=constant_current(d,20*frac,pop=N,edge=matched);pre,post=r['rate_hz'];rows.append(dict(population_size=N,input_fraction=frac,pre_rate_hz=pre,post_rate_hz=post,gain=post/pre,passed=bool(abs(post/pre-1)<=.03)))
    half=rate_matched_synapse(domain,d,.5,160.,.01);fan=fanin(d,half)
    write(out/'PARTIAL_PAPER_DOMAIN.json',dict(parity=checks,clock=clock,range=rows,fanin=fan))
    small=ActivationDomain(100.,5.,1.);ds=derive_neuron(small,0.,1.,.05);sm=rate_matched_synapse(small,ds,1.,40.,.01);small_rows=[]
    for frac in (.25,.5,1.):
        r,_=constant_current(ds,5*frac,edge=sm,dt=.000005);pre,post=r['rate_hz'];small_rows.append(dict(input_fraction=frac,step_s=.000005,pre_rate_hz=pre,post_rate_hz=post,gain=post/pre,passed=bool(abs(post/pre-1)<=.03)))
    report=dict(parity=checks,literal_clock_comparison=clock,literal_mean_G_uS=literal.target_mean_uS,derived_mean_G_uS=matched.target_mean_uS,derived_high_rate_meanfield_G_uS=high_rate_meanfield_conductance_uS(d,100.,160.),literal_peak_G_uS=literal.selected_peak_uS,derived_peak_G_uS=matched.selected_peak_uS,range_results=rows,fanin_results=fan,gait_like_R5_contrast=small_rows,paper_domain_pass=bool(checks['passed'] and all(r['passed'] for r in rows+fan)),gait_like_domain_pass=all(r['passed'] for r in small_rows),training_allowed=False,scope='An electrical function proof is not an actual-family pathway or mechanical proof. NaP/Ca remains excluded.')
    write(out/'RESULTS.json',report);print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();run(a.output)
