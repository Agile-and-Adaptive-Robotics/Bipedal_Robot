"""Actual compiled-family, opt-in electrical port checks; no body or training."""
import os,sys,argparse,json,copy,math,hashlib,time
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).parent;sys.path[:0]=[str(HERE.parent),str(HERE.parent/'w2l_cpg')]
import numpy as np
from graded_profile import default_profile,apply_profile,population_names,field_digest,assert_campaign_allowed
from event_conductance_port import EventConductanceInputPort
from sparse_pairwise_eval import SparsePairwiseForward

def write(p,o):p.write_text(json.dumps(o,indent=2),encoding='utf-8')
def all_equal(a,b):return a.keys()==b.keys() and all(np.array_equal(a[k],b[k]) for k in a)

def run(family,out):
    assert_campaign_allowed('unit_and_individual_pathway_probes');os.environ['AARL_NET']='' if family=='s3k' else family
    import mujoco,build_network as B
    from params import G
    model=HERE.parents[3]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
    m=mujoco.MjModel.from_xml_path(str(model));acts=[m.actuator(i).name for i in range(m.nu)];net=B.build(acts);names=population_names(net);original_digest=field_digest(net.compiled);original=copy.deepcopy(net.compiled);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    profile=default_profile(family);protocol=dict(family=family,mode='Actual compiled neural family only; static MuJoCo compile supplies92 actuator names; no dynamics.',grade_profile=profile,spiking_source_names=['Ia_soleus_r','Ia_soleus_l'],carrier=dict(fmax_physical_hz=500.,R_mV=20.,theta0_mV=1.,m=0.,delta=.01,cells_per_input=16,host_physical_step_s=float(net.compiled.dt),integration='Exact event-time LIF and reset-synapse integral; host retains existing graded/NaP numerical update'),rationale='500Hz is an explicit internal-carrier probe choice exceeding the full-input latency bound from10deg/7.16Hz budget; it is not a valve rate or proof at arbitrarily small signals.',gates=dict(disabled_128_exact=True,outgoing_off_nonencoded_exact=True,shadow_poison_nonencoded_exact=True,complete_state_exact=True,whole_network_pulse_MN_delta_min_mV=1e-6,encoder_phase_lag_max_deg=10.,normalized_conductance_RMSE_max=.05),source_original_digest=original_digest,training_allowed=False)
    write(out/'PROTOCOL.json',protocol)
    a=SparsePairwiseForward(copy.deepcopy(original));b=EventConductanceInputPort(original,names,protocol['spiking_source_names'],enabled=False);rng=np.random.default_rng(1102)
    for _ in range(128):
        x=rng.uniform(0,.2,len(net.inputs));a.forward(x);b.forward(x);assert np.array_equal(a.V,b.V)
    assert field_digest(net.compiled)==original_digest
    ledger=apply_profile(net,profile,family=family);raw=net.compiled;selected=protocol['spiking_source_names'];ids=np.array([names.index(n) for n in selected]);non=np.setdiff1d(np.arange(len(names)),ids);post=names.index('MN_soleus_r');port=int(np.flatnonzero(raw.input_connectivity[ids[0]])[0]);data={};rows={}
    def instance(off=False):return EventConductanceInputPort(raw,names,selected,outgoing_off=off)
    for mode,off,pulse in [('base',False,False),('pulse',False,True),('off_base',True,False),('off_pulse',True,True)]:
        c=instance(off);trace=[]
        for k in range(300):
            x=np.zeros(len(net.inputs));x[port]=.8 if pulse and .2<=k*c.dt<.3 else 0.;c.forward(x);trace.append(c.V.copy())
        data[mode]=np.array(trace);rows[mode]=dict(spikes=c.total_counts.tolist(),consumption=c.verify_consumption());data[mode+'_events']=np.array(c.events).reshape(-1,3)
    off_exact=bool(np.array_equal(data['off_base'][:,non],data['off_pulse'][:,non]));effect=float(np.max(abs(data['pulse'][:,post]-data['base'][:,post])));assert off_exact and effect>1e-6
    a=instance();b=instance()
    for k in range(128):
        x=np.zeros(len(net.inputs));x[port]=.5;a.forward(x);b.inner.c.V[ids]=10.;b.forward(x);assert np.array_equal(a.V[non],b.V[non])
    snap=a.export_state()
    for _ in range(32):a.forward(x)
    expected=a.export_state();a.import_state(snap)
    for _ in range(32):a.forward(x)
    assert all_equal(expected,a.export_state())
    # Carrier bandwidth: measure delivered event conductance, not a Vm proxy.
    c=instance();omega=2*math.pi*profile['context']['kinematic_bandwidth_hz'];wave=[]
    for k in range(1200):
        t=k*c.dt;x=np.zeros(len(net.inputs));x[port]=.5+.4*math.sin(omega*t);c.forward(x);wave.append([t,x[port],c.gs[post]])
    wave=np.array(wave);keep=wave[:,0]>=.4;basis=np.c_[np.sin(omega*wave[keep,0]),np.cos(omega*wave[keep,0]),np.ones(keep.sum())];source=np.linalg.lstsq(basis,wave[keep,1],rcond=None)[0];signal=np.linalg.lstsq(basis,wave[keep,2],rcond=None)[0];phase=float(np.degrees(np.angle(np.exp(1j*(math.atan2(signal[1],signal[0])-math.atan2(source[1],source[0]))))));mean_target=sum(e[4] for e in c.edge_specs if e[1]==post);rmse=float(np.sqrt(np.mean((wave[keep,2]/mean_target-wave[keep,1])**2)));band=bool(abs(phase)<=10 and rmse<=.05)
    report=dict(family=family,default_disabled128_exact=True,original_network_unchanged=True,profile_consumed=True,graded_profile_modified_edges=len(ledger['modified_edges']),replaced_graded_axons=len(c.edge_specs),actual_spiking_cells=32,remaining_graded_cells=len(names)-2,outgoing_off_nonencoded_exact=off_exact,pulse_soleus_MN_max_delta_mV=effect,shadow_poison_nonencoded128_exact=True,complete_state32_exact=True,actual_events=rows,phase_deg=phase,normalized_G_RMSE=rmse,bandwidth_gate=band,individual_port_gate=band,training_allowed=False,scope='New actual-family Ia-only encoder port; other afferents and RG/PF/MN remain graded. Not a legacy standing winner, complete7-step family pass, or body result.')
    write(out/'RESULTS.json',report);np.savez_compressed(out/'SELECTED_TRACES.npz',**data,conductance_bandwidth=wave,names=np.array(names),right_soleus_MN_index=post);print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--family',required=True,choices=['s3k','w2lvar','syn6']);p.add_argument('--output',required=True);a=p.parse_args();run(a.family,a.output)
