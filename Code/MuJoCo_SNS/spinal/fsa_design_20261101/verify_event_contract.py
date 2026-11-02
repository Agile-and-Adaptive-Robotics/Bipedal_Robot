"""Physical-clock, actual public-step, and consumption guards; no body run."""
import os,sys,copy,json,argparse,hashlib
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).parent;sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
import numpy as np
from fsa_design_20261101.event_conductance_port import EventConductanceInputPort,default_port_profile,install_port_profile,steady_carrier_fraction
from graded_profile import default_profile,apply_profile,population_names

def main(family,out):
    os.environ['AARL_NET']='' if family=='s3k' else family
    files=[Path(__file__),HERE/'event_conductance_port.py',HERE/'graded_profile.py',HERE/'table2_glif.py']+[HERE.parent/q for q in ('build_network.py','build_network_w2lvar.py','build_network_syn6.py','params.py','muscle_map.py','synergy_basis.npz','fsa_results/fsa_backsolve.npz')];before={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in files}
    pin=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\gait_table2\source_profiles\syn6_DESIGN_PROFILE.json');expected_basis=json.loads(pin.read_text())['source_NMF_basis']['sha256'];assert before[str(HERE.parent/'synergy_basis.npz')]==expected_basis,'Protected K6 basis changed; do not substitute/re-fit it'
    import mujoco,build_network as B
    m=mujoco.MjModel.from_xml_path(str(HERE.parents[3]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'))
    n=B.build([m.actuator(i).name for i in range(m.nu)]);apply_profile(n,default_profile(family),family=family);raw=n.compiled;names=population_names(n);p=default_port_profile(family,True)
    a=install_port_profile(n,p,family=family);b=EventConductanceInputPort(raw,names,p['input_names']);x=np.zeros(len(n.inputs));x[b.ports[0][0]]=.8
    for _ in range(64):b.forward(x);n.step(x);assert np.array_equal(a.V,b.V)
    g={};v=a.peaks[0];a.peaks[0]*=.99
    try:a.verify_consumption()
    except ValueError:g['changed_peak_rejected']=True
    else:raise AssertionError('Changed peak accepted')
    a.peaks[0]=v;a.verify_consumption();a.inner.c.g_max_non[a.edge_specs[0][1],a.indices[0]]=1e-9
    try:a.forward(x)
    except ValueError:g['graded_axon_reintroduction_rejected']=True
    else:raise AssertionError('Hidden graded axon accepted')
    a.inner.c.g_max_non[a.edge_specs[0][1],a.indices[0]]=0.;a.verify_consumption();n.compiled=raw
    try:install_port_profile(n,dict(p,unused_clock_scale=.2),family=family)
    except ValueError:g['unused_field_rejected']=True
    else:raise AssertionError('Unused key accepted')
    def carrier(dt,mode='pulse'):
        r=copy.deepcopy(raw);r.dt=dt;r.time_factor_membrane*=dt/raw.dt
        c=EventConductanceInputPort(r,names,p['input_names']);gs=[]
        for k in range(round(1/dt)):
            t=k*dt;u=np.zeros(len(n.inputs));u[c.ports[0][0]]=({'zero':0.,'full':1.}.get(mode,0. if t<.2 else (.8 if t<.3 else (1. if t<.8 else 0.))));c.forward(u);gs.append(c.gs.copy())
        e=np.array(c.events).reshape(-1,3)
        if len(e):e=e[np.lexsort((e[:,0],e[:,2],e[:,1]))]
        return c,np.array(gs),e
    a,ga,ea=carrier(.002);b,gb,eb=carrier(.001);err=float(abs(ea-eb).max());gerr=float(abs(ga-gb.reshape(len(ga),2,-1).mean(axis=1)).max());assert err<1e-10 and gerr<1e-10
    z,gz,ez=carrier(.002,'zero');assert len(ez)==0 and not gz.any()
    f,gf,ef=carrier(.002,'full');rate=float(np.count_nonzero((ef[:,0]>=.2)&(ef[:,1]==0))/(f.N*.8));assert abs(rate-500)<1.
    dcx=np.linspace(0,1,1001);dcy=np.array([steady_carrier_fraction(float(x)) for x in dcx]);error=dcy-dcx;assert error.min()>=-.025-1e-12 and error.max()<=.01/.99+1e-12
    report=dict(family=family,actual_SpinalCord_step64_exact=True,mutation_guards=g,carrier_only_partition=dict(coarse_dt_s=.002,fine_dt_s=.001,event_time_max_error_s=err,mean_G_max_error_uS=gerr,whole_neural_clock_parity_claimed=False),full_input_physical_Hz=rate,zero_input_spikes=len(ez),profile=p,graded_profile=default_profile(family),consumption_ledger=n.fsa_spiking_port_ledger,source_files_sha256={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__),HERE/'event_conductance_port.py',HERE/'graded_profile.py',HERE/'table2_glif.py',HERE.parent/'build_network.py',HERE.parent/'params.py']},steady_DC_contract=dict(analytical_error_bounds_full_scale=[-.025,.01/.99],sampled1001_error_range_full_scale=[float(error.min()),float(error.max())],silent_fraction_upper_inclusive=.025,small_signal_transient_bandwidth_verified=False),no_training_or_body=True)
    report['source_files_sha256']={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in files}
    assert before==report['source_files_sha256'];report['source_files_unchanged']=True;report['K6_basis_pin']=dict(path=str(pin),expected_sha256=expected_basis,basis_unchanged=True,no_NMF_execution=True)
    out=Path(out);out.mkdir(parents=True,exist_ok=True);(out/'API_CLOCK_GUARDS.json').write_text(json.dumps(report,indent=2));(out/'PORT_PROFILE.json').write_text(json.dumps(p,indent=2));print(json.dumps({k:report[k] for k in ('family','actual_SpinalCord_step64_exact','mutation_guards','carrier_only_partition','full_input_physical_Hz','steady_DC_contract','source_files_unchanged')}),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--family',required=True,choices=['s3k','w2lvar','syn6']);a.add_argument('--output',required=True);v=a.parse_args();main(v.family,v.output)
