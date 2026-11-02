"""Opt-in Hunt-bank ankle DF burst shaping; original physics is untouched.

Existing fast excitatory and delayed inhibitory bank edges shape a transient
PF command. Only the direct RG_F→ankle-DF branch is removed. The plantarflexor
RG_E branch, PF antagonistic relays, A/B, and MN/reflex weights are preserved.
"""
import os,sys,json,hashlib,argparse,runpy
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
HERE=Path(__file__).resolve().parent;S=HERE.parent
sys.path[:0]=[str(HERE),str(S/'hunt_training'),str(S/'functional_gait2392'),str(S),str(S/'w2l_cpg')]
import numpy as np

def install(fast,slow,slow_index=3,reciprocal_scale=1.):
    import controller as C
    if not .1<=fast<=8. or not .1<=slow<=16.:raise ValueError('Bounded positive conductances required')
    if slow_index not in (2,3,4,5):raise ValueError('Use an existing slow bank')
    if not .1<=reciprocal_scale<=1.:raise ValueError('Shared bilateral ankle reciprocal inhibition scale in[.1,1]')
    original=C.build;records=[]
    def build(*args,**kwargs):
        n=original(*args,**kwargs);c=n.compiled;changes=[];removed=[]
        for edge,(pre,post) in enumerate(zip(c.pre,c.post)):
            if n.names[pre].startswith('PF_IN_FUNC_ANKLE_Z_') and n.names[post].startswith('PF_FUNC_ANKLE_Z_'):
                c.gmax[edge]*=reciprocal_scale;n.original_g[edge]=c.gmax[edge]
        for side in ('r','l'):
            target='PF_FUNC_ANKLE_Z_N_'+side
            for source,value in [('HUNT_BANK_F_0_'+side,fast),('HUNT_BANK_IN_'+str(6+slow_index)+'_'+side,slow)]:
                ids=np.flatnonzero((c.pre==n.idx[source])&(c.post==n.idx[target]))
                if len(ids)!=1 or int(ids[0]) not in n.pfbank_edges:raise RuntimeError('Reserved Hunt ankle bank edge missing')
                edge=int(ids[0]);n.pf_gain[edge]=value;changes.append(dict(source=source,target=target,gmax=value,reversal_mV=float(c.reversal[edge])))
            ids=np.array([i for i in n.ankle_rg_edges if c.post[i]==n.idx[target]],dtype=int)
            if len(ids)!=1:raise RuntimeError('Exactly one original direct dorsiflexor phase branch expected')
            for edge in ids:c.gmax[edge]=0.;removed.append(n.names[c.pre[edge]]+' -> '+target)
            n.ankle_rg_edges=n.ankle_rg_edges[~np.isin(n.ankle_rg_edges,ids)]
        n.sync();records.append(dict(fast_F_bank_tau_s=.02,slow_F_bank_tau_s=[.02,.05,.1,.2,.4,.7][slow_index],ankle_reciprocal_inhibition_scale=reciprocal_scale,fast_excitation=fast,slow_relay_inhibition=slow,modified_edges=changes,removed_direct_DF_edges=removed,plantarflexor_direct_gain=n.options.get('direct_ankle_pf_gain',0.),scope='Existing reserved Hunt PF-bank connections; no new neurons,physics change,force allocation,or stride clock. Shared bilateral coefficients.'))
        return n
    C.build=build;return records

def pulse_probe(out):
    from sns_toolbox.networks import Network
    from sns_toolbox.neurons import NonSpikingNeuron
    import build_network as BN,params as P
    out.mkdir(parents=True,exist_ok=True);records=[];traces={}
    for slow_tau,duration in [(tau,duration) for tau in (.1,.2) for duration in (.6,.9)]:
        net=Network(name='ideal_RG_ankle_bank_probe')
        for name,tau in [('RG',.02),('fast',.02),('slow',slow_tau),('slow_relay',.02),('PF_DF',P.TAU['pf'])]:net.add_neuron(NonSpikingNeuron(membrane_capacitance=tau,membrane_conductance=1.,bias=0.),name=name)
        for pre,post,g,exc in [('RG','fast',1.,True),('RG','slow',1.,True),('slow','slow_relay',1.,True),('fast','PF_DF',1.5,True),('slow_relay','PF_DF',2.7,False)]:net.add_connection(BN._syn(g,exc),pre,post)
        model=net.compile(dt=.002,backend='numpy');t=np.arange(1000)*.002;v=np.zeros((len(t),5))
        for k,time in enumerate(t):model.V[0]=5. if .2<=time<.2+duration else 0.;model.forward();v[k]=model.V
        command=np.clip(v[:,4]/P.E_HI,0.,1.);active=t[(command>=.1)&(t<.2+duration)];records.append(dict(slow_tau_s=slow_tau,ideal_RG_F_duration_s=duration,peak_PF_voltage_mV=float(v[:,4].max()),time_above_10percent_MN_scale_s=float(len(active)*.002),last_100ms_phase_mean_PF_mV=float(v[(t>=.1+duration)&(t<.2+duration),4].mean()),positive_command_integral_s=float(command.sum()*.002)));traces[str(slow_tau)+'_'+str(duration)]=v
    np.savez_compressed(out/'ANKLE_TEMPORAL_PROBE.npz',t=t,**traces);result=dict(scope='Isolated installed SNS-Toolbox circuit,ideal imposed RG_F voltage. No whole-network or physical performance inference.',coefficients=dict(fast=1.5,slow=2.7),cases=records,module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(out/'ANKLE_TEMPORAL_PROBE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

def network_probe(recording,config,out):
    """Expected-input Hunt C probe retaining the actual antagonistic PF circuit."""
    import controller as C
    z=np.load(recording);cells=z['cells'].tolist();t=z['t'];v=z['voltage'];cfg=json.loads(config.read_text());target=np.load(S/'hunt_training/runs/targets/targets.npz');acts=target['act_names'].tolist();traces={};records={};keep=(t>=5)&(t<13)
    native_build=C.build
    for variant,settings in [('baseline',None),('temporal',(1.5,2.7,3,1.)),('strong_temporal',(8.,14.4,3,1.)),('weaker_reciprocal',(1.5,2.7,3,.25))]:
        C.build=native_build
        if settings is not None:install(*settings)
        n=C.build(acts,options=dict(weights_directory=cfg['weights_directory'],balance=None,biarticular_pair_gain=cfg.get('biarticular_pair_gain',.3),use_timing_afferents=False));u=n.make_inputs();names=['PF_FUNC_ANKLE_Z_'+phase+'_'+side for side in ('r','l') for phase in ('P','N')]+['HUNT_BANK_F_0_r','HUNT_BANK_F_3_r','HUNT_BANK_IN_9_r','PF_IN_FUNC_ANKLE_Z_P_r','PF_IN_FUNC_ANKLE_Z_N_r'];indices=[n.idx[name] for name in names];phasecells=['RG_'+phase+'_'+side for side in ('r','l') for phase in ('E','F')];rg=[n.idx[name] for name in phasecells];source=[cells.index(name) for name in phasecells];result=np.zeros((len(t),len(indices)))
        for k,time in enumerate(t):
            n.compiled.V[rg]=v[max(0,k-1),source] if k else 0.;u[n.ports['DRIVE']]=2.5 if 2.<=time<13. else 0.;n.step(u);result[k]=n.compiled.V[indices]
        traces[variant]=result;records[variant]={name:dict(peak_mV=float(result[keep,i].max()),positive_voltage_integral_mVs=float(np.maximum(result[keep,i],0.).sum()*.002),fraction_above_half_mV=float(np.mean(result[keep,i]>.5))) for i,name in enumerate(names)}
    C.build=native_build;out.mkdir(parents=True,exist_ok=True);np.savez_compressed(out/'NETWORK_TEMPORAL_PROBE.npz',t=t,PF_cells=names,**traces);summary=dict(scope='Actual full neural topology and antagonistic PF relays, recorded RG voltages teacher-forced as expected inputs; no physical model or live feedback. Only declared ankle temporal/reciprocal branches change.',recorded_RG=str(recording),recorded_RG_ranges={name:[float(v[keep,cells.index(name)].min()),float(v[keep,cells.index(name)].max())] for name in phasecells},configuration=str(config),coefficients=dict(temporal=[1.5,2.7,3,1.],strong_temporal=[8.,14.4,3,1.],weaker_reciprocal=[1.5,2.7,3,.25]),results=records);(out/'NETWORK_TEMPORAL_PROBE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--probe',type=Path);parser.add_argument('--network-probe',type=Path);parser.add_argument('--config',type=Path);parser.add_argument('--out',type=Path);parser.add_argument('--protocol',choices=['hunt-air1','hunt-air2','hunt-stand','hunt-ground']);parser.add_argument('--fast',type=float,default=1.5);parser.add_argument('--slow',type=float,default=2.7);parser.add_argument('--slow-index',type=int,choices=(2,3,4,5),default=3);parser.add_argument('--reciprocal-scale',type=float,default=1.);a=parser.parse_args()
    if a.probe:pulse_probe(a.probe)
    elif a.network_probe:
        if a.config is None or a.out is None:parser.error('Expected-input probe requires config and output')
        network_probe(a.network_probe,a.config,a.out)
    else:
        if a.config is None or a.out is None or a.protocol is None:parser.error('Explicit protocol,config,and output required')
        a.out.mkdir(parents=True,exist_ok=True)
        if (a.out/'BACKEND_RESULT.json').exists():raise FileExistsError('Preserve completed trial')
        records=install(a.fast,a.slow,a.slow_index,a.reciprocal_scale);manifest=dict(wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),configuration=str(a.config),configuration_sha256=hashlib.sha256(a.config.read_bytes()).hexdigest(),protocol=a.protocol,fast=a.fast,slow=a.slow,slow_index=a.slow_index,reciprocal_scale=a.reciprocal_scale)
        (a.out/'TEMPORAL_INPUT.json').write_text(json.dumps(manifest,indent=2));sys.argv=[str(HERE/'live_tendon_trial.py'),'--protocol',a.protocol,'--physics','compliant','--source-task-basis','--config',str(a.config),'--out',str(a.out)]
        try:runpy.run_path(str(HERE/'live_tendon_trial.py'),run_name='__main__')
        finally:(a.out/'TEMPORAL_BINDINGS.json').write_text(json.dumps(records,indent=2))
