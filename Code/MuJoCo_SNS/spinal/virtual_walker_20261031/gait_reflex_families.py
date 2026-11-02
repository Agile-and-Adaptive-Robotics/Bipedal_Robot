"""Frozen-family extension of the validated selective Ia/Ib reflex experiment.

Each family is a separate process. The original syn6 wrapper/source and its
recorded experiment remain unchanged; no optimization or model copy is made.
"""
import os,sys,json,copy,time,argparse,contextlib,hashlib
from pathlib import Path
import numpy as np
import gait_spiking_reflex as X
F=X.F;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_family_reflex_v1')
PRISTINE=copy.deepcopy(F.P.G);FAMILY=None;OUT=None;LEDGER=None;CONFIG=None;PARAMETERS=None;SOURCE_PREFIX=None
COL={name:i for i,name in enumerate(X.COLS)}

def configure(family):
    global FAMILY,OUT,LEDGER,CONFIG,PARAMETERS,SOURCE_PREFIX
    FAMILY=family;OUT=ROOT/family
    LEDGER=X.REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'/(family+'_graph_ledger.json')
    base=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
    CONFIG=base/('rustuna_s3k_v1/BEST.json' if family=='s3k' else 'group_diagnosis/NORMALIZED_STANDING_CHECKPOINT.json')
    cfg=json.loads(CONFIG.read_text());PARAMETERS=cfg['trial']['params'] if family=='s3k' else cfg['global_parameters']
    SOURCE_PREFIX=base/('other_family_pressure/s3k_bio_prepared_v1' if family=='s3k' else 'normalized_prepared_transfer/w2lvar_bio_graded_v1')
    X.OUT=OUT;X.BEST=CONFIG;X.DIAGRAM=LEDGER;X.settings=lambda:PARAMETERS.copy()
    X.restore_parameters=restore;X.construct=construct

def restore():
    ledger=json.loads(LEDGER.read_text(encoding='utf-8'));F.P.G.clear();F.P.G.update(copy.deepcopy(PRISTINE))
    F.P.G['ii_set_pose']=1.;F.P.G['posture_to_mn']=0.
    if FAMILY=='w2lvar':
        import build_network_w2lvar as W
        F.P.G.update({k:float(v) for k,v in W.VARIANT_G.items()})
    assert F.P.G==ledger['global_G'],'Independent source defaults/explicit family overlay must match saved family metadata'
    return ledger

def scale_network(n,names):
    F.tie_weights(n);p=PARAMETERS
    for post,pre in zip(*np.nonzero(n.compiled.g_max_non)):
        if names[post].startswith('MN_') and names[pre].startswith('PF_'):n.compiled.g_max_non[post,pre]*=p['pf_mn_gain']
        if names[post].startswith('MN_') and names[pre].startswith('IBEXC'):n.compiled.g_max_non[post,pre]*=p['ib_group_gain']
    cfg=json.loads(CONFIG.read_text())
    if cfg.get('pathway_group_scales'):
        from muscle_map import classify
        for post,pre in zip(*np.nonzero(n.compiled.g_max_non)):
            if not names[post].startswith('MN_'):continue
            group=classify(names[post][3:]).groups[0];source=names[pre].split('_')[0]
            for key,value in cfg['pathway_group_scales'].items():
                cls,target=key.split(':')
                if source==cls and target in ('ALL',group):n.compiled.g_max_non[post,pre]*=value

def construct():
    ledger=restore();os.environ['AARL_NET']=FAMILY;m=F.mujoco.MjModel.from_xml_path(str(F.R.MODEL));acts=[m.actuator(i).name for i in range(m.nu)]
    n=F.B.build(acts,dt=F.P.DT,interleg=True);names=F.canonical_neuron_names(n);scale_network(n,names)
    for key,expected in ledger['compiled_matrix_sha256'].items():assert hashlib.sha256(np.asarray(getattr(n.compiled,key)).tobytes()).hexdigest()==expected,key
    return X.SparsePairwiseForward(n.compiled),names

def causal_audit():
    ledger=json.loads(LEDGER.read_text(encoding='utf-8'));edges=ledger['chemical_edges'];paths={}
    for source in X.CHANNELS:
        first=[e for e in edges if e['source']==source];two=[]
        for a in first:
            for b in edges:
                if b['source']==a['target'] and b['target']=='MN_soleus_r':two.append([a,b])
        direct=[e for e in first if e['target']=='MN_soleus_r']
        paths[source]=dict(direct_to_homonymous_MN=direct,two_hop_to_homonymous_MN=two,all_active_first_hops=first,has_active_direct_or_two_hop_MN_path=bool(direct or two))
    inactive=[e for e in ledger['disabled_source_pairs'] if e['target']=='MN_soleus_r' and e['source'].startswith(('PF_','IBEXC'))]
    return dict(family=FAMILY,settings=PARAMETERS,actual_neurons=ledger['census']['neurons'],spiking_cells=32,remaining_graded_cells=ledger['census']['neurons']-2,paths=paths,
      explicitly_zero_PF_or_IBEXC_MN_pairs=inactive,interpretation='Zero aggregate Ib excitation does not imply zero autogenic Ib inhibition. Active directed paths and their exact four synapse fields above are authoritative; no response is assumed from inactive edges.')

def prepare():
    OUT.mkdir(parents=True,exist_ok=True);audit=causal_audit();X.write(OUT/'CAUSAL_GRAPH_AUDIT.json',audit)
    X.prepare() # unchanged wrapper/default/complete-state audit, own construct hook
    protocol=json.loads((OUT/'PROTOCOL.json').read_text());protocol.update(family=FAMILY,standing_configuration=str(CONFIG),configuration_sha256=X.digest(CONFIG),
      remaining_graded_cells=audit['remaining_graded_cells'],driver=str(Path(__file__)),driver_sha256=X.digest(__file__),
      graph_audit=str(OUT/'CAUSAL_GRAPH_AUDIT.json'),group_scales=json.loads(CONFIG.read_text()).get('pathway_group_scales'),
      predictions='Use only the active paths in CAUSAL_GRAPH_AUDIT.json. Both frozen candidates have PF→MN=0 and IBEXC→MN=0; remaining direct/disynaptic Ib inhibition is retained. True zero responses are retained, never refitted.')
    if not list((OUT/'cases').glob('*.result.json')):X.write(OUT/'PROTOCOL.json',protocol)
    (OUT/'EXPERIMENT_RECORD.md').write_text('# '+FAMILY+' selective sensory-spiking reflex tests\n\nMachine: easteregg2. Fixed source candidate,32 LIF cells and'+str(audit['remaining_graded_cells'])+' graded cells. PF→MN and IBEXC→MN outputs remain exactly zero. Active autogenic reflex paths are audited before physical tests. No fitting, new support, or standing promotion.\n\nProtocol: `'+str(OUT/'PROTOCOL.json')+'`.\n',encoding='utf-8')
    print(json.dumps(dict(family=FAMILY,offline_ready=True,Ib_direct=len(audit['paths']['Ib_soleus_r']['direct_to_homonymous_MN']),Ib_two_hop=len(audit['paths']['Ib_soleus_r']['two_hop_to_homonymous_MN']))),flush=True)

def run_case(case,fresh=False):
    tag=case+('_fresh' if fresh else '');prefix=OUT/'cases'/(FAMILY+'_'+tag);reportfile=prefix.with_suffix('.comparison.json')
    if reportfile.exists():return json.loads(reportfile.read_text())
    protocol=json.loads((OUT/'PROTOCOL.json').read_text());assert (OUT/'OFFLINE_AUDIT.json').exists() and protocol['reflex_core_sha256']==X.core_digest()
    assert X.digest(CONFIG)==protocol['configuration_sha256'] and X.digest(LEDGER)==protocol['source_diagram_sha256']
    ledger=restore();objects=[];mode,pulse,off=X.CASES[case]
    def transform(inner,names):
        for key,expected in ledger['compiled_matrix_sha256'].items():assert hashlib.sha256(np.asarray(getattr(inner.c,key)).tobytes()).hexdigest()==expected,key
        obj=X.ReflexForward(inner,names,mode,pulse,off);objects.append(obj);return obj
    cfg=json.loads(CONFIG.read_text());p=PARAMETERS;start=time.perf_counter()
    with (OUT/'physics.log').open('a',encoding='utf-8') as log,contextlib.redirect_stdout(log):
        result=F.run(FAMILY,prefix.parent,duration=X.DURATION,tag=tag,gain=p['pf_mn_gain'],ib_group_gain=p['ib_group_gain'],standing_bias_scale=p['standing_bias_scale'],
          pathway_group_scales=cfg.get('pathway_group_scales'),pinned_actuator_settling=True,compiled_transform=transform,capture_neural_inputs=True)
    objects[0].save(prefix);assessment=X.score(result,prefix)[1]
    report=dict(family=FAMILY,case=tag,configuration_case=case,validation_repeat_of=case if fresh else None,result=str(prefix.with_suffix('.result.json')),reflex=str(prefix.with_suffix('.reflex.npz')),
      assessment=assessment,metadata=objects[0].metadata(),wall_s=time.perf_counter()-start,failure=result['failure'],protocol_sha256=X.digest(OUT/'PROTOCOL.json'),driver_sha256=X.digest(__file__))
    X.write(reportfile,report)
    with (OUT/'CASES.jsonl').open('a',encoding='utf-8') as stream:stream.write(json.dumps(report)+'\n')
    print(json.dumps(dict(family=FAMILY,case=tag,native_failure=result['failure'],completed_samples=result['completed_samples'],wall_s=report['wall_s'])),flush=True)
    return report

def prefix_audit():
    current=OUT/'cases'/(FAMILY+'_graded_baseline.npz');rows={}
    with np.load(SOURCE_PREFIX.with_suffix('.npz')) as old,np.load(current) as new:
        for key in ('t','qfull','q','ctrl','neuro','com'):
            a,b=old[key][:len(new[key])],new[key];rows[key]=dict(exact=bool(np.array_equal(a,b)),maximum_absolute_difference=float(np.max(abs(a-b))))
    report=dict(saved=str(SOURCE_PREFIX.with_suffix('.npz')),current=str(current),arrays=rows,exact=all(v['exact'] for v in rows.values()))
    X.write(OUT/'DEFAULT_PHYSICAL_PREFIX_AUDIT.json',report);assert report['exact'],'Do not continue: saved family physical prefix differs'
    print(json.dumps(dict(family=FAMILY,default_physical_prefix_exact=True)),flush=True)

def load(case):
    prefix=OUT/'cases'/(FAMILY+'_'+case);data={}
    for key,suffix in [('reflex','.reflex.npz'),('motor','.motor_state.npz'),('native','.npz'),('apparatus','.apparatus.npz')]:
        with np.load(prefix.with_suffix(suffix)) as a:
            names=('t','q','qfull','ctrl','neuro','com','key_joints') if key=='native' else a.files
            data[key]={n:a[n] for n in names}
    data['result']=json.loads(prefix.with_suffix('.result.json').read_text());return data

def analyze():
    cases={k:load(k) for k in X.CASES};rows=[];curves={};pairs=[]
    for mode in ('graded','spiking'):
        for channel in ('Ia','Ib'):pairs.append((mode+'_'+channel+'_pulse',mode+'_baseline',channel))
    for channel in ('Ia','Ib'):pairs.append(('spiking_'+channel+'_off_pulse','spiking_'+channel+'_off_baseline',channel))
    for key,baseline,channel in pairs:
        a,b=cases[key],cases[baseline];ra,rb=a['reflex']['rows'],b['reflex']['rows'];t=ra[:,0];pre=t<X.ONSET
        assert np.array_equal(ra[pre],rb[pre],equal_nan=True)
        for n in ('qfull','q','ctrl','neuro','com'):assert np.array_equal(a['native'][n][a['native']['t']<X.ONSET],b['native'][n][b['native']['t']<X.ONSET])
        dose=float(ra[:,COL['injected_current_nA']].sum()*F.P.DT);assert abs(dose-.1)<1e-12
        dm=ra[:,COL['soleus_MN_mV']]-rb[:,COL['soleus_MN_mV']];window=(t>=X.ONSET)&(t<X.ONSET+X.WIDTH+.25)
        midx=list(a['motor']['actuator_names']).index('soleus_r');df=-a['motor']['engine_force_cache_N'][:,midx]+b['motor']['engine_force_cache_N'][:,midx]
        tf=a['motor']['time']+F.R.WARMUP-F.P.DT;fw=(tf>=X.ONSET)&(tf<X.ONSET+X.WIDTH+.25)
        j=list(a['native']['key_joints']).index('ankle_angle_r');dq=a['native']['q'][:,j]-b['native']['q'][:,j]
        row=dict(family=FAMILY,case=key,matched_baseline=baseline,pre_onset_exact=True,dose_nA_s=dose,delta_MN_min_mV=float(dm[window].min()),delta_MN_max_mV=float(dm[window].max()),
          delta_tension_min_N=float(df[fw].min()),delta_tension_max_N=float(df[fw].max()),maximum_abs_ankle_effect_deg=float(abs(dq[window]).max()),source_outputs_off=X.CASES[key][2] is not None)
        if row['source_outputs_off']:
            exact=all(np.array_equal(a['native'][n],b['native'][n]) for n in ('qfull','q','ctrl','neuro','com'));assert exact and np.max(abs(dm))==0. and np.max(abs(df))==0.
            row['downstream_physics_exact']=True
        if X.CASES[key][0]=='spiking':
            r=a['reflex'];ew=(r['event_t']>=X.ONSET)&(r['event_t']<X.ONSET+X.WIDTH-1e-12);k=0 if channel=='Ia' else 1
            row['actual_spike_events_100ms']=int(r['event_spikes'][ew,k].sum());mask=r['event_spikes'].astype(bool)
            assert np.all(r['event_v'][mask]==0.) and np.all(r['event_ref'][mask]==.001)
            row['reset_refractory_exact']=True
        rows.append(row);curves[key]=dict(neural_time_s=t,delta_MN_mV=dm,force_time_neural_s=tf,delta_tension_N=df,delta_ankle_deg=dq)
    repeats={}
    for channel in ('Ia','Ib'):
        key='spiking_'+channel+'_pulse';anchor=OUT/'cases'/(FAMILY+'_'+key);fresh=OUT/'cases'/(FAMILY+'_'+key+'_fresh')
        if not fresh.with_suffix('.comparison.json').exists():continue
        checks={}
        for suffix in ('.npz','.apparatus.npz','.apparatus_extended.npz','.motor_state.npz','.inputs.npz','.reflex.npz','.release_state.npz'):
            with np.load(anchor.with_suffix(suffix)) as a,np.load(fresh.with_suffix(suffix)) as b:
                checks[suffix]={n:bool(np.array_equal(a[n],b[n],equal_nan=True)) if a[n].dtype.kind in 'fc' else bool(np.array_equal(a[n],b[n])) for n in a.files if n!='cfg'}
        exact=all(all(v.values()) for v in checks.values());assert exact;repeats[key]=dict(all_saved_arrays_exact=exact,checks=checks)
    conditions=[]
    for case,a in cases.items():
        r=a['result'];conditions.append(dict(case=case,failure=r['failure'],completed_free_samples=r['completed_samples'],mean_foot_load_bw=r['mean_ground_load_bw'],upward_support_bw=r['mean_upward_rig_support_bw'],downward_support_bw=r['mean_downward_rig_force_bw'],max_knee_angle_deg=r['max_knee_angle_deg']))
    X.write(OUT/'REFLEX_CHECKPOINT.json',dict(family=FAMILY,settings=PARAMETERS,paired_effects=rows,conditions=conditions,fresh_repeats=repeats,
      scope='Short2.4s free-physics functional reflexes, unchanged standing-duration scorer remains unpassed. Only32 sensory LIF cells; remaining source graph graded. No gain fitting, new support, or curriculum promotion.'))
    np.savez_compressed(OUT/'paired_pulse_effects.npz',**{key+'__'+n:v for key,data in curves.items() for n,v in data.items()})
    text='# '+FAMILY+' functional reflex results\n\nMachine: easteregg2. Frozen source candidate, no fitting. Exactly 32 sensory-spiking cells; PF→MN and IBEXC→MN remain 0. Accepted biological default prefix and all within-pair pre-onset traces are exact. Actual source-output-off controls eliminate downstream response while retaining spikes. Existing 0.6 s initialization precedes 2.4 s live physical integration; this does not satisfy the independent long-standing duration gate.\n\n'
    text+='| Case | MN delta (mV) | Tension delta (N) | Maximum ankle effect (deg) |\n| --- | ---: | ---: | ---: |\n'
    for r in rows:text+=f"| {r['case']} | {r['delta_MN_min_mV']:.6g} to {r['delta_MN_max_mV']:.6g} | {r['delta_tension_min_N']:.6g} to {r['delta_tension_max_N']:.6g} | {r['maximum_abs_ankle_effect_deg']:.6g} |\n"
    text+='\nAudit: `'+str(OUT/'CAUSAL_GRAPH_AUDIT.json')+'`. Exact saved physical prefix: `'+str(OUT/'DEFAULT_PHYSICAL_PREFIX_AUDIT.json')+'`. Actual event/MN/native-force/activation/angle traces are in cases; LIF voltage is threshold-normalized, decoded/graded voltages use relative-rest mV. No source model/controller files were changed. Fresh-repeat checks: '+', '.join(repeats)+'.\n'
    (OUT/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8');print(json.dumps(dict(family=FAMILY,paired_effects=rows,fresh_repeats=list(repeats))),flush=True)

def compare():
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none','pdf.fonttype':42})
    families=['s3k','w2lvar','syn6'];table=[];all_data={}
    for family in families:
        directory=ROOT/family if family!='syn6' else Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait2392_reflex_demo_v2')
        cp=json.loads((directory/'REFLEX_CHECKPOINT.json').read_text());ledger=json.loads((X.REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'/(family+'_graph_ledger.json')).read_text(encoding='utf-8'))
        for channel in ('Ia','Ib'):
            for mode in ('graded','spiking'):
                case=mode+'_'+channel+'_pulse';r=next(row for row in cp['paired_effects'] if row['case']==case)
                row=dict(family=family,channel=channel,encoding=mode,PF_MN_multiplier=ledger['settings']['pf_mn_gain'],IBEXC_MN_multiplier=ledger['settings']['ib_group_gain'])
                if family=='syn6':row.update(MN_delta_min_mV=r['delta_MN_mV_min'],MN_delta_max_mV=r['delta_MN_mV_max'],tension_delta_min_N=r['delta_tension_N_min'],tension_delta_max_N=r['delta_tension_N_max'],maximum_abs_ankle_effect_deg=r['maximum_abs_delta_ankle_deg'],actual_pulse_spikes=r.get('actual_pulsed_population_spikes_during_100ms'))
                else:row.update(MN_delta_min_mV=r['delta_MN_min_mV'],MN_delta_max_mV=r['delta_MN_max_mV'],tension_delta_min_N=r['delta_tension_min_N'],tension_delta_max_N=r['delta_tension_max_N'],maximum_abs_ankle_effect_deg=r['maximum_abs_ankle_effect_deg'],actual_pulse_spikes=r.get('actual_spike_events_100ms'))
                table.append(row)
        for case in ('graded_baseline','spiking_baseline','graded_Ia_pulse','graded_Ib_pulse','spiking_Ia_pulse','spiking_Ib_pulse'):
            prefix=directory/'cases'/(family+'_'+case);data={}
            for tag,suffix in [('reflex','.reflex.npz'),('motor','.motor_state.npz'),('native','.npz')]:
                with np.load(prefix.with_suffix(suffix)) as a:
                    selected=('t','q','key_joints') if tag=='native' else a.files
                    data[tag]={n:a[n] for n in selected}
            all_data[(family,case)]=data
    X.write(ROOT/'FAMILY_RESPONSE_TABLE.json',table)
    import csv
    with (ROOT/'FAMILY_RESPONSE_TABLE.csv').open('w',encoding='utf-8',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
    route_titles={'s3k':'s3k\nDirect Ib inhibition','w2lvar':'normalized w2lvar\nIb→IBIN→MN','syn6':'syn6\nIb inhibition +\ngroup excitation'}
    for channel in ('Ia','Ib'):
        fig,axes=plt.subplots(5,3,figsize=(6.5,8.4),sharex=True);handles=[]
        for j,family in enumerate(families):
            for mode,color,style in [('graded','#0072b2','-'),('spiking','#d55e00','--')]:
                a=all_data[(family,mode+'_'+channel+'_pulse')];b=all_data[(family,mode+'_baseline')]
                r=a['reflex'];t=r['rows'][:,0]-X.ONSET;dm=r['rows'][:,COL['soleus_MN_mV']]-b['reflex']['rows'][:,COL['soleus_MN_mV']]
                line,=axes[0,j].plot(t,r['rows'][:,COL['injected_current_nA']],color=color,ls=style,lw=1.,label=mode)
                if j==0:handles.append(line)
                axes[2,j].plot(t,dm,color=color,ls=style,lw=1.)
                aid=list(a['motor']['actuator_names']).index('soleus_r');df=-a['motor']['engine_force_cache_N'][:,aid]+b['motor']['engine_force_cache_N'][:,aid]
                tf=a['motor']['time']+F.R.WARMUP-F.P.DT-X.ONSET;axes[3,j].plot(tf,df,color=color,ls=style,lw=1.)
                jid=list(a['native']['key_joints']).index('ankle_angle_r');dq=a['native']['q'][:,jid]-b['native']['q'][:,jid];axes[4,j].plot(a['native']['t']-X.ONSET,dq,color=color,ls=style,lw=1.)
            r=all_data[(family,'spiking_'+channel+'_pulse')]['reflex'];k=0 if channel=='Ia' else 1;et=r['event_t']-X.ONSET
            for cell in range(16):
                fired=r['event_spikes'][:,k,cell]>0;axes[1,j].scatter(et[fired],np.full(fired.sum(),cell),marker='|',s=11,color='#222222')
            axes[1,j].set_yticks([0,15]);axes[0,j].set_title(route_titles[family] if channel=='Ib' else family+'\nIa homonymous excitation',fontsize=10)
            for row in (2,3,4):axes[row,j].axhline(0.,ls=':',color='#555555',lw=.8)
            for ax in axes[:,j]:ax.set_xlim(-.04,.32);ax.set_xticks([0.,.15,.3]);ax.grid(alpha=.18);ax.axvspan(0.,X.WIDTH,color='#aaaaaa',alpha=.12)
            axes[4,j].set_xlabel('Pulse-relative time (s)',fontsize=10)
        for i,label in enumerate(['Current\n(nA)','LIF events\n(cell index)','Soleus MN Δ\n(mV)','Tension Δ\n(N)','Ankle Δ\n(deg)']):axes[i,0].set_ylabel(label)
        fig.legend(handles,['graded','spiking'],loc='lower center',bbox_to_anchor=(.5,.135),ncol=2,frameon=False)
        fig.subplots_adjust(left=.14,right=.98,top=.9,bottom=.23,hspace=.22,wspace=.42)
        fig.text(.035,.025,channel+' input: +1 nA for 100 ms, two 16-cell sensory populations; other cells graded.\nΔ = pulse minus each family’s own no-pulse case. Output-off controls give zero\ndownstream response (dotted). These are different frozen standing candidates,\nnot a matched-parameter biological comparison or a standing/curriculum pass.',fontsize=10)
        for suffix in ('pdf','svg','png'):fig.savefig(ROOT/('family_'+channel+'_reflex_comparison.'+suffix),dpi=220,bbox_inches='tight',facecolor='white')
        plt.close(fig)
    text='# Frozen-family sensory-spiking reflex comparison\n\nMachine: easteregg2. See `FAMILY_RESPONSE_TABLE.csv/.json` for full numerical results and the per-family records for provenance, actual spikes, guards, and source-output-off controls. All source candidates use the same physical Gait2392/native-muscle framework, but different fixed network parameters and initialization-acquired states; this is not a pure topology-only effect estimate.\n\n'
    text+='| Family | Input | Encoding | MN delta range (mV) | Tension delta range (N) | Ankle effect (deg) |\n| --- | --- | --- | ---: | ---: | ---: |\n'
    for r in table:text+=f"| {r['family']} | {r['channel']} | {r['encoding']} | {r['MN_delta_min_mV']:.5g} to {r['MN_delta_max_mV']:.5g} | {r['tension_delta_min_N']:.5g} to {r['tension_delta_max_N']:.5g} | {r['maximum_abs_ankle_effect_deg']:.5g} |\n"
    text+='\nS3k uses direct inhibitory Ib→MN as a model abstraction; normalized w2lvar uses excitatory Ib→IBIN followed by inhibitory IBIN→MN. Both retain zero PF→MN and aggregate IBEXC→MN outputs. Syn6 retains its fitted nonzero group-excitatory branch, alongside inhibition. A rebound is observed in syn6, but individual inhibitory/excitatory branch ablations were not performed here, so its unique cause is not assigned. No multiplicative stance AND gate is invented. These short 3 s neural-time/2.4 s free-physics trials retain the original long-standing scorer threshold and do not promote walking or reduced support. Native guards and actual support reactions remain recorded.\n\n'
    for family in ('s3k','w2lvar'):text+='- `'+str(ROOT/family/'EXPERIMENT_RECORD.md')+'`\n'
    text+='\nProgram on easteregg2: `'+str(Path(__file__))+'`. Use `D:\\Anaconda\\envs\\myo\\python.exe`, then this full program path and `--family s3k --sequence` or `--family w2lvar --sequence`. Completed cases are reused. `--family NAME --case spiking_Ia_pulse --fresh` and the Ib counterpart repeat the identical condition; `--family NAME --analyze` updates checks; `--compare` makes these figures. No original source/controller/model, GUI, or manuscript file changed.\n'
    (ROOT/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8');print(json.dumps(dict(comparison_rows=len(table),folder=str(ROOT))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--family',choices=['s3k','w2lvar']);p.add_argument('--prepare',action='store_true');p.add_argument('--sequence',action='store_true');p.add_argument('--case',choices=X.CASES);p.add_argument('--fresh',action='store_true');p.add_argument('--analyze',action='store_true');p.add_argument('--compare',action='store_true');a=p.parse_args()
    if a.compare:compare();sys.exit(0)
    if not a.family:p.error('--family is required for model operations')
    configure(a.family)
    if a.prepare:prepare()
    elif a.case:run_case(a.case,a.fresh)
    elif a.sequence:
        for case in X.CASES:
            run_case(case)
            if case=='graded_baseline':prefix_audit()
    elif a.analyze:analyze()
    else:p.error('Select a stage')
