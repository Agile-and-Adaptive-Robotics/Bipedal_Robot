"""Bounded electrical-only heel phase probes of an explicit composed profile.

No parameter fit, sensory oracle, body binding, or mechanical simulation.
"""
import os
for key in('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib,time,sys
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
import gait_design_composition as COM,graded_profile as GP,central_afferent_budget as CB

def crossings(t,y,direction=1):
    y=direction*y;k=np.flatnonzero((y[:-1]<=0)&(y[1:]>0))
    return t[k]+(t[k+1]-t[k])*(-y[k])/(y[k+1]-y[k])

def plot(t,base,pulse,names,input_current,onset,width,period,path,label,pf_names):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    first=max(0.,onset-period);last=min(t[-1],onset+2.3*period);use=(t>=first)&(t<=last)
    pairs=[('RG_E_r','RG_F_r'),('InE_r','InF_r'),tuple(pf_names),('MN_iliacus_r','MN_bifemlh_r'),('MN_vas_int_r','MN_bifemsh_r'),('MN_soleus_r','MN_tib_ant_r')]
    titles=['RG','RG inhibitory relays','Representative synergy PF','Hip MNs','Knee MNs','Ankle MNs']
    fig,axes=plt.subplots(7,1,figsize=(6.5,10),sharex=True)
    axes[0].plot(t[use],input_current[use],color='black');axes[0].set_ylabel('Heel nA');axes[0].set_title(label,loc='left',fontsize=10)
    for ax,pair,title in zip(axes[1:],pairs,titles):
        for j,name in enumerate(pair):
            i=names.index(name)
            ax.plot(t[use],base[use,i],color='.65',ls=('-','--')[j],lw=1.)
            ax.plot(t[use],pulse[use,i],color=('#0072B2','#D55E00')[j],lw=1.,label=name.removesuffix('_r'))
        ax.set_ylabel('Model mV');ax.set_title(title,loc='left',fontsize=10);ax.legend(fontsize=10,ncol=2,loc='upper right');ax.grid(alpha=.2)
    for ax in axes:ax.axvspan(onset,onset+width,color='#f2c27d',alpha=.3);ax.set_xlim(first,last)
    axes[-1].set_xlabel('Neural output time (s); gray: matched unperturbed baseline')
    fig.tight_layout();fig.savefig(path.with_suffix('.png'),dpi=180);fig.savefig(path.with_suffix('.pdf'));plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--profile',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--duration',type=float,default=12.);p.add_argument('--mono-overlay',action='store_true');a=p.parse_args()
    if not 12.<=a.duration<=24.:raise ValueError('Bounded observation duration must be12–24s')
    a.out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();profile=json.loads(a.profile.read_text())
    registration={'profile':str(a.profile.resolve()),'profile_sha256':CB.sha(a.profile),'duration_s':a.duration,'pulse_width_s':.1,'normalized_heel_signal':1.,'phases':[.2,.5,.8],
                  'baseline':'Exact source-recipe steady WALK DRIVE and POSTURE0.6; physical afferent/contact/directPOST/balance inputs zero. This is a declared unbound deafferented electrical fixture, not a body replay.',
                  'phase_definition':'RG-E minus RG-F upward zero crossing is phase0; period from late baseline crossings. Negative next-event time shift means advance.',
                  'policy':'Only HEEL_c_r differs in pulse branches. Existing heel PF branch remains; no invented contralateral reset, parameter change, or initial-state change.',
                  'execution':'Exact sparse compiled forward; distinct input rows prevalidated, full composed-profile checks before/after each case. No per-step full-matrix hashing; equations, step, and states unchanged.',
                  'mono_overlay':a.mono_overlay,'training':False,'physics_steps':0}
    (a.out/'PREREGISTRATION.json').write_text(json.dumps(registration,indent=2),encoding='utf-8')
    if a.mono_overlay:
        import syn6_mono_composed as builder
    else:builder=COM
    net,proof=builder.build_named_profile(a.profile);names=GP.population_names(net);w=net.compiled;initial=CB.clone(w);dt=float(CB.core(net).dt)
    base_profile=json.loads(Path(profile['base_profile']['path']).read_text())if a.mono_overlay else profile
    recipe_path=HERE.parent/'curriculum_syn6_stage5.json'
    if base_profile['source_sha256'][str(recipe_path)]!=CB.sha(recipe_path):raise ValueError('Source WALK recipe no longer bound to selected base profile')
    recipe=json.loads(recipe_path.read_text())['params']
    def verify_case():
        if a.mono_overlay:return builder.verify(net)
        return COM.verify(net,proof['composition']['coincidence_design'])
    u=net.make_inputs();u[net.input_index('DRIVE')]=recipe['drive'];u[net.input_index('POSTURE')]=.6
    u,domain=GP.adapt_inputs(net,u,enforce_domain=True);heel_port=net.input_index('HEEL_c_r');amplitude=float(recipe['heel_rge'])
    upulse=u.copy();upulse[heel_port]=amplitude;upulse,pulse_domain=GP.adapt_inputs(net,upulse,enforce_domain=True)
    if a.mono_overlay:u=builder.validate_inputs(net,u);upulse=builder.validate_inputs(net,upulse)
    heel=net.idx['HEEL_r'];ine=net.idx['InE_r'];inf=net.idx['InF_r'];core=CB.core(net)
    assert core.g_max_non[inf,heel]==0 and core.g_max_non[ine,heel]>0
    parameter_digest=GP.field_digest(core);N=round(a.duration/dt);t=(np.arange(N)+1)*dt
    allbase=np.empty((N,len(names)));parameter_fields=['gmax','reversal','elo','erange'];operators={k:getattr(w,k).copy()for k in parameter_fields}
    verify_case()
    for k in range(N):
        net.compiled.forward(u);allbase[k]=net.compiled.V
        if not np.isfinite(allbase[k]).all():raise ValueError('Nonfinite baseline neural state')
    verify_case()
    diff=allbase[:,net.idx['RG_E_r']]-allbase[:,net.idx['RG_F_r']];events=crossings(t,diff);stable=events[events>=3.]
    report={'registration':registration,'source_recipe_drive_nA':float(recipe['drive']),'posture_input_nA':.6,'heel_pulse_nA':amplitude,'source_heel_normalization':'runner heel_sig=min(heel_force/(0.35BW),1.5); this probe chooses1.0, within that caller range.',
            'afferent_domain_check':domain,'pulse_afferent_domain_check':pulse_domain,'nodes':len(names),'dt_s':dt,'baseline_RG_E_crossings_s':events.tolist(),'baseline_all_neurons_finite':True,'baseline_all_neuron_voltage_range_mV':[float(allbase.min()),float(allbase.max())],
            'heel_to_InF_gmax_uS':0.,'heel_to_InE_gmax_uS':float(core.g_max_non[ine,heel]),'phase_responses':[],'body_result':False,'trained':False}
    save_names=['RG_E_r','RG_F_r','InE_r','InF_r','HEEL_r']+[f'PF_S{i}_r'for i in range(1,7)]+['MN_iliacus_r','MN_bifemlh_r','MN_vas_int_r','MN_bifemsh_r','MN_soleus_r','MN_tib_ant_r']
    ids=np.array([net.idx[n]for n in save_names]);base=allbase[:,ids]
    np.savez_compressed(a.out/'BASELINE.npz',t=t,voltage=base,neuron_names=np.array(save_names),constant_input=u,input_names=np.array(net.inputs))
    if len(stable)<4:
        report['status']='No qualified autonomous baseline: fewer than4 late RG-E transitions. No phase response inferred.'
    else:
        period=float(np.median(np.diff(stable[-5:])));anchor=float(stable[stable>=4.][0]);report['baseline_period_s']=period
        pf_names=[f'PF_S{net.syn_families["r"].index(family)+1}_r'for family in('E','F')]
        for phase in registration['phases']:
            onset=anchor+phase*period;first=round(onset/dt);last=first+round(.1/dt);onset=first*dt
            net.compiled=CB.clone(initial);observed=np.empty((N,len(ids)));pre_exact=True;state_finite=True;heel_in_current=np.zeros(N)
            verify_case()
            for k in range(N):
                active=first<=k<last;current=upulse if active else u;current_w=net.compiled
                before=current_w.c.V;slots=np.flatnonzero((current_w.post==ine)&(current_w.pre==heel))
                slot=int(slots[0]);fraction=np.clip((before[heel]-current_w.elo[slot])/current_w.erange[slot],0,1)
                heel_in_current[k]=current_w.gmax[slot]*fraction*(current_w.reversal[slot]-before[ine])
                current_w.forward(current);values=current_w.V;observed[k]=values[ids]
                if k<first:pre_exact=pre_exact and np.array_equal(values,allbase[k])
                state_finite=state_finite and bool(np.isfinite(values).all())
            verify_case()
            pert_diff=observed[:,save_names.index('RG_E_r')]-observed[:,save_names.index('RG_F_r')]
            shifts={}
            for direction,name in((1,'E_onset'),(-1,'F_onset')):
                be=crossings(t,diff,direction);pe=crossings(t,pert_diff,direction);ordinal=np.flatnonzero(be>onset+.1)
                matched=[(float(be[j]),float(pe[j]),float(pe[j]-be[j]))for j in ordinal[:3]if j<len(pe)]
                shifts[name]=dict(matched_events_baseline_pulse_shift_s=matched,next_event_phase_advance_cycles=-matched[0][2]/period if matched else None)
            pulse_current=np.zeros(N);pulse_current[first:last]=amplitude
            prefix=f'HEEL_PHASE_{str(phase).replace(".","p")}'
            np.savez_compressed(a.out/(prefix+'.npz'),t=t,voltage=observed,neuron_names=np.array(save_names),heel_applied_current_nA=pulse_current,heel_to_InE_actual_chemical_nA=heel_in_current)
            row=dict(phase=phase,onset_s=onset,width_s=.1,current_nA=amplitude,all_other_inputs_and_initial_state_identical=True,all_neural_states_exact_before_pulse=pre_exact,all_neurons_finite=state_finite,transitions=shifts,
                     peak_heel_InE_current_nA=float(heel_in_current.max()),maximum_heel_InE_current_nA=float(heel_in_current.max()),minimum_heel_InE_current_nA=float(heel_in_current.min()))
            report['phase_responses'].append(row)
            plot(t,base,observed,save_names,pulse_current,onset,.1,period,a.out/prefix,f'Syn6 single-InE heel probe | baseline phase{phase:.1f}',pf_names)
        report['status']='Three finite-pulse electrical phase probes completed; not a full PRC or a neuromechanical result.'
    report['parameters_unchanged']=GP.field_digest(CB.core(net))==parameter_digest
    report['sparse_operators_unchanged']=all(np.array_equal(getattr(net.compiled,k),value)for k,value in operators.items())
    if not report['parameters_unchanged']or not report['sparse_operators_unchanged']:raise AssertionError('Phase observation changed model parameters')
    report['wall_s']=time.perf_counter()-start
    (a.out/'PHASE_RESPONSE.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items()if k in('status','baseline_period_s','phase_responses','parameters_unchanged','nodes','wall_s')},indent=2))
if __name__=='__main__':main()
