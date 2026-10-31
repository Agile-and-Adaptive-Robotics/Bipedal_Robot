from source_inputs import historical_source
"""Bounded tonic-current A/B with existing source phase-init pulses, no mechanics."""
import json,math
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from import_virtual_walker import HERE,REPO,build
from native_kernel import NativeController
outdir=HERE/'runs/newer_rg_isolation_v1';outdir.mkdir(exist_ok=True)
controllers={}
for kind,stem in [('GUI1nA','W2L_VirtualWalker_Newer_GUI1nA'),('20nA','W2L_VirtualWalker_Newer_SIContact')]:
 source=historical_source(stem);path=build(False,source);meta=json.loads((path.parent/'IMPORT_AUDIT.json').read_text());controllers[kind]=(source,meta['muscles'],NativeController(source,meta['muscles'],neural_time_scale=2.))
a=controllers['GUI1nA'][2];b=controllers['20nA'][2];assert np.array_equal(a.params,b.params) and np.array_equal(a.edges,b.edges) and np.array_equal(a.tonic,b.tonic) and a.stimuli==b.stimuli
rows=[]
for label,tonic in [('source_2_0',None),('tied_0_0',0.),('tied_0_0_V3_output_off',0.),('source_2_0_V3_output_off',None),('original80_positive_control',None)]:
 source,bank,_=controllers['GUI1nA']
 if label=='original80_positive_control':source=REPO/'Neuromechanical_Models/W2L_VirtualWalker_20261031/W2L_VirtualWalker_Standalone.asim'
 c=NativeController(source,bank,symmetric_rg_ext_tonic=tonic,neural_time_scale=2.);edited=[]
 if 'V3_output_off' in label:
  for k,e in enumerate(c.edges):
   if (c.names[int(e[0])],c.names[int(e[1])]) in [('L V3','R RG ext IN'),('R V3','L RG ext IN')]:edited.append(dict(source=c.names[int(e[0])],target=c.names[int(e[1])],source_gmax=float(e[2]),new_gmax=0.));c.edges[k,2]=0.
  assert len(edited)==2
 dt=.002;times=np.arange(0,8.,dt);v=np.empty((len(times),len(c.v)));terms=[];rg=[c.names.index(n) for n in ('L RG ext','L RG flx','R RG ext','R RG flx') if n in c.names];zeros=np.zeros(len(bank))
 assert np.all(c.params[:,2]==1.),'Current-term nA label below assumes source Size1'
 for k,t in enumerate(times):c.step(t,zeros,zeros,zeros,{},dt,feedback=False);v[k]=c.v
 # At the final state, physical driving terms—not the rest-referenced exponential numerator.
 for post in rg+[c.names.index(n) for n in ('L RG ext IN','R RG ext IN','L V3','R V3') if n in c.names]:
  incoming=[]
  for edge in c.edges:
   pre,target=int(edge[0]),int(edge[1])
   if target!=post:continue
   conductance=edge[2]*np.clip((c.v[pre]-edge[3])/(edge[4]-edge[3]),0,1);current=float(conductance*(edge[5]-c.v[post]));incoming.append(dict(source=c.names[pre],source_id=next(g for g,i in c.ids.items() if i==pre),g_active= float(conductance),reversal_mV=float(edge[5]),current_nA=current))
  leak=float(c.params[post,0]-c.v[post]);ca=float(c.params[post,3]*c.m[post]*c.h[post]*(c.params[post,10]-c.v[post]));tonic_current=float(c.tonic[post]);syn=sum(x['current_nA'] for x in incoming);terms.append(dict(target=c.names[post],voltage_mV=float(c.v[post]),leak_nA=leak,calcium_nA=ca,tonic_nA=tonic_current,afferent_nA=0.,synaptic_nA=syn,net_nA=leak+ca+tonic_current+syn,incoming=sorted(incoming,key=lambda x:abs(x['current_nA']),reverse=True)))
 metrics={}
 for i in rg:
  mask=times>=1.;trace=v[mask,i];peaks,_=find_peaks(trace,prominence=1.,distance=round(.15/dt));metrics[c.names[i]]=dict(range_mV=float(np.ptp(trace)),prominent_peaks=int(len(peaks)),median_period_s=float(np.median(np.diff(times[mask][peaks]))) if len(peaks)>1 else None)
 stimuli=[dict(target=c.names[i],start_s=s,end_s=e,amplitude_nA=a) for i,s,e,a in c.stimuli];record=dict(label=label,source=str(source),census=[len(c.v),len(c.edges)],tonic_RGext_nA=[float(c.tonic[c.names.index(s+' RG ext')]) for s in ('L','R') if s+' RG ext' in c.names],source_existing_startup_stimuli=stimuli,extra_pulse_added=False,neural_time_scale=2.,edited_edges=edited,mechanics='None; isolated neural core with all physical afferents disabled',GUI1nA_and_20nA_candidates_have_identical_intrinsic_core=True,metrics=metrics,oscillation_gate=all(m['prominent_peaks']>=4 and m['range_mV']>=1 for m in metrics.values()),final_current_terms=terms,current_term_scope='Native-equivalent nA for these Size1 neurons, evaluated at final voltage; calcium gating is final state. All startup pulses are source-loaded, not duplicated.');np.savez_compressed(outdir/(label+'.npz'),t=times,neuro=v,neuron_names=c.names);rows.append(record);(outdir/'RESULT.json').write_text(json.dumps(rows,indent=2))
print(json.dumps([{k:r[k] for k in ('label','source_existing_startup_stimuli','extra_pulse_added','metrics','oscillation_gate')} for r in rows],indent=2))
