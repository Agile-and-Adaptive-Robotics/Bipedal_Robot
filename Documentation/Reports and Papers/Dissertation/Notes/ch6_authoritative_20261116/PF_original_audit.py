from pathlib import Path
import json, hashlib, sys, importlib.util
import numpy as np

root=Path(r'D:\GitHub\Bipedal_Robot')
spinal=root/'Code/MuJoCo_SNS/spinal'
sys.path.insert(0,str(spinal))
p=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104\source_snapshot_pf_design_20261105\pf_timing_audit.py')
spec=importlib.util.spec_from_file_location('old_audit',p)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
fsa_path=spinal/'fsa_results/fsa_backsolve.npz'
fsa=np.load(fsa_path)
raw=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104')
portable=root/'Documentation/Reports and Papers/Dissertation/Defense/Session1_Handoff_20261104/experiments_20261104/circuit_trace_20261105'
cases=[('s3k',raw/'s3k_relay_output_20261104'),('syn6',raw/'waveform_edge_D_20261104/syn6')]
out={'method_file':str(p),'method_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'fsa_sha256':hashlib.sha256(fsa_path.read_bytes()).hexdigest(),'cases':{}}
for family,fit in cases:
    old=json.loads((portable/family/'timing_audit/PF_TIMING_AUDIT.json').read_text())
    T=old['reference_period_s']; z=np.load(fit/'BEST_TRACE.npz'); names=z['pf_names'].tolist()
    # Same declared original clock, 200 phase samples for each of cycles 2..10.
    phase=np.arange(200)/200; ts=phase*T
    cf=Path(json.loads((fit/'MANIFEST.json').read_text())['fit'])
    cfg=json.loads((cf/'MANIFEST.json').read_text())['configuration']
    rows=[];rank={};max_delta=0
    for side in ('r','l'):
        ys=[];wanted=[]
        for col,name in enumerate(names):
            if not name.endswith('_'+side):continue
            cy=np.array([np.interp(ts+n*T,z['t'],np.clip(z['pf'][:,col]/5,0,1)) for n in range(2,11)])
            y=cy.mean(0); af=mod.features(y); row={'cell':name,'actual':af};ys.append(y)
            if family=='syn6':
                target=mod.primitive_target(fsa,cfg,side,ts,int(name.split('_')[1][1:])-1)
                tf=mod.features(target); wanted.append(target)
                lag=(tf['fundamental_phase_deg']-af['fundamental_phase_deg']+180)%360-180
                count=len(tf['onset_phase_percent'])==len(af['onset_phase_percent']) and len(tf['offset_phase_percent'])==len(af['offset_phase_percent'])
                row.update(target=tf,lag=lag,passed=bool(abs(lag)<=30 and abs(af['duty_fraction']-tf['duty_fraction'])<=.15 and count and af['excursion']>=.2))
            oldrow=next(r for r in old['rows'] if r['cell']==name)
            for key in ('minimum','maximum','excursion','duty_fraction','fundamental_phase_deg'):
                max_delta=max(max_delta,abs(af[key]-oldrow['actual'][key]))
            rows.append(row)
        rank[side]={'actual':mod.effective_rank(np.array(ys).T)}
        if wanted:rank[side]['target']=mod.effective_rank(np.array(wanted).T)
    out['cases'][family]={'raw':str(fit/'BEST_TRACE.npz'),'raw_sha256':hashlib.sha256((fit/'BEST_TRACE.npz').read_bytes()).hexdigest(),'window_seconds':[2*T,(10+199/200)*T], 'T':T,'rows':rows,'rank':rank,'max_feature_delta_from_retained_audit':max_delta}

# Recompute repaired trace under same declaration without re-running neural model.
fit=raw/'syn6_PF_input_repair_mirrored_20261105'
cfg=json.loads((fit/'MANIFEST.json').read_text())['configuration'];T=cfg['target_period_s'];z=np.load(fit/'BEST_TRACE.npz')
rows=[]
for col,name in enumerate(z['pf_names'].tolist()):
    phase=np.arange(200)/200
    y=np.array([np.interp((phase+n)*T,z['t'],np.clip(z['pf'][:,col]/5,0,1)) for n in range(2,11)]).mean(0)
    target=mod.primitive_target(fsa,cfg,name[-1],phase*T,int(name.split('_')[1][1:])-1)
    af,tf=mod.features(y),mod.features(target)
    lag=(tf['fundamental_phase_deg']-af['fundamental_phase_deg']+180)%360-180
    count=len(tf['onset_phase_percent'])==len(af['onset_phase_percent']) and len(tf['offset_phase_percent'])==len(af['offset_phase_percent'])
    rows.append({'cell':name,'lag':lag,'duty_error':abs(af['duty_fraction']-tf['duty_fraction']),'passed':bool(abs(lag)<=30 and abs(af['duty_fraction']-tf['duty_fraction'])<=.15 and count and af['excursion']>=.2)})
out['repaired']={'raw':str(fit/'BEST_TRACE.npz'),'raw_sha256':hashlib.sha256((fit/'BEST_TRACE.npz').read_bytes()).hexdigest(),'passed_count':sum(r['passed'] for r in rows),'rows':rows}
destination=Path(__file__).resolve().parent/'PF_ORIGINAL_RESULTS.json'
destination.write_text(json.dumps(out,indent=2),encoding='utf-8')
print('Saved',destination)
