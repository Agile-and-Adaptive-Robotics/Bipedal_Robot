"""Saved-array and primary-record audit; no neural/mechanical integration."""
from pathlib import Path
import json,hashlib,sys,importlib.util,math
import numpy as np
from scipy.signal import find_peaks
ROOT=next(p for p in Path(__file__).resolve().parents if p.name=='Bipedal_Robot')
P=ROOT/'Code/MuJoCo_SNS/spinal'
H=P/'reports_spiking_20261002'
OUT=Path(__file__).parent
hashes={}
def record(p):
    p=Path(p);hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return p
result={'scope':'Read-only saved arrays/records; no simulations or generator reruns.','raw':{}}
spec=importlib.util.spec_from_file_location('historical_kine_readonly',P/'kine_ref.py');kine=importlib.util.module_from_spec(spec);spec.loader.exec_module(kine)
record(P/'kine_ref.py')
z=np.load(record(H/'goal4_gate2_traces.npz'),allow_pickle=False)
t=z['t'];dt=float(np.median(np.diff(t)))
v=z['v_ns'];pks,_=find_peaks(v,prominence=.3);ns_period=float(np.mean(np.diff(t[pks])))
win=int(round(.1/dt));env=np.convolve(np.diff(z['cnt_e'],prepend=0)-np.diff(z['cnt_f'],prepend=0),np.ones(win)/win,mode='same')*1000*dt/.1
pks,_=find_peaks(env,prominence=1,distance=int(.15/dt));sp_period=float(np.mean(np.diff(t[pks])))
result['rhythm']={'source':str(H/'goal4_gate2_traces.npz'),'window_s':[float(t[0]),float(t[-1])],'ns_period_s':ns_period,'spiking_period_s':sp_period,'period_increase_percent':100*(sp_period-ns_period)/ns_period,'definition':'Whole20s interpeak mean; graded E-F voltage prominence0.3mV; spiking E-F count envelope100ms, prominence1 and minimum peak distance150ms; same scored filters as source gate.'}
for label in ['s3k','spiking']:
    fp=record(P/f'spinal_run_spkbase_{label}.npz')
    z=np.load(fp,allow_pickle=True)
    rec={'keys':{k:list(z[k].shape)for k in z.files},'time_s':[float(z['t'][0]),float(z['t'][-1])]}
    names=list(z['key_joints']);q=z['q'];rec['joint_names']=names
    for k in ['tilt_max','bal_tilt_max','nan','fell','walk_start','key_joints']:
        if k in z.files and k!='key_joints':rec[k]=np.asarray(z[k]).tolist()
    if 'pelvis_tilt'in names:rec['maximum_abs_pelvis_tilt_deg']=float(np.max(abs(q[:,names.index('pelvis_tilt')])))
    if 'contact'in z.files:
        rec['contact_shape']=list(z['contact'].shape)
        for start in [2.,4.5,5.]:
            w=z['t']>=start
            rec[f'contact_fraction_after_{start}s_gt20N']=np.mean(z['contact'][w]>20,axis=0).tolist()
    rec['recomputed_cycles']={}
    for side in ['r','l']:
        cyc=kine.sim_side(z['t'],z['q'],z['neuro'],2.,side,z['contact'])
        rec['recomputed_cycles'][side]={'cycles':int(cyc[5]),'contact_fraction':float(cyc[3]),'used_contact':bool(cyc[4]),'mean_cycle_knee_min_deg':float(np.min(cyc[0]['knee']))if cyc[0]is not None else None}
    try:
        cmp=kine.compare(z['t'],z['q'],z['neuro'],2.,contact=z['contact'])
        rec['current_source_rescore']={k:v for k,v in cmp.items()if isinstance(v,(float,int,bool))}
        record(kine.IK_MOT);record(kine.GRF_MOT)
    except Exception as exc:rec['rescore_error']=repr(exc)
    result['raw'][label]=rec
for n in ['goal1_walkers_spiking.md','goal3_animatlab_spiking.md','goal4_dissertation_spiking.md','tools/check_selfsustain_spiking.py','tools/goal4_dump_gate2_traces.py','tools/baseline_eval.py','tools/analyze.py','tools/make_spiking.py','tools/sameleg_xcorr.py','logs/baseline_s3k.log','logs/ground_eval_spiking.log','logs/gate2_rhythm2.log']:
    record(H/n)
result['animatlab_metrics']={}
for name in ['metrics_baseline_W2L_modern.json','metrics_baseline_BipRG_modern.json','metrics_W2L_modern_spiking.json','metrics_BipRG_modern_spiking.json','metrics_BipRG_modern_spiking_v2.json','metrics_Biped_standalone_stale_spiking.json','metrics_bisect_W2L_neuronsonly.json']:
    d=json.loads(record(H/name).read_text());result['animatlab_metrics'][name]=d
    asim=Path(d['asim'])
    if asim.exists():record(asim)
    d['independent_raw_joint_ranges']={}
    d['independent_raw_rg_metrics']={}
    for ch,entry in d.get('charts',{}).items():
        if not any(c.get('type')=='angle_deg'for c in entry.get('cols',[])):continue
        fp=asim.parent/entry['file']
        if not fp.exists():d['independent_raw_joint_ranges'][ch]={'missing':str(fp)};continue
        record(fp)
        raw=np.genfromtxt(fp,delimiter='\t',skip_header=1)
        while len(raw)and np.all(raw[-1,2:]==0):raw=raw[:-1]
        vals={}
        for i,c in enumerate(entry['cols'],start=2):
            if c.get('type')!='angle_deg':continue
            value=raw[:,i]
            if not c['name'].lower().endswith('_deg'):value=np.rad2deg(value)
            vals[c['name']]={'min':float(np.min(value)),'max':float(np.max(value)),'sample_window_s':[float(raw[0,1]),float(raw[-1,1])],'matches_saved_rounding':bool(round(float(np.min(value)),1)==c['min']and round(float(np.max(value)),1)==c['max'])}
        d['independent_raw_joint_ranges'][ch]=vals
    for ch,entry in d.get('charts',{}).items():
        if not any('RG'in c.get('name','')and'IN'not in c.get('name','')and'mean_mV'in c for c in entry.get('cols',[])):continue
        fp=asim.parent/entry['file']
        if not fp.exists():continue
        record(fp);raw=np.genfromtxt(fp,delimiter='\t',skip_header=1)
        while len(raw)and np.all(raw[-1,2:]==0):raw=raw[:-1]
        t=raw[:,1];times={}
        for i,c in enumerate(entry['cols'],start=2):
            cell=c.get('name','')
            if 'RG'not in cell or 'IN'in cell or 'mean_mV'not in c:continue
            v=raw[:,i];threshold=-.04 if d['mode']=='spiking'else-.06
            ix=np.flatnonzero((v[1:]>threshold)&(v[:-1]<=threshold))+1
            event=t[ix-1]+(threshold-v[ix-1])*(t[ix]-t[ix-1])/(v[ix]-v[ix-1])
            kept=[]
            for x in event:
                if not kept or x-kept[-1]>= (.001 if d['mode']=='spiking'else.05):kept.append(float(x))
            burst=kept if d['mode']=='graded'else([kept[0]]+[b for a,b in zip(kept[:-1],kept[1:])if b-a>.05]if kept else[])
            period=sorted(np.diff(burst))[len(burst[1:])//2]if len(burst)>=2 else None
            measured={'crossings_or_spikes':len(kept),'burst_count':len(burst),'median_period_s':None if period is None else float(period),'source_window_s':[float(t[0]),float(t[-1])]}
            if d['mode']=='spiking':measured['rate_Hz']=float(len(kept)/(t[-1]-t[0]))
            d['independent_raw_rg_metrics'][ch+'/'+cell]=measured;times[cell]=kept
        if 'same_leg_xcorr'in d and 'L RG ext'in times and'L RG flx'in times:
            count=[]
            for cell in ['L RG ext','L RG flx']:
                hist=np.zeros(int((t[-1]-t[0])/.02))
                for x in times[cell]:
                    j=int((x-t[0])/.02)
                    if 0<=j<len(hist):hist[j]+=1
                count.append(hist)
            d['independent_raw_rg_same_leg_zero_lag_correlation']=float(np.corrcoef(*count)[0,1])
result['hashes']=hashes
(OUT/'HISTORICAL_SPIKING_RESULTS.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({'rhythm':result['rhythm'],'raw':{k:{a:b for a,b in v.items()if a not in ['keys','current_source_rescore']}for k,v in result['raw'].items()}},indent=2))
for name,d in result['animatlab_metrics'].items():
    print(name,json.dumps({k:v for k,v in d.items()if k in ['rg_bursts','antiphase','same_leg_xcorr','independent_raw_joint_ranges']},indent=1))
    for ch,c in d.get('charts',{}).items():
        vals=[v for v in c.get('cols',[])if ('RG' in v.get('name','')and 'rate_Hz'in v)]
        if vals:print('RG RATES',vals)
