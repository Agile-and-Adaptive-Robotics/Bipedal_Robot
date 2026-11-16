"""Read-only recomputation of saved historical human-model IK comparisons."""
from pathlib import Path
import hashlib
import importlib.util
import json
import time
import numpy as np

ROOT=Path(r'D:\GitHub\Bipedal_Robot')
S=ROOT/'Code/MuJoCo_SNS/spinal'
O=S/'reports_20261028_publication_flow/human_ik'
OUT=Path(__file__).resolve().parent
sources={}
def bind(p):
    p=Path(p);sources[str(p)]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};return p
start=time.perf_counter()
kpath=bind(S/'kine_ref.py');spec=importlib.util.spec_from_file_location('historical_kine_ref',kpath);K=importlib.util.module_from_spec(spec);spec.loader.exec_module(K)
z=np.load(bind(O/'human_fresh.npz'),allow_pickle=False);n=int(np.sum(z['t']>0))+1;t=z['t'][:n];q=z['q'][:n];contact=z['contact'][:n];neuro=z['neuro'][:n]
bind(K.IK_MOT);bind(K.GRF_MOT)
ref=K.load_reference();mask=t>=2
result={'score_scope':'Composite kinematic score, including shape, range, mean, phase, period, contact and bilateral terms','sample_window_s':[float(t[mask][0]),float(t[mask][-1])],'whole_recording_s':[float(t[0]),float(t[-1])],'cycle_grid_samples':K.NPHASE,'LOAD_N':K.LOAD_N,'sides':{},'sources':sources}
saved=json.loads(bind(O/'overlay_metrics.json').read_text());savedscore=json.loads(bind(O/'human_metrics.json').read_text())
csv=np.genfromtxt(bind(O/'gait2392_ik_overlay.csv'),delimiter=',',names=True)
for side,col in [('r',0),('l',1)]:
    mean,periods,duty,cfrac,used,ncyc,on=K.sim_side(t,q,neuro,2.,side,contact)
    cols=(3,4,5)if side=='r'else(8,9,10);cycs,ps,raw=K._cycles_from_onsets(t[mask],t[mask],on,cols,q[mask])
    ti,names,vals=K._read_mot(K.IK_MOT);tg,gn,gv=K._read_mot(K.GRF_MOT);force=gv[:,gn.index('ground_force_vy'if side=='r'else'1_ground_force_vy')-1];beg,end=K._first_pair_in_ik(K._loading_onsets(tg,force),ti)
    values=dict(cycles=ncyc,model_period_mean_s=float(np.mean(periods)),model_contact_duty=duty,reference_period_s=ref['T_'+side],reference_duty=ref['duty_'+side])
    csv_error={}
    for j in ['hip','knee','ankle']:
        values['rmse_'+j+'_deg']=float(np.sqrt(np.mean((cycs[j]-ref[side][j])**2)))
        csv_error[j]={'mean':float(abs(cycs[j]-csv[f'{side}_{j}_model_mean_deg']).max()),'sd':float(abs(np.std(raw[j],axis=0)-csv[f'{side}_{j}_model_sd_deg']).max()),'reference':float(abs(ref[side][j]-csv[f'{side}_{j}_reference_deg']).max())}
    result['sides'][side]={'metrics':values,'saved_metric_differences':{k:float(values[k]-v)for k,v in saved[side].items()},'onsets_s':on.tolist(),'retained_periods_s':list(periods),'reference_loading_cycle_s':[beg,end],'csv_max_absolute_errors':csv_error,'contact_based':bool(used)}
compare=K.compare(t,q,neuro,2.,ref=ref,contact=contact)
result['score_recomputed']=compare['kine_score'];result['saved_score']=savedscore['kine_score'];result['score_difference']=compare['kine_score']-savedscore['kine_score'];result['score_terms']=compare
baseline=S/'reports_20261016/wean_human/finalist_eval.npz'
if baseline.exists():
    bz=np.load(bind(baseline),allow_pickle=False)
    result['original_finalist_numeric_comparison']={k:bool(z[k].dtype==bz[k].dtype and z[k].shape==bz[k].shape and z[k].tobytes()==bz[k].tobytes())for k in z.files if k in bz.files and np.issubdtype(z[k].dtype,np.number)}
for p in [ROOT/'Documentation/Reports and Papers/Dissertation/Figures/60-sns/Components/gait2392_ik_overlay/publication_ik_overlay_20261028.py',ROOT/'Documentation/Reports and Papers/Dissertation/Figures/95-AppendixD/Components/generate_gait2392_ik_overlay.py',ROOT/'Documentation/Reports and Papers/Dissertation/Figures/95-AppendixD/Components/gait2392_ik_overlay_provenance.json',S/'reports_20261016/wean_human/round8_best3.json',S/'reports_20261016/wean_human/finalist_spec.json']:
    if p.exists():bind(p)
result['elapsed_s']=time.perf_counter()-start
(OUT/'HISTORICAL_IK_RESULTS.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items()if k not in ['sources','score_terms']},indent=2))
