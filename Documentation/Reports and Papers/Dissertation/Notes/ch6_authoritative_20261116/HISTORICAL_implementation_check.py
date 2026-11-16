"""Read-only checks of retained implementation recordings; no models are run."""
from pathlib import Path
import hashlib,json,numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks
R=next(p for p in Path(__file__).resolve().parents if p.name=='Bipedal_Robot')
S=R/'Code/Matlab/SNS_Simscape';W=R/'Code/MuJoCo_SNS/spinal/campaigns/20260930'
sources={}
def hashed(p):
    p=Path(p);sources[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return p
bridge=loadmat(hashed(S/'mujoco_bridge/logs/bridge_tests_abc.mat'),simplify_cells=True)
gt=loadmat(hashed(S/'mujoco_bridge/logs/gt_const.mat'),simplify_cells=True)
fields=list(bridge['a']['vals']);mat=np.column_stack([bridge['a']['vals'][k]for k in fields])
assert np.array_equal(mat[1:],gt['sens'])
result={'bridge':{'fields':fields,'simulink_samples':len(bridge['a']['time']),'scored_steps':len(gt['sens']),'scored_times_s':[float(bridge['a']['time'][1]),float(bridge['a']['time'][-1])],'maximum_sensor_errors':np.max(abs(mat[1:]-gt['sens']),axis=0).tolist(),'exact':True,'pass_a':bool(bridge['pass_a']),'pass_b':bool(bridge['pass_b']),'pass_c':bool(bridge['pass_c']),'sensor_dt_s':float(np.mean(np.diff(bridge['b']['time']))),'two_rate_samples':len(bridge['b']['time'])}}
sl=loadmat(hashed(S/'results/SNS_W2L_CPG_run_20260930.mat'),simplify_cells=True)
with np.load(hashed(W/'w2l_numpy_ref.npz'),allow_pickle=False)as npz:
    vn=npz['V'];tn=npz['t'];names=npz['names'].tolist()
ts=sl['t'];wm=ts>=2;wn=tn>=2
tr={'L_RG_ext':'vLE','L_RG_flx':'vLF','R_RG_ext':'vRE','R_RG_flx':'vRF'}
comparison={}
for name,f in tr.items():
    j=names.index(name.replace('_',' '));a=sl[f][wm];b=vn[wn,j];n=min(len(a),len(b));a=a[:n];b=b[:n]
    aa=a-a.mean();bb=b-b.mean();lags=range(-100,101)
    cc=[float(np.corrcoef(aa[max(0,-l):n-max(0,l)],bb[max(0,l):n-max(0,-l)])[0,1])for l in lags]
    best=int(np.nanargmax(cc));on_a=a>.5*np.max(a);on_b=b>.5*np.max(b)
    ta=ts[wm][:n][np.flatnonzero(on_a[1:]&~on_a[:-1])+1]
    tb=tn[wn][:n][np.flatnonzero(on_b[1:]&~on_b[:-1])+1]
    comparison[name]={'best_lag_ms':(best-100)*2.,'best_lag_correlation':cc[best],'unshifted_correlation':float(np.corrcoef(a,b)[0,1]),'rmse_mV':float(np.sqrt(np.mean((a-b)**2))),'bursts_sim':len(ta),'bursts_numpy':len(tb),'period_sim_s':float(np.mean(np.diff(ta))),'period_numpy_s':float(np.mean(np.diff(tb)))}
result['w2l']={'scoring_window':'t>=2s through last native sample; truncate Simulink6001point grid to matching NumPy6000point grid, 5000 scoring points','source_grid_s':.002,'comparison':comparison,'saved_verdict':bool(sl['pass']),'interpretation':'Shared RG signals, representative Simulink core versus full NumPy network; best-lag correlation is not zero-lag agreement. Reference JSON protocol label is stale: code sends alternating heel trains.'}
deng=loadmat(hashed(S/'results/deng_cpg_ref.mat'),simplify_cells=True)
td=deng['t']/1000. # Original exported clock is milliseconds, as deng_cpg_ode.py specifies.
peaks,_=find_peaks(deng['traces'][:,0],prominence=1.)
peaks=peaks[td[peaks]>.3]
result['deng_numpy']={'samples':len(td),'window_s':[float(td[0]),float(td[-1])],'native_clock_units':'ms; seconds=t/1000','bursts_after_0p3s':len(peaks),'period_s':float(np.mean(np.diff(td[peaks]))),'extensor_flexor_correlation':float(np.corrcoef(deng['traces'][:,0],deng['traces'][:,1])[0,1]),'tau_h_mode':deng['tau_h_mode'],'names':deng['names'].tolist(),'scope':'Raw NumPy reference only; Simulink pointwise arrays not retained with this run, so Simulink1.938s/-0.857 is checked via source script and primary log. sns_run_deng_demo.m compares seconds directly to rawreference milliseconds, so its printed maxdeltaV is not a valid matched-clock measure.'}
for p in [W/'w2l_trace_compare.json',W/'w2l_trace_compare_20260930.py',W/'w2l_numpy_ref_20260930.py',W/'w2l_numpy_ref.json',S/'dev/build_sns_w2l_cpg_20260930.m',S/'results/SNS_W2L_CPG.slx',S/'mujoco_bridge/matlab/run_bridge_tests.m',S/'mujoco_bridge/matlab/bridge_groundtruth.py',S/'mujoco_bridge/matlab/e2_cpg_open_loop.m',S/'mujoco_bridge/logs/test_abc.log',S/'mujoco_bridge/logs/test_abc_laptop.log',S/'mujoco_bridge/logs/e2_cpg.log',S/'mujoco_bridge/matlab/exp_out/e2_cpg.mat',S/'demos/sns_run_deng_demo.m',S/'demos/SNS_Deng_CPGDemo.slx',S/'logs/sns_run_deng2.log',S/'logs/sns_verify.log',S/'logs/desktop_capability_20261024_b.log',S/'sns_verify_from_json.m',S/'results/verify_ref.mat',S/'results/SNS_SpinalNetwork.slx',S/'results/SNS_SpinalNetwork.slx.r2025a',R/'Code/MuJoCo_SNS/spinal/deng_cpg_ode.py']:
    if p.exists():hashed(p)
result['source_sha256']=sources
out=Path(__file__).parent/'HISTORICAL_IMPLEMENTATION_RESULTS.json';out.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items()if k!='source_sha256'},indent=2))
