"""Read-only audit of saved timing pairs and MATLAB numeric recordings.

Run with D:/Anaconda/python.exe. No benchmark or scientific simulation is run.
The bounded MCOS extraction reads the saved FileWrapper property table in its
stored order, verifies timeseries object IDs and sample counts, and checks the
spiking file's independently saved final/max scalars. No time interpolation is
used for the cup statistics.
"""
from pathlib import Path
import hashlib, io, json, time, sys, subprocess
import numpy as np
from scipy.io import loadmat
from scipy.io.matlab._mio5 import MatFile5Reader

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if p.name=='Bipedal_Robot')
ACC=Path(r'D:\Users\Ben Bolen\NeuroAcceleration_20261104')
SRC=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261102_9h')
HASHES={}

def hashed(p):
    p=Path(p);HASHES[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return p

def js(p):return json.loads(hashed(p).read_text())

def compare_npz(a,b):
    a,b=hashed(a),hashed(b)
    out={'a':str(a),'b':str(b),'numeric':{},'nonnumeric':[]}
    with np.load(a,allow_pickle=False)as x,np.load(b,allow_pickle=False)as y:
        assert set(x.files)==set(y.files)
        for k in x.files:
            try:xa,ya=x[k],y[k]
            except ValueError:
                out['nonnumeric'].append(k);continue
            if xa.dtype.kind not in 'biufc':out['nonnumeric'].append(k);continue
            same=np.array_equal(xa,ya,equal_nan=True)
            out['numeric'][k]={'shape':list(xa.shape),'exact':bool(same)}
            assert same,(a,k)
        for k in ['t','t_state','t_neural','physics_clock_s']:
            if k in x.files:out[k+'_range']=[float(x[k][0]),float(x[k][-1])]
    return out

start=time.perf_counter()
result={'execution':'Read-only saved-array audit; no simulations or timing reruns.','python':sys.executable,'numpy':np.__version__}
trace_arrays={}
source_pairs=[]
for suffix in ['.npz','.assistance.npz','.contact_events.npz','.regional_contacts.npz','.mechanics.npz']:
    source_pairs.append(compare_npz(SRC/('configurations/source_baseline_s1'+suffix),SRC/('configurations/source_native_s1'+suffix)))
source_pairs.append(compare_npz(SRC/'figures/source_baseline_s1/aligned_trace.npz',SRC/'figures/source_native_s1/aligned_trace.npz'))
source_runs={}
for name in ['source_native_s1','source_baseline_s1']:
    d=js(SRC/f'configurations/{name}.result.json')
    c=js(SRC/f'configurations/{name}.json')
    source_runs[name]={'wall_s':d.get('wall_s'),'simulated_s':d.get('simulated_s'),'compiled_scalars':{k:v for k,v in d.get('compiled',{}).items()if isinstance(v,(str,int,float,bool))},'configuration':c}
    source_runs[name]['external_process_timing']=js(SRC/f'timings/{name}/RESULT.json')
    meta=js(SRC/f'figures/{name}/ALIGNED_TRACE_METADATA.json')
    source_runs[name]['capture_metadata']=meta
    model=SRC/f'figures/{name}/model.mjb'
    if model.exists():hashed(model)
result['source_native_sparse']={'runs':source_runs,'comparisons':source_pairs,'numeric_array_count':sum(len(p['numeric'])for p in source_pairs),'repeat_count':1}
compact_pairs=[]
for suffix in ['.npz','.assistance.npz','.balance.npz','.mechanics.npz','.regional_contacts.npz','.contact_events.npz','/capture/aligned_trace.npz']:
    compact_pairs.append(compare_npz(ACC/('body/sparse_final'+suffix),ACC/('body/compact_final'+suffix)))
compact_runs={}
for name in ['sparse_final','compact_final']:
    d=js(ACC/f'body/{name}.result.json')
    compact_runs[name]={'wall_s':d['wall_s'],'simulated_s':d['simulated_s'],'accepted':d['assessment']['accepted'],'compiled_scalars':{k:v for k,v in d['compiled'].items()if isinstance(v,(str,int,float,bool))},'configuration':js(ACC/f'body/{name}.json')}
    compact_runs[name]['model_sha256']=hashlib.sha256(hashed(ACC/f'body/{name}/capture/model.mjb').read_bytes()).hexdigest()
result['compact_body']={'runs':compact_runs,'comparisons':compact_pairs,'numeric_array_count':sum(len(p['numeric'])for p in compact_pairs),'repeat_count':1}
assert compact_runs['sparse_final']['model_sha256']==compact_runs['compact_final']['model_sha256']
physical_eval=hashed(ROOT/'Code/MuJoCo_SNS/spinal/hunt_training/physical_eval.py')
result['physical_eval_current_sha256']=HASHES[str(physical_eval)]
result['physical_eval_saved_sha256']={name:js(ACC/f'body/{name}.result.json')['source_sha256'].get('physical_eval.py')for name in ['sparse_final','compact_final']}
recorded_commit='9d7b8d4d50d6e7434eb2bd9743bf160578231de3'
recorded_bytes=subprocess.check_output(['git','show',recorded_commit+':Code/MuJoCo_SNS/spinal/hunt_training/physical_eval.py'],cwd=ROOT)
recorded_bytes=recorded_bytes.replace(b'\r\n',b'\n').replace(b'\n',b'\r\n')
recorded_hash=hashlib.sha256(recorded_bytes).hexdigest()
assert all(h==recorded_hash for h in result['physical_eval_saved_sha256'].values())
recorded_lines=recorded_bytes.decode('utf-8').splitlines()
result['compact_recorded_timer_source']={'git_commit':recorded_commit,'path':'Code/MuJoCo_SNS/spinal/hunt_training/physical_eval.py','line_endings':'CRLF (exact saved hash restored from LF git blob)','sha256':recorded_hash,'timer_lines':[{'line':i+1,'text':line}for i,line in enumerate(recorded_lines)if 'start=time.perf_counter()'in line or 'wall_s=time.perf_counter()-start'in line or 'write_text(json.dumps(result'in line], 'scope':'Evaluator timer starts before R.main, includes capture save and assessment, excludes interpreter imports/pre-run setup and final result/sidecar serialization.'}
result['source_compiled_models_equal']=HASHES[str(SRC/'figures/source_baseline_s1/model.mjb')]==HASHES[str(SRC/'figures/source_native_s1/model.mjb')]
assert result['source_compiled_models_equal']
result['kernels']={}
for family in ['hunt','s3k','w2lvar','syn6']:
    d=js(ACC/f'results/KERNEL_BENCHMARK_{family}.json')[family]
    for timed in d['timings_3000steps'].values():assert timed['median_s']==float(np.median(timed['repeats_s']))
    result['kernels'][family]=d
hashed(ACC/'source/benchmark.py')
hashed(ACC/'source/summarize_verification.py')
old=ROOT/'Code/MuJoCo_SNS/spinal/reports_20261027_codex'
result['historical_w2l']={'benchmark':js(old/'benchmark_exact.json'),'script':str(hashed(old/'benchmark_exact.py')),'dense_fast_pair':compare_npz(old/'trial458_d30.npz',old/'trial458_exact_fast_d30.npz'),'timing_repeats':1,'timed_steps_per_backend':10000,'electrical_comparison_steps':5000}

def cup_timeseries(p):
    d=loadmat(hashed(p),simplify_cells=True)
    raw=d['__function_workspace__'].tobytes()
    # The subsystem uses an abbreviated 8-byte MAT header. Supply only the
    # standard text/header padding; preserve every saved matrix byte.
    stream=io.BytesIO(b' '*124+raw[:4]+raw[8:])
    reader=MatFile5Reader(stream,simplify_cells=True)
    reader.initialize_read();reader.read_file_header()
    header,endpos=reader.read_var_header();wrapper=reader.read_var_array(header)
    values=wrapper.MCOS['arr'][0]
    fields=['th','mCup','F','A','Ft','At']
    if 'spk'in d['runs'][0]:fields+=['spk']
    runs=[]
    ordinal=0
    for r in d['runs']:
        rec={'label':r['label'],'k':int(r['k']),'series':{}}
        for f in fields:
            object_id=int(r[f]['arr'][0][4])
            assert object_id==ordinal*5+1,(f,object_id,ordinal)
            ti=10+18*ordinal;yi=ti+9
            t=np.asarray(values[ti]).reshape(-1);y=np.asarray(values[yi]).reshape(-1)
            assert t.shape==y.shape and np.all(np.diff(t)>=0)
            assert t[0]==0 and t[-1]==10
            rec['series'][f]=(t,y)
            ordinal+=1
        t,th=rec['series']['th'];m=rec['series']['mCup'][1]
        deg=np.rad2deg(th);w=t>2
        wi=np.where(w)[0][np.argmax(abs(deg[w]))]
        stats={'samples':len(t),'final_time_s':float(t[-1]),'final_load_kg':float(m[-1]),'endpoint_sag_deg':float(deg[-1]),'maximum_abs_sag_after_2s_deg':float(abs(deg[wi])),'maximum_abs_sag_after_2s_at_s':float(t[wi]),'maximum_abs_sag_after_2s_load_kg':float(np.interp(t[wi],*rec['series']['mCup'])),'timeseries_object_count':len(fields),'load_ramp_max_abs_error_kg':float(np.max(abs(m-.0625*rec['series']['mCup'][0])))}
        if 'final'in r:
            assert abs(float(r['final'])-stats['endpoint_sag_deg'])<1e-12
            assert abs(float(r['sagAfter2'])-stats['maximum_abs_sag_after_2s_deg'])<1e-12
            stats['saved_scalar_checks_exact_within_1e_12']=True
        rec['stats']=stats
        runs.append(rec)
    return runs

result['cup']={}
for name,p in [('graded',ROOT/'Code/Matlab/SNS_Simscape/results/sns_beer_results.mat'),('spiking',ROOT/'Code/MuJoCo_SNS/spinal/reports_spiking_20261002/goal2_beer_spiking.mat')]:
    runs=cup_timeseries(p)
    result['cup'][name]={'source':str(p),'runs':[{k:v for k,v in r.items()if k!='series'}for r in runs]}
    for r in runs:
        condition='on'if r['k']else'off'
        for field,(t,y)in r['series'].items():
            trace_arrays[f'cup_{name}_{condition}_{field}_t']=t
            trace_arrays[f'cup_{name}_{condition}_{field}_y']=y
for f in ['sns_build_beer_demo.m','sns_run_beer_demo.m','sns_build_beer_spiking.m','sns_run_beer_spiking.m']:
    hashed(ROOT/'Code/Matlab/SNS_Simscape/demos'/f)

knee_path=ROOT/'Code/Matlab/SNS_Simscape/results/knee_transient_reflex_20261028/replay_summary.mat'
knee=loadmat(hashed(knee_path),simplify_cells=True)
t=knee['t'];c=knee['C'];p=knee['P'];f=knee['F'];names=knee['names'].tolist()
trace_arrays.update(knee_t=t,knee_C=c,knee_P=p,knee_F=f,knee_names=np.asarray(names,dtype='U'))
diff=-np.rad2deg(p[:,0]-c[:,0]);window=t>=30
peak=np.where(window)[0][np.argmax(diff[window])]
# Saved column0 is internal angle; retain explicit paper sign conversion.
base=float(np.mean(c[(t>=25)&(t<30),0]))
baseline=float(-np.rad2deg(base))
after=np.where(t>=30)[0]
displacement=-np.rad2deg(p[:,0]-base)
good=(abs(displacement)<=.01)&(abs(p[:,1])<=1e-4)
bad=np.where((t>=30.03)&~good)[0]
recover=bad[-1]+1
result['knee']={'source':str(knee_path),'names':names,'samples':len(t),'window_s':[float(t[0]),float(t[-1])],'pulse_window_s':[30.,30.03],'baseline_paper_deg':baseline,'extension_peak_deg':float(diff[peak]),'peak_latency_s':float(t[peak]-30),'recovery_s':float(t[recover]-30),'saved_summary':knee['summary'],'late_residual_paper_deg':float(np.mean(displacement[t>=45])),'fine_max_angle_error_deg':float(np.max(abs(p[:,0]-f[:,0]))*180/np.pi),'criteria':'Displacement from mean control baseline25<=t<30 <=0.01deg; stimulated absolute angular velocity<=1e-4rad/s; every subsequent resampled point through50s.'}
hashed(knee_path.parent/'run_transient.m')
hashed(ROOT/'Code/Matlab/SNS_Simscape/demos/sns_build_demo.m')
result['elapsed_read_audit_s']=time.perf_counter()-start
result['source_sha256']=HASHES
trace_file=ROOT/'tmp/ch6_authoritative_20261116/TIMING_trace_arrays.npz'
trace_file.parent.mkdir(parents=True,exist_ok=True)
np.savez_compressed(trace_file,**trace_arrays)
result['derived_untransformed_trace_arrays']={'path':str(trace_file),'sha256':hashlib.sha256(trace_file.read_bytes()).hexdigest(),'scope':'Exact extracted saved cup samples and knee resampled recording arrays; no additional interpolation, smoothing, signs, rescaling, or decimation.'}
out=HERE/'TIMING_RESULTS.json';out.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({'source_arrays':result['source_native_sparse']['numeric_array_count'],'compact_arrays':result['compact_body']['numeric_array_count'],'source_runs':{k:{z:v[z]for z in ['wall_s','simulated_s','compiled_scalars']}for k,v in source_runs.items()},'compact_runs':{k:{z:v[z]for z in ['wall_s','simulated_s','compiled_scalars']}for k,v in compact_runs.items()},'cup':result['cup'],'knee':result['knee'],'elapsed_s':result['elapsed_read_audit_s']},indent=2))
