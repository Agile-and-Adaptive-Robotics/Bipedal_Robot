"""Read-only native Di Russo baseline check; no model is run."""
from pathlib import Path
import hashlib,json,numpy as np
R=Path(__file__).resolve()
while R.name!='Bipedal_Robot':R=R.parent
P=R/'Code/MuJoCo_SNS/spinal/reports_20261027_codex/native_baseline'
E=Path(r'D:\GitHub\_zcode_build_tools\dirusso_full')
sources={}
def h(p):
    b=p.read_bytes();sources[str(p)]={'sha256':hashlib.sha256(b).hexdigest(),'md5':hashlib.md5(b).hexdigest()};return b
a=P/'eval_spinal_native.sto.sto';b=E/'eval_repeat.sto.sto'
aa,bb=h(a),h(b)
same_after_header=aa.splitlines()[1:]==bb.splitlines()[1:]
with a.open()as f:
    lines=[f.readline()for _ in range(7)]
names=lines[6].strip().split('\t');data=np.loadtxt(a,skiprows=7);rep=np.loadtxt(b,skiprows=7)
col={n:i for i,n in enumerate(names)}
t=data[:,col['time']];tx=data[:,col['pelvis_tx']]
assert np.array_equal(data,rep)
summary={'rows':len(t),'columns':len(names),'data_exact_repeat':bool(np.array_equal(data,rep)),'text_exact_after_filename_header':same_after_header,'window_s':[float(t[0]),float(t[-1])],'dt_storage_s':float(t[1]-t[0]),'progression_m':float(tx[-1]-tx[0]),'speed_linear_fit_t_ge1_mps':float(np.polyfit(t[t>=1],tx[t>=1],1)[0]),'pelvis_height_range_m':[float(np.min(data[:,col['pelvis_ty']])),float(np.max(data[:,col['pelvis_ty']]))]}
for s in ['l','r']:
    y=data[:,col['leg0_l.grf_norm_y'if s=='l'else'leg1_r.grf_norm_y']]
    summary['foot_'+s]={'stance_fraction_gt0p1':float(np.mean(y>.1)),'peak_load_bw':float(np.max(y))}
for rel in ['NATIVE_BASELINE_20261027.md','native_build_manifest.txt','native_baseline_summary.json','summarize_native.py','eval_run.log','eval_repeat.log']:h(P/rel)
for rel in ['spinal_controller_main.par','spinal_controller_main.scone','import/initializations/spinal.par','import/initializations/f0918_init_walk_20170205.sto','import/models/gait9dof18musc_Thelen_20170320.osim']:h(E/rel)
expected={'import/initializations/spinal.par':'733a0355bd0376125a056c4cf75986a9','import/initializations/f0918_init_walk_20170205.sto':'d1b4898f74c5477aca3c74c87c2f3662','import/models/gait9dof18musc_Thelen_20170320.osim':'d958a0b9e7a32c6a2e3e8eeeefb88d5d'}
summary['input_md5_manifest_matches']={k:sources[str(E/k)]['md5']==v for k,v in expected.items()}
assert all(summary['input_md5_manifest_matches'].values())
out={'summary':summary,'sources':sources,'no_simulations':True}
(Path(__file__).parent/'HISTORICAL_NATIVE_RESULTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
