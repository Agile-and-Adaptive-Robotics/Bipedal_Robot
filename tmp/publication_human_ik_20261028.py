from pathlib import Path
import os,sys,json
R=Path('D:/GitHub/Bipedal_Robot');S=R/'Code/MuJoCo_SNS/spinal';O=S/'reports_20261028_publication_flow/human_ik';O.mkdir(exist_ok=True)
for n in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[n]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent);os.environ.pop('AARL_NET',None);os.chdir(S);sys.path.insert(0,str(S))
import matplotlib
matplotlib.use('Agg')
import matplotlib.figure
matplotlib.figure.Figure.savefig=lambda *args,**kw:None
import wean_tune_20261016 as W
winner=json.loads((S/'reports_20261016/wean_human/round8_best3.json').read_text());fb={k:float(winner['params'].get(k,0)) for k in W.FB_KEYS}
v=W.run_eval('publication_fresh',winner['knobs'],winner['scale'],fb,npz=str(O/'human_fresh.npz'))
(O/'human_metrics.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
