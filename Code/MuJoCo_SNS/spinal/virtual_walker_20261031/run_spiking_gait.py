"""Fresh external Gait2392 corrected-pose source vs pure-input spiking mirror."""
import argparse,contextlib,hashlib,json,runpy,sys,os
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['AARL_NET']='w2lvar';os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from import_virtual_walker import HERE
from spiking_sns_inputs import SpikeInputForward
H=HERE.parent/'hunt_training';sys.path[:0]=[str(H),str(HERE.parent),str(HERE.parent/'functional_gait2392')]
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['graded','hybrid'],required=True);p.add_argument('--output-root',default=r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\gait2392');p.add_argument('--population-size',type=int,default=16);p.add_argument('--tag',default='corrected');a=p.parse_args();out=Path(a.output_root);out.mkdir(parents=True,exist_ok=True)
stem=a.tag+'_'+a.mode
if (out/(stem+'.result.json')).exists():raise FileExistsError('Preserving completed result; choose a new --tag for a fresh repeat')
source=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/aligned_support_ladder_v1/repeat_0p5.json';cfg=json.loads(source.read_text());cfg['trace_cells']=['RG_E_r','RG_F_r','RG_E_l','RG_F_l','MN_soleus_r','MN_soleus_l'];path=out/(stem+'.json');path.write_text(json.dumps(cfg,indent=2));saved=[]
import controller as C
original=C.build
def build(*args,**kwargs):
 net=original(*args,**kwargs)
 if a.mode=='hybrid':
  names=['']*len(net.compiled.V)
  for name,i in net.idx.items():
   if i<len(names) and not names[i]:names[i]=name
  net.compiled=SpikeInputForward(net.compiled,names,a.population_size);saved.append(net.compiled)
 return net
C.build=build;sys.argv=[str(H/'physical_eval.py'),str(path)]
with (out/(stem+'.log')).open('w',encoding='utf-8') as log,contextlib.redirect_stdout(log):runpy.run_path(str(H/'physical_eval.py'),run_name='__main__')
result=json.loads(path.with_suffix('.result.json').read_text());record=dict(mode=a.mode,source_configuration=str(source),source_configuration_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),result=result,scope='Same corrected pose,92-entry human model, trained weights, force law and measured support. Spiking condition modifies pure-input neural encoding only; no spiking central-controller claim.')
if saved:
 c=saved[-1];record['spiking']=c.metadata();np.savez_compressed(out/(stem+'.spikes.npz'),spike_counts=np.array(c.logs),names=record['spiking']['spiking_source_names'])
(out/(stem+'.SPIKING.json')).write_text(json.dumps(record,indent=2));print(json.dumps(dict(mode=a.mode,accepted=result['assessment']['accepted'],simulated_s=result['simulated_s'],wall_s=result['wall_s'],spiking_cells=record.get('spiking',{}).get('actual_spiking_cells',0)),indent=2))
