"""Decompose saved MN conductance currents; no new physical simulation."""
import os,sys,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;SPINAL=HERE.parents[1]/'spinal'
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
sys.path[:0]=[str(SPINAL),str(SPINAL/'w2l_cpg')];os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,mujoco
import params as P,build_network as B
from gait_standing_family import tie_weights,canonical_neuron_names
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');rows=[]
# These are the explicit --reflex-stand changes applied before each recorded build.
P.G['ii_set_pose']=1.;P.G['posture_to_mn']=0.
model=mujoco.MjModel.from_xml_path(str(SPINAL.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'))
acts=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
for family in ('s3k','w2lvar','syn6'):
 os.environ['AARL_NET']=family;n=B.build(acts,dt=.002,interleg=True);tie_weights(n);c=n.compiled;names=canonical_neuron_names(n)
 for condition in ('baseline_v2','afferent0_v1'):
  p=ROOT/('standing/'+family+'_'+condition);z=np.load(p.with_suffix('.apparatus.npz'));v=z['neuro'][-250:];assert v.shape[1]==len(names),(family,condition,v.shape,len(names),len(c.V))
  currents={};positive={};negative={}
  for post,pre in zip(*np.nonzero(c.g_max_non)):
   if not names[post].startswith('MN_'):continue
   span=c.e_hi[post,pre]-c.e_lo[post,pre];g=c.g_max_non[post,pre]*np.clip((v[:,pre]-c.e_lo[post,pre])/span,0,1);current=g*(c.del_e[post,pre]-v[:,post]);key=names[pre].split('_')[0]
   currents[key]=currents.get(key,0.)+float(np.mean(current))/len(acts);positive[key]=positive.get(key,0.)+float(np.mean(np.maximum(current,0)))/len(acts);negative[key]=negative.get(key,0.)+float(np.mean(np.minimum(current,0)))/len(acts)
  result=dict(family=family,condition=condition,last_samples=len(v),mean_per_MN_signed_current_nA=currents,mean_positive_current_nA=positive,mean_negative_current_nA=negative,scope='Analytic synaptic-current decomposition at saved late-standing neuron states. External POST_ currents and membrane leak are separate; no current cancellation applied. Rebuilt graph uses explicit reflex-stand flags and same mirrored pooling.')
  rows.append(result)
(ROOT/'standing/MN_CURRENT_DECOMPOSITION.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
