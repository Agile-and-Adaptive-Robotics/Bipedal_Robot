"""Canonical-right and jointly fitted bilateral static surrogate diagnostics.

Uses fresh OpenSim equilibrium data on a functional pose grid, rather than
assuming old Step3 transient samples are static truth. No production mutation.
"""
import os,sys,json,pickle,math,time
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parents[1];ROOT=H.parents[2];sys.path[:0]=[str(H),str(H/'functional_gait2392')]
import numpy as np,mujoco,opensim as osim,bsolve_ik as B
from scipy.optimize import least_squares
from mechanics import CoupledArms,JOINTS
OUT=Path(__file__).resolve().parent/'runs/audit_20261030';OUT.mkdir(parents=True,exist_ok=True)
source=ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'
om=osim.Model(str(source));om.setUseVisualizer(False);s=om.initSystem();cs=om.getCoordinateSet();oms=om.getMuscles()
m,d=B.build_air_model(False);ar=CoupledArms(m)
signs={r['coordinate']:r['sign'] for r in json.loads((H/'hunt_training/runs/ik_audit/fk_signs.json').read_text())}
families=('bifemsh','vas_med','vas_int','vas_lat','soleus','tib_ant');activations=(.05,.3,1.)
# Train on separated functional levels; interpolation at intermediate poses is held out.
poses=[dict(hip=h,knee=k,ankle=a,heldout=False) for h,k,a in ((0,0,0),(0,-30,0),(30,-60,0),(30,-90,0),(0,0,-20),(0,0,20),(30,-60,-20),(30,-60,20))]
poses += [dict(hip=15,knee=-45,ankle=10,heldout=True),dict(hip=15,knee=-75,ankle=-10,heldout=True)]
profile='base_v1';output_name='SHARED_FORCE_SURROGATE.json'
if '--hip-v2' in sys.argv:
 profile='hip_extension_v2';output_name='SHARED_FORCE_SURROGATE_HIP_V2.json'
 families=('psoas','iliacus','glut_max1','glut_max2','glut_max3','med_gas','lat_gas')
 poses=[dict(hip=h,knee=k,ankle=a,heldout=False) for h,k,a in
  ((-15,0,-10),(0,0,0),(15,-30,10),(30,-60,0),(45,-60,-10),(60,-90,20),(0,-60,-20),(30,-30,20),(0,-90,20),(45,0,-20))]
 poses += [dict(hip=h,knee=k,ankle=a,heldout=True) for h,k,a in
  ((-7.5,-15,5),(22.5,-45,10),(52.5,-75,-10),(15,-75,-5))]
if '--biarticular-v3' in sys.argv:
 profile='biarticular_extension_v3';output_name='SHARED_FORCE_SURROGATE_BIARTICULAR_V3.json'
 families=('semimem','semiten','bifemlh','rect_fem','tfl');activations=(.01,.05,.15,.3,1.)
 poses=[dict(hip=h,knee=k,ankle=0.,adduction=ad,heldout=False) for h,k,ad in
  ((-15,0,-10),(0,0,0),(0,-30,10),(15,-60,-10),(30,-90,0),(45,-30,10),(60,-60,-10),(30,0,10),(-15,-90,10),(60,0,0),(15,-90,-10),(45,-60,0))]
 poses += [dict(hip=h,knee=k,ankle=0.,adduction=ad,heldout=True) for h,k,ad in
  ((-7.5,-15,5),(22.5,-45,-5),(52.5,-75,5),(7.5,-75,5),(37.5,-15,-5))]
data={f:[] for f in families};failures=[]
for pose in poses:
 d.qpos[ar.addresses]=0.
 for name in JOINTS:
  base=name.rsplit('_',1)[0];angle=pose['hip'] if base=='hip_flexion' else pose['knee'] if base=='knee_angle' else pose['ankle'] if base=='ankle_angle' else pose.get('adduction',0.) if base=='hip_adduction' else 0.
  cs.get(name).setValue(s,math.radians(angle),False);d.qpos[m.joint(name).qposadr[0]]=math.radians(angle)*signs[name]
 ar.project(d);mujoco.mj_forward(m,d);om.realizePosition(s)
 lengths={f+'_'+side:float(d.actuator_length[m.actuator(f+'_'+side).id]) for f in families for side in ('r','l')}
 for activation in activations:
  for i in range(oms.getSize()):oms.get(i).setActivation(s,activation)
  try:om.equilibrateMuscles(s);om.realizeDynamics(s)
  except Exception as ex:failures.append(dict(pose=pose,activation=activation,error=str(ex)));continue
  for f in families:
   for side in ('r','l'):
    name=f+'_'+side;data[f].append(dict(side=side,pose=pose,activation=activation,length_m=lengths[name],opensim_tension_N=float(oms.get(name).getTendonForce(s))))

def tension(L,a,gp,bp,lr):
 return -a*mujoco.mju_muscleGain(L,0.,lr,1.,gp[:9])-mujoco.mju_muscleBias(L,lr,1.,bp[:9])

rows={};start=time.perf_counter()
for family,points in data.items():
 ids={side:m.actuator(family+'_'+side).id for side in ('r','l')}
 original={side:dict(gain=m.actuator_gainprm[i].copy(),bias=m.actuator_biasprm[i].copy(),dynamics=m.actuator_dynprm[i].copy(),lengthrange=m.actuator_lengthrange[i].copy()) for side,i in ids.items()}
 template=original['r']['gain'].copy();lr=original['r']['lengthrange'].copy();f0=float(oms.get(family+'_r').getMaxIsometricForce())
 train=np.array([not x['pose']['heldout'] for x in points]);held=~train
 target=np.array([x['opensim_tension_N'] for x in points]);baseline=[];canonical=[]
 for p in points:
  o=original[p['side']];baseline.append(tension(p['length_m'],p['activation'],o['gain'],o['bias'],o['lengthrange']))
  o=original['r'];canonical.append(tension(p['length_m'],p['activation'],o['gain'],o['bias'],o['lengthrange']))
 def predict(x):
  gp=template.copy();gp[[0,1,7]]=[x[0],x[1],x[2]];gp[2]=f0*x[3]
  return np.array([tension(p['length_m'],p['activation'],gp,gp,lr) for p in points])
 def residual(x):return (predict(x)[train]-target[train])/f0
 starts=[np.array([template[0],template[1],template[7],template[2]/f0]),np.array([.1,1.6,1.,1.]),np.array([.5,1.5,.5,1.])]
 fits=[least_squares(residual,np.clip(x,[.010001,1.000001,.010001,.500001],[.999999,1.989999,4.999999,1.999999]),bounds=([.01,1.,.01,.5],[1.,1.99,5.,2.]),max_nfev=1000) for x in starts]
 best=min(fits,key=lambda r:np.mean(r.fun**2));joint=predict(best.x)
 def stats(pred):
  pred=np.asarray(pred)
  return dict(train_RMSE_N=float(np.sqrt(np.mean((pred[train]-target[train])**2))),heldout_RMSE_N=float(np.sqrt(np.mean((pred[held]-target[held])**2))),
   train_normalized_RMSE=float(np.sqrt(np.mean((pred[train]-target[train])**2)))/f0,heldout_normalized_RMSE=float(np.sqrt(np.mean((pred[held]-target[held])**2)))/f0,
   maximum_bilateral_force_difference_N=float(np.max(abs(pred[0::2]-pred[1::2]))))
 gp=template.copy();gp[[0,1,7]]=best.x[:3];gp[2]=f0*best.x[3]
 rows[family]=dict(source_Fmax_N=f0,source_Fmax_left_N=float(oms.get(family+'_l').getMaxIsometricForce()),
  original_parameters={side:{k:v.tolist() for k,v in o.items()} for side,o in original.items()},shared_lengthrange=lr.tolist(),
  canonical_right=stats(canonical),baseline=stats(baseline),joint_fit=stats(joint),joint_fit_four_parameters=best.x.tolist(),joint_gainprm=gp.tolist(),joint_biasprm=gp.tolist(),joint_dynprm=original['r']['dynamics'].tolist(),
  points=[dict(p,baseline_tension_N=float(b),canonical_right_tension_N=float(c),joint_tension_N=float(j)) for p,b,c,j in zip(points,baseline,canonical,joint)])
 print(family,'baseline',rows[family]['baseline'],'joint',rows[family]['joint_fit'],flush=True)
report=dict(profile=profile,scope='Static equilibrium only; source Thelen tendon compliance causes activation-dependent fiber states. Shared rigid surrogate may not fit all activations exactly. Functional grid and interpolation holdouts; not full-ROM or dynamic validation. Canonical-right variant copies complete static law and lengthrange, not only Fmax; dynamics unchanged/unassessed. No production edits.',source_opensim=str(source),activations=activations,poses=poses,failures=failures,fit_wall_s=time.perf_counter()-start,muscles=rows)
(OUT/output_name).write_text(json.dumps(report,indent=2))
