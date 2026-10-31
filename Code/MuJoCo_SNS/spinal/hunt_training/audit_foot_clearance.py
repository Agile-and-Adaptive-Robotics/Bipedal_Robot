"""Replay foot-mesh clearance, independently of contact-force unloading."""
import os,sys,json,argparse
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H.parent)]
import numpy as np,mujoco,bsolve_ik as B
from scipy.spatial import ConvexHull
p=argparse.ArgumentParser();p.add_argument('recordings',nargs='+',type=Path);a=p.parse_args()
m,d=B.build_air_model(False);plane=m.geom('ground-plane').id;geoms=[]
for s in ('r','l'):
 gs=[]
 for g in range(m.ngeom):
  if m.body(int(m.geom_bodyid[g])).name not in ('calcn_'+s,'toes_'+s) or m.geom_type[g]!=mujoco.mjtGeom.mjGEOM_MESH:continue
  mesh=m.geom_dataid[g];start=m.mesh_vertadr[mesh];n=m.mesh_vertnum[mesh];v=m.mesh_vert[start:start+n].copy()
  if len(v)>4:v=v[ConvexHull(v).vertices]
  gs.append((g,v))
 geoms.append(gs)
for path in a.recordings:
 path=path.resolve();z=np.load(path,allow_pickle=True);result=json.loads(path.with_suffix('.result.json').read_text());assist=np.load(path.with_suffix('.assistance.npz'))
 idx=np.flatnonzero((z['t']>=5)&(z['t']<13))[::10];height=[]
 for i in idx:
  d.qpos[:]=z['qfull'][i];mujoco.mj_kinematics(m,d);normal=d.geom_xmat[plane].reshape(3,3)[:,2].copy()
  if normal[2]<0:normal=-normal
  height.append([min(float(np.min(v@(d.geom_xmat[g].reshape(3,3).T@normal))+(d.geom_xpos[g]-d.geom_xpos[plane])@normal) for g,v in gs) for gs in geoms])
 height=np.asarray(height);bw=result['measured']['body_weight_N'];force=assist['samples'][idx,7:9]/bw;rows=[]
 for side in (0,1):
  unloaded=force[:,side]<.025;peaks=[]
  bounds=np.r_[0,np.flatnonzero(unloaded[1:]!=unloaded[:-1])+1,len(unloaded)]
  for start,end in zip(bounds[:-1],bounds[1:]):
   if unloaded[start] and end-start>=4:peaks.append(float(height[start:end,side].max()))
  rows.append(dict(unloaded_samples=int(unloaded.sum()),median_clearance_when_unloaded_m=float(np.median(height[unloaded,side])) if unloaded.any() else None,swing_peak_clearances_m=peaks,swings_clearing_10mm=sum(v>=.01 for v in peaks),maximum_clearance_m=float(height[:,side].max())))
 report=dict(complete=result['assessment'].get('complete'),feet=rows,scope='Minimum visual foot-mesh signed plane distance, every20ms; force<.025BW for>=80ms. 10mm is an engineering diagnostic, not a published biological threshold. Does not replace movement/contact/frontal gates.')
 path.with_suffix('.clearance.json').write_text(json.dumps(report,indent=2));np.savez_compressed(path.with_suffix('.clearance.npz'),t=z['t'][idx],height_m=height,load_bw=force);print(path.stem,json.dumps(report),flush=True)
