"""Fresh fixed5% CAD repeats with purely observational mesh/state capture."""
import argparse,json,sys,subprocess,hashlib,shutil
from pathlib import Path
import numpy as np
import cad_mixed_coactivation05 as J
L=J.L;C=J.C;BASE=J.OUT
def read(p):return json.loads(p.read_text())
def child(name):
 parent=BASE/name;parent.mkdir(exist_ok=True)
 for file in ('ASSUMED_PAIR_CONFIGURATION.json','STATIC_AUDIT.json'):shutil.copyfile(BASE/file,parent/file)
 (parent/'source_snapshot').mkdir(exist_ok=True);C.OUT=parent;L.OUT=parent;holder=[];original_build=C.R.bn.build
 def observe_build(*args,**kwargs):
  n=original_build(*args,**kwargs);holder.append(n);return n
 C.R.bn.build=observe_build;original_apply=L.Mixed.apply
 def observe(self,data):
  original_apply(self,data)
  if self.logs:
   if not hasattr(self,'foot_vertices'):
    self.foot_vertices=[]
    for g in range(self.m.ngeom):
     mid=self.m.geom_dataid[g]
     if (self.m.geom(g).name or '').startswith('foot_') and self.m.geom_type[g]==C.mujoco.mjtGeom.mjGEOM_MESH and mid>=0:
      vert=self.m.mesh_vert[self.m.mesh_vertadr[mid]:self.m.mesh_vertadr[mid]+self.m.mesh_vertnum[mid]].copy();self.foot_vertices.append((g,vert))
   assert self.foot_vertices,'Actual visualfoot meshes required'
   clearance=min(float(np.min(v@data.geom_xmat[g].reshape(3,3)[1]+data.geom_xpos[g,1])) for g,v in self.foot_vertices);self.logs[-1].update(actual_foot_mesh_min_y_m=np.array(clearance),qpos_full=data.qpos.copy(),qvel_full=data.qvel.copy())
 L.Mixed.apply=observe
 summary=J.live('mixed');a=L.ACTIVE;states={};net=holder[0]
 for prefix,obj in [('sparse_',net.compiled),('backend_',net.compiled.c),('bank_',a.bank),('mechanics_',a.d)]:
  keys=obj.__dict__.keys() if hasattr(obj,'__dict__') else ('qpos','qvel','qacc_warmstart','act','ctrl','qfrc_applied','xfrc_applied','time','userdata','mocap_pos','mocap_quat')
  for k in keys:
   value=getattr(obj,k)
   if isinstance(value,np.ndarray) and value.dtype.kind!='O':states[prefix+k]=value.copy()
   elif isinstance(value,(bool,int,float)):states[prefix+k]=np.asarray(value)
 np.savez_compressed(parent/'mixed/FINAL_COMPLETE_STATE.npz',**states);print(json.dumps(summary,indent=2))
def compare(a,b,subset=False):
 checks={}
 for p in sorted(a.glob('*.npz')):
  other=b/p.name;assert other.exists(),str(other)
  with np.load(p) as x,np.load(other) as y:
   assert set(x.files)<=set(y.files) if subset else set(x.files)==set(y.files)
   checks[p.name]={k:bool(np.array_equal(x[k],y[k])) for k in x.files}
 assert checks and all(ok for r in checks.values() for ok in r.values()),'CAD fresh mismatch'
 return checks
def main():
 initial=read(BASE/'mixed/BRIDGE_SUMMARY.json');assert not initial['failure'] and initial['duration_s']>=2 and max(initial['knee_max_rad'])<=0 and all(r[1]>100 for r in initial['pair_force_range_N'])
 manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(J.__file__),Path(L.__file__),Path(C.__file__),BASE/'ASSUMED_PAIR_CONFIGURATION.json',BASE/'COACTIVATION_REGISTRATION.json')};path=BASE/'FRESH_SOURCE_MANIFEST.json'
 if path.exists():assert read(path)==manifest
 else:path.write_text(json.dumps(manifest,indent=2))
 for name in ('fresh_a','fresh_b'):
  case=BASE/name/'mixed/BRIDGE_SUMMARY.json'
  if not case.exists():
   with (BASE/(name+'.log')).open('w') as log:p=subprocess.run([sys.executable,str(Path(__file__)),'--child',name],stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
   if p.returncode:raise RuntimeError(str(BASE/(name+'.log')))
  r=read(case);assert not r['failure'] and r['duration_s']>=2 and max(r['knee_max_rad'])<=0
 old=compare(BASE/'mixed',BASE/'fresh_a/mixed',subset=True);fresh=compare(BASE/'fresh_a/mixed',BASE/'fresh_b/mixed');mesh=[]
 with np.load(BASE/'fresh_a/mixed/ACTUATOR_BRIDGE.npz') as z:mesh=[float(z['actual_foot_mesh_min_y_m'].min()),float(z['actual_foot_mesh_min_y_m'].max())]
 result=dict(original_physical_arrays_exact=old,independent_fresh_all_arrays_and_states_exact=fresh,actual_mesh_min_y_range_m=mesh,actual_mesh_sample_dt_s=.002,scope='Samefixed5% engineeringhybrid; observationonly additions. Strongapparatus/reference52/assumedBPA remain explicit. No independentbalanceorautonomousneuralclaim.');(BASE/'FRESH_AND_MESH_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(fresh_exact=True,mesh_range_m=mesh),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--child');a=p.parse_args()
 if a.child:child(a.child)
 else:main()
