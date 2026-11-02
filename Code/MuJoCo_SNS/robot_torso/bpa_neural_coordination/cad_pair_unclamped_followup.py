"""Short frozen reference/pressure repeats, live SNS state audit, actual mesh."""
import json,sys,argparse,subprocess,hashlib
from pathlib import Path
import numpy as np
import cad_pair_unclamped as U
BASE=U.OUT;C=U.C;L=U.L
def read(p):return json.loads(p.read_text())
def child(mode):
 gate=U.checks();assert gate['passed'];U.checks=lambda:gate;U.OUT=BASE/'fresh';U.OUT.mkdir(exist_ok=True);C.OUT=U.OUT;L.OUT=U.OUT;original=C.R.bn.build;holder=[];traces=[];verified=[]
 def build(*args,**kwargs):
  n=original(*args,**kwargs);holder.append(n);names=sorted(n.idx,key=n.idx.get);indices=[names.index(p+'_'+s) for s in ('r','l') for p in ('PF_E1','PF_E2','IBEXC_knee_ext')];step=n.step
  def observe(x):
   if L.ACTIVE is not None and L.ACTIVE.d.time>.004 and not verified:
    start=n.compiled.export_state();one=step(x).copy();end=n.compiled.export_state();n.compiled.import_state(start);two=step(x).copy();again=n.compiled.export_state();assert np.array_equal(one,two) and all(np.array_equal(end[k],again[k]) for k in end),'LiveSNScontinuation mismatch';verified.append(True);out=two
   else:out=step(x)
   traces.append((float(L.ACTIVE.d.time),n.compiled.V[indices].copy()));return out
  n.step=observe;return n
 C.R.bn.build=build;oldapply=L.Mixed.apply
 def apply(self,data):
  oldapply(self,data)
  if self.logs:
   if not hasattr(self,'mesh_vertices'):
    self.mesh_vertices=[]
    for g in range(self.m.ngeom):
     mid=self.m.geom_dataid[g]
     if (self.m.geom(g).name or '').startswith('foot_') and self.m.geom_type[g]==C.mujoco.mjtGeom.mjGEOM_MESH and mid>=0:self.mesh_vertices.append((g,self.m.mesh_vert[self.m.mesh_vertadr[mid]:self.m.mesh_vertadr[mid]+self.m.mesh_vertnum[mid]].copy()))
   assert self.mesh_vertices;mesh=min(float((v@data.geom_xmat[g].reshape(3,3)[1]+data.geom_xpos[g,1]).min()) for g,v in self.mesh_vertices);self.logs[-1].update(foot_mesh_min_y_m=np.array(mesh),qpos_full=data.qpos.copy(),qvel_full=data.qvel.copy())
 L.Mixed.apply=apply;r=U.run(mode);assert verified;folder=U.OUT/mode;np.savez_compressed(folder/'PRESYNAPTIC_TRACE.npz',time=[r[0] for r in traces],voltage=[r[1] for r in traces],names=np.array([p+'_'+s for s in ('r','l') for p in ('PF_E1','PF_E2','IBEXC_knee_ext')]));a=L.ACTIVE;states={}
 for prefix,obj,keys in [('bank_',a.bank,a.bank.__dict__.keys()),('mechanics_',a.d,('qpos','qvel','qacc_warmstart','act','ctrl','qfrc_applied','xfrc_applied','time','userdata','mocap_pos','mocap_quat'))]:
  for key in keys:
   value=getattr(obj,key)
   if isinstance(value,np.ndarray) and value.dtype.kind!='O':states[prefix+key]=value.copy()
   elif isinstance(value,(bool,int,float)):states[prefix+key]=np.asarray(value)
 np.savez_compressed(folder/'CAPTURED_ACTUATOR_MECHANICS_STATE.npz',**states);(folder/'LIVE_SNS_STATE_AUDIT.json').write_text(json.dumps(dict(exact=True,scope='One liveSNSforward snapshot→step→restore→repeat; all core/sparse/FreePair arrays includingneuralhistory exact. DoesnotrestorewholeCADrunnerlocals.'),indent=2));print(json.dumps(r,indent=2))
def compare(a,b):
 out={}
 for p in sorted(a.glob('*.npz')):
  with np.load(p) as x,np.load(b/p.name) as y:
   assert set(x.files)<=set(y.files);out[p.name]={k:bool(np.array_equal(x[k],y[k])) for k in x.files}
 assert out and all(v for row in out.values() for v in row.values()),'Freshpressure/reference trace mismatch';return out
def main():
 for mode in ('reference','mixed'):
  r=read(BASE/mode/'UNCLAMPED_SUMMARY.json');assert not r['failure'] and r['duration_s']>=2 and r['actual_POST_current_max_nA']==0 and r['equation_error_max_mV']<=1e-12 and max(r['knee_max_rad'])<=0
 sources=[Path(__file__),Path(U.__file__),BASE/'RECRUITMENT_PREREGISTRATION.json'];manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};path=BASE/'FRESH_SOURCE_MANIFEST.json'
 if path.exists():assert read(path)==manifest
 else:path.write_text(json.dumps(manifest,indent=2))
 for mode in ('reference','mixed'):
  done=BASE/'fresh'/mode/'LIVE_SNS_STATE_AUDIT.json'
  if not done.exists():
   with (BASE/('fresh_'+mode+'.log')).open('w') as log:p=subprocess.run([sys.executable,str(Path(__file__)),'--child',mode],stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
   if p.returncode:raise RuntimeError(str(BASE/('fresh_'+mode+'.log')))
  assert read(done)['exact']
 result=dict(matched_reference_fresh=compare(BASE/'reference',BASE/'fresh/reference'),mixed_pressure_fresh=compare(BASE/'mixed',BASE/'fresh/mixed'),live_neural_continuation_both_exact=True,whole_runner_restart_tested=False)
 with np.load(BASE/'fresh/mixed/ACTUATOR_BRIDGE.npz') as z:result['actual_pressure_model_foot_mesh_y_range_m']=[float(z['foot_mesh_min_y_m'].min()),float(z['foot_mesh_min_y_m'].max())]
 (BASE/'FRESH_REPEAT_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(reference_exact=True,pressure_exact=True,live_SNS_state_exact=True,mesh=result['actual_pressure_model_foot_mesh_y_range_m']),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--child',choices=['reference','mixed']);a=p.parse_args()
 if a.child:child(a.child)
 else:main()
