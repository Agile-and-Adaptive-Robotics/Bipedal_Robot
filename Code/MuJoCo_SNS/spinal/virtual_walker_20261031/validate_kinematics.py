"""Check moving source poses, not just a visually plausible initial rendering."""
import argparse, json, hashlib
from pathlib import Path
import numpy as np
from import_virtual_walker import HERE,REPO,build,safe,mujoco,FRAME
def run(newer=False):
 path=build(newer);folder=path.parent;meta=json.loads((folder/'IMPORT_AUDIT.json').read_text());source=Path(meta['source']).parent/'verification.npz'
 raw=source.read_bytes();sha=hashlib.sha256(raw).hexdigest();(folder/f'native_verification_{sha[:12]}.npz').write_bytes(raw)
 a=np.load(source);z=a['data'];labels=a['labels'].tolist();m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m)
 bodies=['foot_L','toe_L','foot_R','toe_R','Left hip wing','Right hip wing'];ids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,safe(s)) for s in bodies]
 results={}
 for sign in [1,-1]:
  errors=[];root_errors=[]
  for k in range(0,len(z),50):
   row=z[k];d.qpos[:]=m.qpos0;d.qpos[:7]=[0,0,0,1,0,0,0]
   for name in ['hip_L','knee_L','ankle_L','toe_L','hip_R','knee_R','ankle_R','toe_R']:
    ji=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,name);d.qpos[m.jnt_qposadr[ji]]=sign*row[labels.index(name+'_angle')]
   mujoco.mj_forward(m,d);X=d.xpos[ids].copy();Y=np.array([[row[labels.index(b+'_'+q)] for q in 'xyz'] for b in bodies])@FRAME.T;xc=X-X.mean(0);yc=Y-Y.mean(0);U,s,Vt=np.linalg.svd(xc.T@yc);R=Vt.T@U.T
   if np.linalg.det(R)<0:Vt[-1]*=-1;R=Vt.T@U.T
   shift=Y.mean(0)-R@X.mean(0);err=(X@R.T+shift)-Y;errors.append(np.sqrt(np.mean(err**2)));root_errors.append(np.linalg.norm(shift-FRAME@np.array([row[labels.index('root_'+q)] for q in 'xyz'])))
  results['source_coordinate_sign_'+str(sign)]=dict(marker_rms_max_m=float(max(errors)),marker_rms_mean_m=float(np.mean(errors)),root_position_error_max_m=float(max(root_errors)),samples=len(errors))
 best=results['source_coordinate_sign_1'];passed=best['marker_rms_max_m']<.005 and best['root_position_error_max_m']<.005
 out=dict(machine='easteregg2',source_native_trace=str(source),source_native_trace_sha256=sha,method='Apply source relative joint angles; fit one rigid global root transform to six physical markers after explicit source→anatomical Z reflection. Correct source joint sign must preserve limb shape, not just match the assembled pose.',axis_convention='Source→target S=diag(1,1,−1). Rotations S R S; hinge axial vectors transform det(S) S a; slide directions S a. Source coordinates are preserved.',passed=passed,results=results)
 (folder/'KINEMATIC_PARITY.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));return passed
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--newer',action='store_true');args=p.parse_args();raise SystemExit(0 if run(args.newer) else 1)
