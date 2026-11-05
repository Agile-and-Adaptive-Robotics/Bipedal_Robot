"""Compare native recorded geometry against imported FK, independent of dynamics.

Fit one proper root rotation per sample, since native root orientation was not
recorded. Compare both signs of native revolute displacement, without changing
the model or fitting any neural parameter.
"""
import argparse,json
from pathlib import Path
import numpy as np
import mujoco
from import_virtual_walker import FRAME,safe
def main():
 p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 root=Path(r'D:\GitHub\Bipedal_Robot');native=np.load(root/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/verification.npz');labels=native['labels'].tolist();values=native['data'];m=mujoco.MjModel.from_xml_path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\LiDeng_W2L_20261104\fresh_repeat\W2L_VirtualWalker_Newer\repeat\model.xml');d=mujoco.MjData(m)
 bodies=['foot_L','toe_L','foot_R','toe_R','Left hip wing','Right hip wing'];joints=['hip_L','knee_L','ankle_L','toe_L','hip_R','knee_R','ankle_R','toe_R'];samples=np.linspace(0,len(values)-1,120,dtype=int);rows=[]
 for sign in (1.,-1.):
  residuals=[];maxima=[]
  for k in samples:
   mujoco.mj_resetData(m,d);d.qpos[:3]=0.;d.qpos[3:7]=[1.,0.,0.,0.]
   for name in joints:d.qpos[int(m.joint(name).qposadr[0])]=sign*values[k,labels.index(name+'_angle')]
   mujoco.mj_forward(m,d);P=np.array([d.xpos[m.body(safe(name)).id]for name in bodies]);origin=values[k,[labels.index('root_'+axis)for axis in 'xyz']];Q=np.array([(values[k,[labels.index(name+'_'+axis)for axis in 'xyz']]-origin)@FRAME for name in bodies])
   U,_,V=np.linalg.svd(P.T@Q);correction=np.eye(3);correction[2,2]=np.linalg.det(U@V);rotation=U@correction@V;error=np.linalg.norm(P@rotation-Q,axis=1);residuals.append(float(np.sqrt(np.mean(error**2))));maxima.append(float(error.max()))
  rows.append(dict(native_joint_displacement_sign=sign,mean_registered_position_RMS_m=float(np.mean(residuals)),maximum_registered_position_error_m=max(maxima),sample_RMS_m=residuals))
 result=dict(comparison=rows,bodies=bodies,samples=len(samples),scope='Native recorded positions versus imported FK; free root rotation fitted, root translation removed; no neural/body simulation or repair applied');(a.output/'RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps([{k:v for k,v in row.items()if k!='sample_RMS_m'}for row in rows],indent=2))
if __name__=='__main__':main()
