import os,sys,pickle
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import mujoco,numpy as np
root=Path(r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\Step3_muscleKinetics')
for name in ('bifemsh_r','soleus_r'):
 d=pickle.loads((root/(name+'.pkl')).read_bytes());print(name)
 for k,v in d.items():
  if isinstance(v,dict):print(k,{kk:(np.shape(vv),str(vv)[:100]) for kk,vv in v.items()})
  else:print(k,np.shape(v),str(v)[:140])
print(mujoco.mju_muscleGain.__doc__);print(mujoco.mju_muscleBias.__doc__)
