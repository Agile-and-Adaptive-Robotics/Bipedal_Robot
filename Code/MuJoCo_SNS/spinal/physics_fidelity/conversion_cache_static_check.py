"""Check old conversion samples against fresh source static equilibrium."""
import os,sys,pickle,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,opensim as osim
ROOT=Path(__file__).resolve().parents[4];source=ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'
m=osim.Model(str(source));m.setUseVisualizer(False);s=m.initSystem();cs=m.getCoordinateSet();muscles=m.getMuscles();rows={}
for name in ('bifemsh_r','soleus_r'):
 cache=pickle.loads((source.parent/'Step3_muscleKinetics'/(name+'.pkl')).read_bytes());points=[]
 for values,cached in zip(cache['jit_list_set'],np.asarray(cache['mtu_force_osim'])[0]):
  for k in range(cs.getSize()):
   c=cs.get(k)
   if not c.getLocked(s):c.setValue(s,c.getDefaultValue(),False)
  for joint,value in zip(cache['jit_uniq'],values):cs.get(joint).setValue(s,float(value),False)
  for k in range(muscles.getSize()):muscles.get(k).setActivation(s,1.)
  try:
   m.equilibrateMuscles(s);m.realizeDynamics(s);fresh=float(muscles.get(name).getTendonForce(s))
   points.append(dict(coordinates=dict(zip(cache['jit_uniq'],values)),cached_tension_N=float(cached),fresh_tension_N=fresh,difference_N=fresh-float(cached)))
  except Exception as ex:points.append(dict(error=str(ex)))
 rows[name]=dict(points=points,max_absolute_difference_N=max(abs(p['difference_N']) for p in points if 'difference_N'in p))
out=Path(__file__).resolve().parent/'runs/audit_20261030/CONVERSION_CACHE_STATIC_CHECK.json';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
