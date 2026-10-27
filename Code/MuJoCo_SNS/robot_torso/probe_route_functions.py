import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import opensim as osim
import numpy as np
HERE=Path(__file__).resolve().parent
source=HERE.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim'
model=osim.Model(str(source)); state=model.initSystem()
point=next(p for i in range(model.getMuscles().getSize()) for p in [model.getMuscles().get(i).getGeometryPath().getPathPointSet().get(j) for j in range(model.getMuscles().get(i).getGeometryPath().getPathPointSet().getSize())] if p.getConcreteClassName()=='MovingPathPoint')
moving=osim.MovingPathPoint.safeDownCast(point)
print('MOVING_API', [x for x in dir(moving) if 'Function' in x or 'location' in x])
print('SIMMSPLINE_API',[x for x in dir(osim.SimmSpline) if 'Derivative' in x or 'Value' in x or 'Point' in x])
f=moving.get_x_location(); x=osim.Vector(1,0.); derivative=osim.StdVectorInt(); derivative.append(0)
print('FUNCTION_VALUE_DERIVATIVE',f.calcValue(x),f.calcDerivative(derivative,x))
foot=json.loads((HERE.parent/'spinal/reports_20261027_robot_control/foot_inventory.json').read_text())
for comp in foot['components']:
 t=np.array(comp['transform']); r=t[:9].reshape(3,3).T
 if comp['mass_properties_density_one']:
  raw=np.array(comp['mass_properties_density_one'][0]); center=r@raw[:3]+t[9:12]
  print(comp['name'],'COM',np.round(center,5),'volume',raw[3])
