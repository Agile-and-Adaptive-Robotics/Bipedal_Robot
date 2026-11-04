"""Coarse OpenSim length/MIF screening; not routing optimization."""
from pathlib import Path
import itertools,json,math,csv
import opensim as osim
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parents[1]/'packaging/hardware'
OUT.mkdir(parents=True,exist_ok=True)
mapping=json.loads((ROOT/'Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json').read_text())
files=list((ROOT/'Solid_Models/OpenSim').rglob('gait2392_simbody.osim'))
assert files
osim.Logger.setLevelString('error')
m=osim.Model(str(files[0]));s=m.initSystem();muscles=m.getMuscles();coords=m.getCoordinateSet()
def fmax(L,d):
 a,b=(.4895,.03068) if d==10 else (1.4877,.0248)
 return 620*a*math.atan(b*(L-.0075)*620)
rows=[]
for item in mapping['actuators']:
 names=item['human_group'];names=[n for n in names if muscles.contains(n)]
 assert names,item
 dofs=[item['primary_dof']]+item['secondary_dofs']
 ranges=[]
 for dof in dofs:
  bounds=mapping['rom_table_deg'][dof]
  if 'knee' in dof:bounds=[-120,0]
  ranges.append([bounds[0]+(bounds[1]-bounds[0])*i/4 for i in range(5)])
 maxlen=0
 for angles in itertools.product(*ranges):
  for i in range(coords.getSize()):
   c=coords.get(i);c.setValue(s,c.getDefaultValue(),False)
  for dof,angle in zip(dofs,angles):coords.get(dof).setValue(s,math.radians(angle),False)
  m.realizePosition(s)
  maxlen=max(maxlen,*[muscles.get(n).getLength(s) for n in names])
 mif=sum(muscles.get(n).getMaxIsometricForce() for n in names)
 f10=fmax(maxlen,10);f20=fmax(maxlen,20)
 rows.append(dict(id=item['id'],group=item['master_name'],human_group=';'.join(names),mif_N=mif,coarse_Lmax_m=maxlen,Fmax10_N=f10,Fmax20_N=f20,ten_mm_screen_pass=mif<=.5*f10,initial_n20_half_Fmax=math.ceil(mif/(.5*f20)),sampled_dofs=';'.join(dofs)))
with (OUT/'bpa_force_screen.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(OUT/'force_screen_provenance.json').write_text(json.dumps(dict(model=str(files[0]),grid_points_per_coordinate=5,meaning='Coarse sampled maximum path length, not a certified maximum. Existing map; gluteal alternatives not adopted. No transmission losses or torque matching in this initial force screen.',rows=rows),indent=2))
for r in rows:
 if r['group'] in ['Gracilis','Vastus Intermedius','Semimembranosus','Gluteus Maximus','Adductor Magnus']:
  print(r)
