from pathlib import Path
import json,numpy as np
p=Path(__file__).resolve().parents[3];a=p/'Code/MuJoCo_SNS/spinal/reports_20261027_robot_control'
cs=json.loads((a/'robot_mass_geometry.json').read_text())['components']
pel=next(c for c in cs if '/02_01_PE_' in c['name']);R=np.array(pel['transform'])[:9].reshape(3,3)
for s in ['r','l']:
 f=next(c for c in cs if ('03_01_FB_R' if s=='r' else '03_03_FB_L') in c['name'])
 t=next(c for c in cs if ('05_01_TI_R' if s=='r' else '05_02_TI_') in c['name'])
 def axis(c):
  d=(R@np.array(c['transform'])[:9].reshape(3,3).T)[:,1]
  return d if d[1]<0 else -d
 u,v=axis(f),axis(t);angle=np.degrees(np.arctan2(np.cross(u,v)[2],u[:2]@v[:2]))
 print(s,'femur',u,'tibia',v,'neutral physical knee',angle)
