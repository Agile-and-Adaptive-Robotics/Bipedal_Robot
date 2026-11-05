"""09_BA_003's resolved link holes and named Theta1 coordinate system."""
from pathlib import Path
import json,numpy as np,math
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'09_BA_003'
ASSEMBLY=json.loads((DATA/'assembly_geometry.json').read_text())
DATUMS=json.loads((DATA/'datums.json').read_text())
T1=np.array(DATUMS['Theta1 Coordinate System'][9:12])
def holes(prefix):
 row=next(r for r in ASSEMBLY['components'] if r['name'].startswith(prefix))
 cs=np.array([c['center_assembly'] for c in row['cylinders'] if abs(c['local'][6]-.00315)<1e-6])
 cs=np.unique(np.round(cs[:,:2],10),axis=0)
 assert len(cs)==2
 return cs[np.argmax(cs[:,1])],cs[np.argmin(cs[:,1])]
A,C0=holes('04_06_FL_001');B,D0=holes('04_05_BL_001')
assert np.linalg.norm(D0-T1[:2])<1e-8
LF=np.linalg.norm(C0-A);LB=np.linalg.norm(D0-B);V=C0-D0
def rot(q):return np.array([[math.cos(q),-math.sin(q)],[math.sin(q),math.cos(q)]])
def candidates(q):
 v=rot(q)@V;E=A-v;delta=E-B;d=np.linalg.norm(delta);u=delta/d
 a=(LB*LB-LF*LF+d*d)/(2*d);h2=LB*LB-a*a
 if h2 < -1e-10:raise ValueError('CAD linkage cannot close at this angle')
 h=math.sqrt(max(0,h2));base=B+a*u;perp=np.array([-u[1],u[0]])
 return [base+h*perp,base-h*perp]
BRANCH=int(np.argmin([np.linalg.norm(p-D0) for p in candidates(0)]))
def pose(q):
 Q=rot(q);D=candidates(q)[BRANCH];C=D+Q@V
 alpha=np.linalg.solve(np.c_[C-A,-(D-B)],B-A)[0];I=A+alpha*(C-A)
 return {'A':A,'B':B,'C':C,'D':D,'I':I,'Q':Q,'q':q}
if __name__=='__main__':
 errors=[]
 for q in np.linspace(np.deg2rad(-120),np.deg2rad(10),101):
  p=pose(q);errors.extend([abs(np.linalg.norm(p['C']-A)-LF),abs(np.linalg.norm(p['D']-B)-LB),abs(np.linalg.norm(p['C']-p['D'])-np.linalg.norm(V))])
 assert max(errors)<1e-10
 print('09_BA_003 link lengths mm:',LF*1000,LB*1000,np.linalg.norm(V)*1000,'closure residual m:',max(errors),'t1 coincidence m:',np.linalg.norm(D0-T1[:2]))
