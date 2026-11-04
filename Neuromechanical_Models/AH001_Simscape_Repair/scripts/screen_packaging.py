"""Conservative neutral-pose capsule screening for saved BPA bundles.

No optimizer inputs are changed. Capsule end-rounding and polyline corners are
conservative preview envelopes. Signed gaps are capsule-to-capsule gaps, not
CAD-validated manufactured clearances.
"""
from pathlib import Path
import json,csv,math
import numpy as np
import trimesh
HERE=Path(__file__).resolve().parents[1]/'packaging/output'
D=json.loads((HERE/'scene_data.json').read_text()); W=D['W0']
def apply(T,p):
 T=np.asarray(T);return np.asarray(p)@T[:3,:3].T+T[:3,3]
def slice_path(P,start,end):
 out=[];s=0
 for A,B in zip(P[:-1],P[1:]):
  length=np.linalg.norm(B-A)
  if length<1e-10:continue
  a=max(start,s);b=min(end,s+length)
  if b>=a:
   if not out:out.append(A+(B-A)*(a-s)/length)
   out.append(A+(B-A)*(b-s)/length)
  s+=length
 return np.asarray(out)
def distance(a,b,c,d):
 # Exact closest points between two finite segments, including parallel cases.
 u=b-a;v=d-c;w=a-c;A=u@u;B=u@v;C=v@v;E=u@w;F=v@w
 cand=[];den=A*C-B*B
 if den>1e-18:
  s=(B*F-C*E)/den;t=(A*F-B*E)/den
  if 0<=s<=1 and 0<=t<=1:cand.append(np.linalg.norm(w+s*u-t*v))
 for p,x,y in ((a,c,d),(b,c,d),(c,a,b),(d,a,b)):
  z=y-x;t=np.clip((p-x)@z/max(z@z,1e-30),0,1);cand.append(np.linalg.norm(p-x-t*z))
 return min(cand)
bundles=[];route_audit=[]
scene=trimesh.load(HERE/'AH001_CAD.blend_import.glb',force='scene')
palette=[[18,120,98,255],[182,79,38,255],[91,74,163,255],[16,99,166,255],[171,51,99,255],[122,110,20,255],[46,125,140,255]]
def add_tube(P,start,end,radius,color,name):
 sleeve=slice_path(P,start,end)
 for k,(a,b) in enumerate(zip(sleeve[:-1],sleeve[1:])):
  if np.linalg.norm(b-a)<1e-8:continue
  mesh=trimesh.creation.cylinder(radius=radius,segment=np.array([a,b]),sections=16)
  mesh.visual.face_colors=color;scene.add_geometry(mesh,node_name=f'{name}_{k}',geom_name=f'{name}_{k}')
for ri,r in enumerate(D['routes']):
 i=int(np.argmin(np.abs(r['angles'][0])));j=int(np.argmin(np.abs(r['angles'][1])))
 L=r['rest']-r['contraction'][i][j];P=np.array([apply(W[b],p) for b,p in zip(r['bodies'],r['local'])])
 total=np.linalg.norm(np.diff(P,axis=0),axis=1).sum()
 minL=r['rest']*(1-{10:.17,20:.25,40:.32}[round(r['diameter']*1000)])
 maxD={10:.0225,20:.038,40:.075}[round(r['diameter']*1000)]
 piN=math.sqrt((r['rest']**2-minL**2)/(maxD**2-r['diameter']**2))
 dia=math.sqrt(r['rest']**2+(piN*r['diameter'])**2-np.clip(L,minL,r['rest'])**2)/piN
 for n,off in enumerate(r['offsets']):
  Q=P+np.asarray(W[r['offset_body']])[:3,:3]@off
  sleeve=slice_path(Q,r['fit'],r['fit']+L)
  bundles.append((r['site'],n,sleeve,dia/2,maxD/2))
  prefix=f"{r['site']}_BPA{n+1}"
  add_tube(Q,r['fit'],r['fit']+L,dia/2,palette[ri%len(palette)],prefix+'_sleeve')
  add_tube(Q,0,r['fit'],r['diameter']*.53,[85,90,95,255],prefix+'_origin_fitting')
  add_tube(Q,r['fit']+L,r['fit']*2+L,r['diameter']*.53,[85,90,95,255],prefix+'_distal_fitting')
  add_tube(Q,r['fit']*2+L,total,.0018,[50,60,65,255],prefix+'_tendon')
 route_audit.append(dict(site=r['site'],bpas=r['n'],gain=r['gain'],rest_mm=r['rest']*1000,current_mm=L*1000,diameter_mm=dia*1000,max_diameter_mm=maxD*1000,path_mm=total*1000,mount_length_mm=(L+2*r['fit'])*1000,sleeve_path_shortfall_mm=max(0,L+2*r['fit']-total)*1000,sleeve_corners=max(0,len(slice_path(P,r['fit'],r['fit']+L))-2),calibration_valid=bool(minL<=L<=r['rest'])))
rows=[]
for ai,(sa,na,PA,ra,ma) in enumerate(bundles):
 for sb,nb,PB,rb,mb in bundles[ai+1:]:
  if len(PA)<2 or len(PB)<2:continue
  # Broad-phase AABB with the maximum calibrated radii.
  if np.any(PA.max(axis=0)+ma<PB.min(axis=0)-mb) or np.any(PB.max(axis=0)+mb<PA.min(axis=0)-ma):continue
  dist=min(distance(a,b,c,d) for a,b in zip(PA[:-1],PA[1:]) for c,d in zip(PB[:-1],PB[1:]))
  gap=(dist-ra-rb)*1000;maxgap=(dist-ma-mb)*1000
  if maxgap<5:rows.append(dict(route_a=sa,bpa_a=na+1,route_b=sb,bpa_b=nb+1,current_capsule_gap_mm=gap,max_OD_capsule_gap_mm=maxgap))
for name,records in [('route_packaging_screen.csv',route_audit),('bundle_overlap_candidates.csv',sorted(rows,key=lambda r:r['max_OD_capsule_gap_mm']))]:
 with (HERE/name).open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(records[0]));wr.writeheader();wr.writerows(records)
summary=dict(bundle_count=len(bundles),max_OD_overlap_candidates=sum(r['max_OD_capsule_gap_mm']<0 for r in rows),current_overlap_candidates=sum(r['current_capsule_gap_mm']<0 for r in rows),sleeve_shortfall_routes=[r['site'] for r in route_audit if r['sleeve_path_shortfall_mm']>1e-6],sharp_sleeve_corner_routes=[r['site'] for r in route_audit if r['sleeve_corners']>0],scope='Neutral right-side only; preview attachment registration; capsule-to-capsule screening, not CAD surface/contact certification.')
(HERE/'screen_summary.json').write_text(json.dumps(summary,indent=2,default=int));print(json.dumps(summary,indent=2,default=int))
scene.export(HERE/'AH001_Packaging.glb')
