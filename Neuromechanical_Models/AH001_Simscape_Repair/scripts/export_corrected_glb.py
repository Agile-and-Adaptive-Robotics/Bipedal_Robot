"""Static import companion using the same geometry as the corrected .blend."""
from pathlib import Path
import json,base64,numpy as np,trimesh
from packaging_geometry import parts
P=Path(__file__).resolve().parents[1]/'packaging/output'
D=json.loads((P/'scene_data_corrected.json').read_text());s=trimesh.Scene();pose=D['poses'][0];W=pose['W']
for m in D['meshes']:
 v=np.frombuffer(base64.b64decode(m['v']),dtype='<f4').reshape(-1,3);f=np.frombuffer(base64.b64decode(m['f']),dtype='<u4').reshape(-1,3)
 t=trimesh.Trimesh(v,f,process=False);t.apply_transform(W[m['body']]);t.visual.face_colors=[110,125,135,255];s.add_geometry(t,node_name=m['name'])
palette=[[8,92,74,255],[168,61,26,255],[89,56,140,255],[13,84,145,255],[153,38,79,255],[120,102,15,255],[31,107,120,255]]
for i,r in enumerate(D['routes']):
 geo,_=parts(r,W,pose['knee_deg'])
 for name,g in geo.items():
  faces=[]
  for poly in g['faces']:
   faces.extend((poly[0],poly[k],poly[k+1]) for k in range(1,len(poly)-1))
  t=trimesh.Trimesh(g['vertices'],faces,process=False)
  t.visual.face_colors=palette[i%len(palette)] if g['kind']=='sleeve' else [20,25,30,255] if g['kind'] in ['cable','bowden'] else [100,110,115,255]
  s.add_geometry(t,node_name=r['site']+' / '+name)
result=s.export(P/'AH001_Packaging_Corrected.glb')
print('Exported corrected GLB:',len(s.geometry),'components')
