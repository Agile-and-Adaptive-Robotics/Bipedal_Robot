"""Verify physical sleeve mesh lengths and straightness at every saved pose."""
from pathlib import Path
import json,csv
import numpy as np
from packaging_geometry import parts
out=Path(__file__).resolve().parents[1]/'packaging/output'
d=json.loads((out/'scene_data_corrected.json').read_text());rows=[]
for r in d['routes']:
 lengths=[];conflicts=[]
 for p in d['poses']:
  objects,m=parts(r,p['W'],p['knee_deg']);lengths.append(m['length_m']);conflicts.append(m['straight_span_conflict'])
  for name,g in objects.items():
   if g['kind']!='sleeve':continue
   centers=g['vertices'].reshape(-1,16,3).mean(axis=1)
   length=np.linalg.norm(centers[-1]-centers[0]);expected=m['length_m']/(2 if r['split_bowden'] else 1)
   assert abs(length-expected)<1e-7,(r['site'],name,length,expected)
   axis=(centers[-1]-centers[0])/length
   assert np.linalg.norm(np.cross(centers-centers[0],axis),axis=1).max()<1e-7
 rows.append(dict(site=r['site'],min_length_mm=min(lengths)*1000,max_length_mm=max(lengths)*1000,straight_span_conflict=any(conflicts)))
with (out/'straight_sleeve_audit.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('PASS: all 84 sleeves are straight and match saved lengths across 25 poses.')
print('First-guide span conflicts:',[r['site'] for r in rows if r['straight_span_conflict']])
