import json,numpy as np,base64
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'packaging/output'
d=json.loads((p/'scene_data.json').read_text());W={k:np.array(v) for k,v in d['W0'].items()}
def apply(T,v):return np.asarray(v)@T[:3,:3].T+T[:3,3]
for r in d['routes']:
 if r['site'] not in ['soleus_r','tib_ant_r','tib_post_r','flex_dig_r','ext_hal_r','med_gas_r']:continue
 C=np.array(r['neutral_chain']);F=W[r['offset_body']]@C[2]
 P=np.array([apply(W[b],v) for b,v in zip(r['bodies'],r['local'])])
 print(r['site'],'frame-origin',np.round(F[:3,3],5),'R',np.round(F[:3,:3],3).tolist(),'end',np.round(P[-1],5),'footlocalend',np.round(apply(np.linalg.inv(W['foot_r']),P[-1]),5))
print('CADfoot origin',W['foot_r'][:3,3])
for m in d['meshes']:
 if m['body']=='foot_r':
  v=np.frombuffer(base64.b64decode(m['v']),dtype='<f4').reshape(-1,3)
  print(m['name'],'local box',np.round([v.min(0),v.max(0)],5))
