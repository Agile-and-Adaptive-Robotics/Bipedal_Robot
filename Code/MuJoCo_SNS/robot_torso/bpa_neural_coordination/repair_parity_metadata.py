"""Repair name/count diagnostics only; never advance or change saved neurons."""
import os,json
from pathlib import Path
import gait_family_parity as P
import numpy as np,mujoco

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\family_parity')

def mirror(n):
 for a,b in (('_r','_l'),('_l','_r'),('_R','_L'),('_L','_R')):
  if n.endswith(a):return n[:-2]+b
 return n

def main():
 model=P.REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml';m=mujoco.MjModel.from_xml_path(str(model));actors=[m.actuator(i).name for i in range(m.nu)];changes=[]
 for family in ('s3k','w2lvar','syn6'):
  os.environ['AARL_NET']=family;n=P.B.build(actors,dt=.002,interleg=True);names=[None]*len(n.compiled.V)
  for label,index in n.idx.items():
   if names[index] is None:names[index]=label
  path=ROOT/(family+'_parity.npz');z=np.load(path);data={k:z[k] for k in z.files};z.close();before=data['neuron_names'].tolist()
  assert before[:len(names)]==names and data['reference'].shape[1]==len(names)
  aliases={label:int(index) for label,index in n.idx.items() if label!=names[index]}
  data['neuron_names']=np.array(names);data['logging_aliases_json']=json.dumps(aliases)
  numeric={k:v.copy() for k,v in data.items() if k in ('reference','fast','inputs')}
  np.savez_compressed(path,**data)
  with np.load(path) as check:
   for key,value in numeric.items():assert np.array_equal(value,check[key])
  matrix=n.compiled.g_max_non;lookup={s:i for i,s in enumerate(names)};seen=set();unequal=[]
  for post,pre in zip(*np.nonzero(matrix)):
   a,b=names[pre],names[post];ma,mb=mirror(a),mirror(b)
   if (ma,mb)==(a,b) or ma not in lookup or mb not in lookup:continue
   mp,mq=lookup[mb],lookup[ma];pair=tuple(sorted(((int(post),int(pre)),(int(mp),int(mq)))))
   if pair in seen:continue
   seen.add(pair);g,h=float(matrix[post,pre]),float(matrix[mp,mq])
   if abs(g-h)>1e-12:unequal.append(dict(source=a,target=b,conductance=g,mirror_conductance=h))
  meta=ROOT/(family+'_parity.json');r=json.loads(meta.read_text());previous=r.get('legacy_weight_audit',dict(count=r['mirrored_nonzero_weight_differences'],examples=r['first_weight_differences']))
  r.update(legacy_weight_audit=previous,mirrored_nonzero_weight_differences=len(unequal),first_weight_differences=unequal[:12],weight_audit_definition='Unequal mirrored pairs counted once, including one-sided zero pathways; both suffix directions swapped.',logging_aliases=aliases,neuron_name_count=len(names),metadata_repair_scope='Canonical names and mirror diagnostics only; all saved numerical inputs/reference/fast traces unchanged. No neural or physical integration rerun.')
  meta.write_text(json.dumps(r,indent=2),encoding='utf-8');changes.append(dict(family=family,old_name_count=len(before),canonical_name_count=len(names),logging_aliases=aliases,legacy_unequal_count=previous['count'],corrected_unequal_pair_count=len(unequal)))
 (ROOT/'METADATA_REPAIR.json').write_text(json.dumps(changes,indent=2),encoding='utf-8');print(json.dumps(changes,indent=2))

if __name__=='__main__':main()
