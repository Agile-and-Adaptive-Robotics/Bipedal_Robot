from pathlib import Path
import numpy as np,json,hashlib
P=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261023\Revision_E\review\easteregg2_followup')
def counts(t,f):
 s=t[1:][(f[1:]>20)&(f[:-1]<=20)];a=[]
 for v in s:
  if not a or v-a[-1]>.45:a.append(float(v))
 return len(a)
r={}
for p in sorted(P.glob('rawlog*.npz')):
 with np.load(p,allow_pickle=False) as n:
  t=n['t'];h=n['heelF'];q=n['q']
  r[p.stem]={'keys':n.files,'heel_counts_full_window':[counts(t,h[:,i]) for i in range(2)],'knee_ranges_deg':[float(np.ptp(q[:,i])) for i in (1,4)],'tilt_max_deg':float(np.max(n['tilt'])),'t_end':float(t[-1])}
j=json.loads((P/'changed_body_TRUEWINNER_w25.json').read_text());rows=j['rows']
r['mass_changes']={'count':len(rows),'survive_sum':sum(x['survive'] for x in rows),'tilt_over_45':sum(x['tilt']>45 for x in rows),'tilt_range':[min(x['tilt'] for x in rows),max(x['tilt'] for x in rows)],'minimum_pelvis_height':min(x['pelz'] for x in rows),'cycle_range':[min(min(x['cycL'],x['cycR']) for x in rows),max(max(x['cycL'],x['cycR']) for x in rows)]}
r['snapshot_manifest']=[{'path':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(P.iterdir()) if p.is_file() and p.name!='independent_audit.json']
(P/'independent_audit.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk!='keys'} if isinstance(v,dict) else v for k,v in r.items() if k!='snapshot_manifest'},indent=2))
