from pathlib import Path
import json
R=Path(r'D:\Github\Bipedal_Robot');p=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/review_matrix.json'
rows=json.loads(p.read_text(encoding='utf-8'))
out=[]
for r in rows:
 if r['kind']!='Comment thread':continue
 msgs=r.get('messages',[]); advisor=[m['text'] for m in msgs if m.get('author')=='ajh26']
 if r['file']=='60-sns.tex' or any(any(x in t.lower() for x in ['huh','confus','understand','awkward','gap in knowledge']) for t in advisor):
  out.append({k:r.get(k) for k in ['id','file','line','section','anchor_status','highlight','context','messages']})
(R/'tmp/clarity_catalog.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
for r in out:
 print(r['id'],r['file'],r['line'],r['section'])
 print('COMMENTS:', ' | '.join(m['text'] for m in r['messages'] if m['author']=='ajh26'))
 print('ANCHOR:',r['context'][:1100]);print()
