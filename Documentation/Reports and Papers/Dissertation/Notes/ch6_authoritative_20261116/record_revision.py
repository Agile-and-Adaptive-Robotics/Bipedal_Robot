"""Record full source changes, numerical tokens, conservation, and commit paths."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import subprocess

HERE=Path(__file__).resolve().parent
DISS=HERE.parent.parent
ROOT=DISS.parents[2]
FINAL=DISS/'Final'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
number=re.compile(r'(?<![A-Za-z])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?')
ledger=[]
diffs=[]
content={}
for name in ['60-sns.tex','95-AppendixD.tex']:
    before=(HERE/'before'/name).read_text(encoding='utf-8-sig').splitlines()
    after=(FINAL/'chapters'/name).read_text(encoding='utf-8-sig').splitlines()
    diffs.extend(difflib.unified_diff(before,after,fromfile='before/'+name,tofile='Final/chapters/'+name,lineterm=''))
    for tag,i,j,k,l in difflib.SequenceMatcher(None,before,after,autojunk=False).get_opcodes():
        if tag=='equal':continue
        old='\n'.join(before[i:j]); new='\n'.join(after[k:l])
        on=number.findall(old); nn=number.findall(new)
        nums=[]
        for nt,oi,oj,ni,nj in difflib.SequenceMatcher(None,on,nn,autojunk=False).get_opcodes():
            if nt!='equal':nums.append({'operation':nt,'old':on[oi:oj],'new':nn[ni:nj]})
        ledger.append({'file':str(FINAL/'chapters'/name),'operation':tag,
            'old_start_line':i+1,'new_start_line':k+1,'old_text':old,'new_text':new,
            'all_old_numerical_tokens':on,'all_new_numerical_tokens':nn,'numerical_token_changes':nums,
            'interpretation_and_source_mapping':'NUMERICAL_CHANGE_LEDGER.md; AIR/PF/CAD/TIMING evidence reviews'})
    labels=lambda lines:set(re.findall(r'\\label\{([^}]+)\}','\n'.join(lines)))
    graphics=lambda lines:re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}','\n'.join(lines))
    content[name]={'before_sha256':sha(HERE/'before'/name),'after_sha256':sha(FINAL/'chapters'/name),
        'labels_removed':sorted(labels(before)-labels(after)),'labels_added':sorted(labels(after)-labels(before)),
        'before_graphics':graphics(before),'after_graphics':graphics(after)}
(HERE/'MANUSCRIPT.diff').write_text('\n'.join(diffs)+'\n',encoding='utf-8')
(HERE/'SOURCE_CHANGE_LEDGER.json').write_text(json.dumps(ledger,indent=2)+'\n',encoding='utf-8')
baseline=json.loads((HERE/'FINAL_BEFORE_HASHES.json').read_text(encoding='utf-8-sig'))
allowed={str(FINAL/'chapters'/n).casefold() for n in ['60-sns.tex','95-AppendixD.tex']}
changes=[]; unchanged=[]; missing=[]
for row in baseline:
    p=Path(row['Path'])
    if not p.exists():missing.append(str(p));continue
    if sha(p).lower()!=row['SHA256'].lower():changes.append(str(p))
    else:unchanged.append(str(p))
newfiles=[str(p) for p in FINAL.rglob('*') if p.is_file() and str(p) not in {x['Path'] for x in baseline}]
assert not missing and not newfiles
assert all(x.casefold() in allowed for x in changes)
conservation={'changed_final_files':changes,'unchanged_final_files':unchanged,
    'missing_final_files':missing,'unexpected_new_final_files':newfiles,'chapter_content_changes':content,
    'model_changes':'None. No bracket points, transforms, stiffness ordering, bounds, flags, force law, controller weights, or recordings changed.'}
(HERE/'CONTENT_CONSERVATION.json').write_text(json.dumps(conservation,indent=2)+'\n')
raw=subprocess.check_output(['git','-c','core.quotepath=false','status','--porcelain=v1','-z','-uall'],cwd=ROOT)
parts=raw.decode('utf-8').split('\0'); items=[];i=0
while i<len(parts):
    p=parts[i];i+=1
    if not p:continue
    code=p[:2];rel=p[3:]
    item={'status':code,'relative_path':rel,'absolute_path':str(ROOT/rel)}
    if 'R' in code or 'C' in code:
        item['original_relative_path']=parts[i];i+=1
    items.append(item)
resume_path=HERE/'RESUME_VERIFICATION.json'
resume=json.loads(resume_path.read_text()) if resume_path.exists() else None
# Both self-describing manifests will change on this completion pass.
for name in ['COMMIT_FILES.json','COMMIT_FILES.txt']:
    p=HERE/name;rel=p.relative_to(ROOT).as_posix()
    if not any(x['relative_path']==rel for x in items):
        tracked=bool(subprocess.check_output(['git','ls-files','--',rel],cwd=ROOT).strip())
        items.append({'status':' M' if tracked else '??','relative_path':rel,'absolute_path':str(p)})
items.sort(key=lambda x:x['relative_path'])
committed=[]
if resume:
    prefixes=['Documentation/Reports and Papers/Dissertation/Notes/ch6_authoritative_20261116/',
        'Documentation/Reports and Papers/Dissertation/Figures/60-sns/',
        'Documentation/Reports and Papers/Dissertation/Figures/95-AppendixD/']
    exact=['CHATGPT_REPORT.md','Documentation/Reports and Papers/Dissertation/Final/chapters/60-sns.tex',
        'Documentation/Reports and Papers/Dissertation/Final/chapters/95-AppendixD.tex']
    committed=[x for x in resume['committed_revision_paths']
        if x['relative_path'] in exact or any(x['relative_path'].startswith(k) for k in prefixes)]
payload={'preexisting_changes':[],
    'baseline':'Initial revision and resumed working trees were both clean before agent writes.',
    'resume_head':resume['head_at_resume'] if resume else None,
    'already_committed_revision_files':committed,'session_changes':items,
    'scope':'session_changes are the remaining completion files to commit now; earlier sources/figures already in resume_head.'}
(HERE/'COMMIT_FILES.json').write_text(json.dumps(payload,indent=2)+'\n')
header='Pre-existing uncommitted changes at initial start/resume: NONE.\nNo agent commit or push performed.\n'
if resume:header+='\nEarlier revision already present in commit '+resume['head_at_resume']+'.\n'
body='\nCOMPLETION FILES TO COMMIT NOW:\n\n'+'\n'.join(x['status']+' '+x['absolute_path'] for x in items)
if committed:body+='\n\nEARLIER REVISION VERSIONS ALREADY COMMITTED (any completion changes listed above still require commit):\n\n'+'\n'.join(x['absolute_path'] for x in committed)
(HERE/'COMMIT_FILES.txt').write_text(header+body+'\n')
(HERE/'PREEXISTING_STATUS.txt').write_text('Initial git status --short was empty, before session file creation.\nNo pre-existing working-tree changes were present.\n')
print(json.dumps({'source_hunks':len(ledger),'changed_final_files':len(changes),
    'preserved_final_files':len(unchanged),'session_commit_paths':len(items)}))
