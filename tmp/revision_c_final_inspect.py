from pathlib import Path
import re,json,shutil
R=Path(r'D:\Github\Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';C=D/'Notes/advisor_revision_20261022/Revision_C'
S=C/'review/final_sweep_20261023';S.mkdir(exist_ok=True)
for p in (C/'chapters').glob('*.tex'):
    q=S/'before'/p.name;q.parent.mkdir(exist_ok=True);shutil.copy2(p,q)
shutil.copy2(R/'CHATGPT_REPORT.md',S/'CHATGPT_REPORT_before.md')
text=(R/'CHATGPT_REPORT.md').read_text(encoding='utf-8-sig');entries=[];file=''
for block in re.split(r'(?=^### )',text,flags=re.M):
    m=re.match(r'### ([A-Z]+-\d+)\s*·\s*line (\d+)\s*·\s*(comment|tracked edit)',block)
    if not m:continue
    ident,ln,kind=m.groups();quotes=[]
    for q in re.finditer(r'\*\*(Advisor[^\n]*|Pending advisor edit[^\n]*|Ben.s reply[^\n]*)\*\*:?\s*\n(.*?)(?=\n\*\*|\n###|\n## |\Z)',block,re.S):
        quotes.append({'speaker':q.group(1),'text':q.group(2).strip()})
    entries.append({'id':ident,'line':int(ln),'kind':kind,'quotes':quotes,'raw':block})
(S/'recorded_comments.json').write_text(json.dumps(entries,indent=2),encoding='utf-8')
print('LEDGER',len(entries),'comments',sum(x['kind']=='comment' for x in entries),'edits',sum(x['kind']=='tracked edit' for x in entries))
for e in entries:
    if e['kind']=='comment':
        q=next((q['text'] for q in e['quotes'] if q['speaker'].startswith('Advisor')),'')
        print(e['id'],q.replace('\n',' ')[:230])
B=Path((C/'review/latest_build.txt').read_text());aux=(B/'chapters/15-background.aux').read_text()
print('FIGURES 2.8/2.9',*[m.group(0)[:180] for m in re.finditer(r'\\newlabel\{[^}]+\}\{\{2\.[89]\}[^\n]+',aux)])
for name in ['15-background.tex','30-torque.tex']:
    s=(C/'chapters'/name).read_text(encoding='utf-8-sig')
    print('\nBLOCKS',name)
    for m in re.finditer(r'\\begin\{(table|figure|sidewaysfigure)\}.*?\\end\{\1\}',s,re.S):
        label=re.search(r'\\label\{([^}]+)\}',m.group());images=re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',m.group())
        print('line',s[:m.start()].count('\n')+1,m.group(1),label.group(1) if label else '?','chars',len(m.group()),'images',images)
