import json, re, hashlib, csv
from pathlib import Path
from html.parser import HTMLParser
from datetime import datetime
from zoneinfo import ZoneInfo

ROOT=Path(r'D:\Github\Bipedal_Robot')
BASE=ROOT/'Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023'
SOURCE=BASE/'matrix_source'
entries=json.loads((BASE/'matrix_entries.json').read_text(encoding='utf-8'))
old=json.loads((ROOT/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/review/final_sweep_20261023/recorded_comments.json').read_text(encoding='utf-8'))

class Node:
    def __init__(self,tag='',attrs=()): self.tag=tag; self.attrs=dict(attrs); self.children=[]
    def text(self): return ''.join(c.text() if isinstance(c,Node) else c for c in self.children)
    def find(self,cls):
        out=[]
        if cls in self.attrs.get('class','').split(): out.append(self)
        for c in self.children:
            if isinstance(c,Node): out+=c.find(cls)
        return out
class Tree(HTMLParser):
    def __init__(self,s): super().__init__(convert_charrefs=True); self.root=Node(); self.stack=[self.root]; self.feed(s)
    def handle_starttag(self,t,a):
        n=Node(t,a); self.stack[-1].children.append(n)
        if t not in ['br','img','input','hr','meta','link']: self.stack.append(n)
    def handle_endtag(self,t):
        for i in range(len(self.stack)-1,0,-1):
            if self.stack[i].tag==t: self.stack=self.stack[:i]; break
    def handle_data(self,s): self.stack[-1].children.append(s)

def norm(s): return re.sub(r'\s+',' ',s).strip().strip('>').strip()
def esc(s): return s.replace('|','\\|').replace('\n','<br>')
def link(p,line=None): return f'[{p.name}{":"+str(line) if line else ""}](<{p.as_posix()}{":"+str(line) if line else ""}>)'
def sections(s,pos):
    names=[]; nums=[]; counters=[0]*4
    chapter={'10-introduction.tex':1,'15-background.tex':2,'30-torque.tex':4,'50-database.tex':5,'60-sns.tex':6}
    for m in re.finditer(r'\\(chapter|section|subsection|subsubsection)\*?\{([^}\n]*)\}',s[:pos+1]):
        level=['chapter','section','subsection','subsubsection'].index(m[1])
        if level==0: names=[m[2]]; counters=[0]*4
        else:
            counters[level]+=1
            for j in range(level+1,4): counters[j]=0
            names=names[:level]+[m[2]]
    return names,counters
def classify(t):
    t=t.lower()
    if 'overlapping' in t or 'figure is confusing' in t or 'motor center' in t: return 'Figure/layout','Check the actual rendered figure and the reader-facing purpose; no model experiment required.'
    if 'we?' in t or 'for each leg' in t or 'double check' in t: return 'Attribution/scope','Verify contributors and unilateral/bilateral scope from the source record; make the smallest wording correction.'
    if 'results? lessons?' in t or 'audit' in t: return 'Evidence/account','Locate the existing test and its result. State what was checked, found, and learned. Request a new test only if the supporting evidence is absent.'
    if 'libraries myself' in t: return 'Evidence/comparison','Identify the actual earlier libraries, compare functions and limitations; do not substitute an unrelated citation or performance claim.'
    if 'gap' in t or 'goal' in t or 'big picture' in t or 'organization' in t or 'awkward' in t or 'jumping' in t: return 'Substantial explanation','Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.'
    if 'connecting feedback' in t: return 'Future-work scope','Check existing schema links; specify the missing paper/pathway/review links without presenting finished links as future work.'
    return 'Clarification; scope to verify','Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.'

rows=[]; counts={}
for i,e in enumerate(entries,1):
    ident=f'M-{i:03d}'; tree=Tree(e['html']).root
    users=tree.find('review-panel-entry-user'); times=tree.find('review-panel-entry-time')
    author=users[0].text() if users else 'UNKNOWN'; date=times[0].text() if times else ''
    msgs=[]
    for c in tree.find('review-panel-comment'):
        u=c.find('review-panel-entry-user'); dt=c.find('review-panel-entry-time'); body=c.find('review-panel-comment-body')
        if body: msgs.append({'author':u[0].text() if u else '', 'date':dt[0].text() if dt else '', 'text':body[0].text()})
    clean='\n'.join(l for l in e['text'].splitlines() if l not in ['check','Resolve comment','more_vert','More options','show less','edit','add_circle','delete'])
    change=clean.split('\n',2)[2] if e['type']=='edit' else ''
    raw=(SOURCE/e['file']).read_text(encoding='utf-8')
    # CodeMirror offsets count JavaScript UTF-16 code units, not Python code points.
    pos16=int(e['pos']); prefix=raw.encode('utf-16-le')[:pos16*2].decode('utf-16-le',errors='ignore'); pos=len(prefix)
    line=prefix.count('\n')+1; names,nums=sections(raw,pos)
    ch={'10-introduction.tex':1,'15-background.tex':2,'30-torque.tex':4,'50-database.tex':5,'60-sns.tex':6}[e['file']]
    secnum='.'.join(map(str,[ch]+[n for n in nums[1:] if n]))
    section=secnum+' '+(' / '.join(names[1:]) or names[0] if names else 'Chapter')
    lines=raw.splitlines(); current=lines[min(line-1,len(lines)-1)] if lines else ''
    collapsed=pos>=len(raw.rstrip()) or (not current.strip() and pos>len(raw)-5)
    anchor='Collapsed at end of replaced file; current original location unavailable' if collapsed else 'Current review start mapped to frozen source; highlight end unavailable'
    highlight='Not exposed by current review UI; source line below is context, not an exact highlight.'
    historic=None
    if msgs:
        for o in old:
            if o['kind']=='comment' and any(norm(q.get('text',''))==norm(msgs[0]['text']) for q in o.get('quotes',[])):
                # Avoid matching ubiquitous short comments to the wrong old thread.
                if date not in o['raw'] or not o['id'].startswith({'10-introduction.tex':'INTRO','15-background.tex':'BG','50-database.tex':'DB','60-sns.tex':'SNS','30-torque.tex':'TORQUE'}[e['file']]): continue
                hh=re.search(r'\*\*Highlighted source:\*\*\s*```latex\n(.*?)\n```',o['raw'],re.S)
                if hh:
                    highlight='Historical exact highlight (22 Oct. capture): '+hh[1]; historic={'id':o['id'],'line':o['line'],'highlight':hh[1]}
                break
    category,action=classify(' '.join(m['text'] for m in msgs if m['author']=='ajh26')) if msgs else ('Tracked edit review','Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.')
    kind='Comment thread' if msgs else ('Advisor tracked edit' if author=='ajh26' else 'Ben tracked edit')
    row=dict(id=ident,file=e['file'],author=author,date=date,kind=kind,line=line,section=section,source= str(SOURCE/e['file']),anchor_status=anchor,collapsed=collapsed,context=current,highlight=highlight,historical=historic,messages=msgs,change=change,category=category,next_action=action,self_score=None,score_status='Unassessed: no new revision was made or tested in this pass',offset_utf16=pos16)
    rows.append(row)
    counts[kind]=counts.get(kind,0)+1

for r in rows:
    advisor_text=' '.join(m['text'] for m in r['messages'] if m['author']=='ajh26')
    if any(t in advisor_text for t in ['organization is still unclear','This all still feels a bit awkward.','This paragraph confuses me.']):
        r['self_score']=0
        r['score_status']='0/4 on the prior agent clarity revision: explicit repeated reader criticism; not a claim that every sentence is factually wrong'
    if "what a paper cites" in advisor_text:
        r['next_action']='Define both directions: references are works this paper cites; citing papers are works that cite it. Then verify which relationship the app displays. No new model experiment.'
    if 'Results? Lessons?' in advisor_text:
        r['next_action']='State the observed air-stepping result, the failure to achieve weight-bearing walking, and the specific next capability needed. Check existing traces before requesting any new simulation.'
    if 'What BPA interface?' in advisor_text:
        r['next_action']='Introduce inputs (simulated length and pressure), evaluated force model, output/coupling, and the valve-dynamics limitation in one local passage; verify this account against actual implementation.'
    if 'audit' in advisor_text.lower():
        r['next_action']='Replace generic audit language with actual checks and findings from artifacts. Hypothetical wrong origins or handedness in the advisor comment must not become claimed findings.'

stamp=datetime.now(ZoneInfo('America/Los_Angeles')).strftime('%Y-%m-%d %H:%M:%S Pacific')
intro=f'''# Advisor review matrix and execution contract

Updated {stamp}. This is the primary navigation and scoring view. Older captures and accountability findings remain verbatim below. **Overleaf is unchanged, and scheduled monitoring remains canceled.**

Ben's controlling instruction: section, source line, highlighted/changed text, author, and comment must be recorded; advisor edits without comments and Ben's edits count too. This matrix contains {len(rows)} current review entries: {counts}. All messages and replies within each thread are retained.

## How to use the references

Line links refer to **frozen read-only copies of the current Overleaf editor source** in Notes/advisor_review_20261023/matrix_source, not stale ProofFinal chapters or PDF pages. Section numbers are calculated from that snapshot's LaTeX headings. Pending text is included as shown by the editor; deleted wording remains in the exact change record. Source hashes are in matrix_manifest.json. Raw character offsets are retained only in JSON for reproducibility.

**Correction to my earlier explanation:** the values were review-card start offsets in JavaScript UTF-16 units, not reliable human navigation references. They give no highlight length. Some old threads collapsed to the end of a replaced chapter. Those rows explicitly say the original current location is unavailable; an end-of-file line is not presented as the original passage. Historical exact highlights are labeled historical. For other comments, the source line is context, not a fabricated exact selection. Exact before/after tracked edits are recorded independently of this limitation.

## Self-scoring rules, used before claiming an item is handled

Ben already explained the unclear prose and gaps. He does not need to restate them more forcefully. The obligation is on the editing agents to demonstrate compliance.

| Score | Required evidence |
|---|---|
| Unassessed | No scored revision/evidence review yet. Not a zero and not a pass. This capture does not award itself quality credit. |
| 0 | Revised passage still fails the specific request, contradicts evidence, or repeats the same documented defect. |
| 1 | Some wording changed, but the referent, gap, causal connection, attribution, or result remains unclear. |
| 2 | Direct answer is present, but context or factual verification is incomplete. |
| 3 | Comment intent, surrounding argument, attribution, and cited evidence checked; the local revision is reviewable. Reader acceptance remains unverified. |
| 4 | Ben or the advisor confirms that the concern is resolved. Agents cannot award themselves this score. |

Every future revision must populate: proposed/actual before and after, evidence path, exact checks, score, remaining limitation, and reviewer disposition. A build checks LaTeX; it does not establish clarity or model validity. A repeated substantive concern reopens the row and caps its prior editorial score at 1 pending reassessment. Do not deduct points from Ben for his edits. His edits and the advisor's edits provide examples for evaluating the agents' prose; authorship in Overleaf is an account label, not proof of which agent wrote it.

**Initial scored failures:** rows with “organization is still unclear”, “This all still feels a bit awkward”, and “This paragraph confuses me” receive 0/4 for the earlier agent clarity revision. This is scored from explicit reader feedback, not inferred from the mere presence of an unresolved thread. Other rows remain unassessed until their actual proposed revision and evidence are checked. This makes no claim of a whole-dissertation quality percentage.

Gap check for each affected section: (1) what was unknown or unavailable, (2) why it mattered, (3) what Ben did, (4) what observed result establishes, and (5) what remains unresolved. A list of software names or missing combinations does not pass. Keep knowledge gaps, engineering capability gaps, completed results, and proposed integration distinct.

## Work allocation and minimum additional evidence

Ben's preference is recorded: this chat specifies minor local edits and review criteria; GLM-5.3 and GLM-5.3-Flash handle additional testing, modeling, and optimization through his ZCode setup; the requested final polishing model is “6.1 Astra”, whose exact available identifier must be confirmed before any launch. No models or other chats were invoked in this pass, and no paid work was launched.

| Work | Smallest useful output | Gate before expansion |
|---|---|---|
| Definitions, pronouns, citation direction, paragraph order | One targeted replacement with before/after and preserved citations. No new simulation. | Reader can identify subject, contribution, and link to the argument. |
| “Results? Lessons?” / conversion “audit” | Existing experiment, specific observation, actual defect/fix if one existed, and limitation, with trace/log evidence. | Search existing records before requesting a rerun. The advisor's hypothetical wrong origins/handedness are not findings. |
| Accurate connectome transcription | Source-backed node and edge tables: source, target, sign, weight/convention, delay, cell model, source page/figure/table, and ambiguity. Overlay/render for visual checking plus machine-readable comparison against the built SNS network. | Resolve ambiguous edges; verify every edge and neuron, including omissions and extra edges. A matching census or a plausible diagram alone cannot pass. Never tune weights to hide uncertain wiring. |
| Reference neural replication | One source model and one published stimulus/deletion protocol first, matching equations, inputs, initial state, units, solver/time step, and scoring window; reference/replica traces and quantitative error. | Predefine observable metrics and tolerances justified by the source. Report transcription fidelity separately from model reduction or adaptation. Rhythm alone is not a complete replication. |
| Coupled neuromechanical replication | Original controller AND source mechanical/actuator/contact assumptions; baseline and the chosen perturbation/deletion result. | A controller port on different mechanics is an adaptation until the difference is characterized. Frozen/NaN runs fail. Support-assisted stepping is not unsupported walking. |
| Optimization | Only after wiring and baseline protocol pass; state which parameters may vary, objective, constraints, fixed controls, seeds, and independent validation. | Do not use optimization as the first response to a transcription or explanation failure. |
| Final polish | Only approved affected paragraphs plus their evidence matrix, one chapter at a time. Preserve numbers, labels, citations, and factual limits. | Fix evidence first; constrain context and output size. No automatic whole-dissertation rewrite or model switching. |

This queue is a proposed execution specification, not a claim that new experiments have been completed. The local replication README still labels its contents drafts; an existing verification script is not execution evidence. CROSS_PLATFORM_REPORT.md contains prior matched tests and explicit limits; those should be examined before commissioning duplicate runs. Replication model status needs a separate artifact audit before being upgraded.

## Matrix index

“Anchor at EOF” means the old location was lost in chapter replacement. Open the card for historical highlights and messages. All 137 rows have an author, type, source context or explicit lost-anchor status, and an unassessed score field.

| ID | Author / type | Section | Source line | Text / action | Self-score |
|---|---|---|---|---|---|
'''
for r in rows:
    summary=r['messages'][0]['text'] if r['messages'] else r['change']
    short=norm(summary)
    if len(short)>150: short=short[:147]+'...'
    loc='Anchor at EOF; see card' if r['collapsed'] else link(Path(r['source']),r['line'])
    score='0/4: repeated clarity failure' if r['self_score']==0 else 'Unassessed'
    intro+=f"| [{r['id']}](#{r['id'].lower()}) | {r['author']} / {r['kind']} | {esc(r['section']) if not r['collapsed'] else 'Original current section unavailable'} | {loc} | {esc(short)} | {score} |\n"

intro+='\n## Full review cards\n\n'
for r in rows:
    intro+=f"<a id=\"{r['id'].lower()}\"></a>\n\n### {r['id']} · {r['file']} · {r['kind']}\n\n"
    intro+=f"**Author/date:** {r['author']} · {r['date']}  \n"
    intro+=f"**Section:** {r['section'] if not r['collapsed'] else 'Original current section unavailable: old anchor collapsed at end of replaced file.'}  \n"
    intro+=f"**Frozen source line:** {link(Path(r['source']),r['line'])} ({'end-of-file anchor only, not original location' if r['collapsed'] else 'review start; not a confirmed highlight range'}).  \n"
    intro+=f"**Anchor confidence:** {r['anchor_status']}.\n\n"
    if r['change']:
        intro+=f"**Exact tracked change, including insertions/deletions without comments:**\n\n~~~text\n{r['change']}\n~~~\n\n"
    else:
        intro+=f"**Highlighted text:** {r['highlight']}\n\n"
        if r['historical']:
            intro+=f"**Earlier locator:** {r['historical']['id']}, source line {r['historical']['line']} in the 22 October capture; not the current source line.\n\n"
        for m in r['messages']:
            intro+=f"**{'Advisor' if m['author']=='ajh26' else 'Ben'} ({m['author']}) · {m['date']}:**\n\n~~~text\n{m['text']}\n~~~\n\n"
    if not r['collapsed']:
        intro+=f"**Current anchored source line (context, not asserted exact highlight):**\n\n~~~latex\n{r['context']}\n~~~\n\n"
    intro+=f"**Work class:** {r['category']}. **Next check:** {r['next_action']}\n\n**Self-score:** {r['score_status']}. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.\n\n"

intro+='\n---\n\n# Earlier snapshots and accountability record (preserved verbatim)\n\n'
report=ROOT/'CHATGPT_REPORT.md'
backup=BASE/'CHATGPT_REPORT_before_matrix.md'
previous=backup.read_text(encoding='utf-8') if backup.exists() else report.read_text(encoding='utf-8')
if not backup.exists(): backup.write_text(previous,encoding='utf-8')
report.write_text(intro+previous,encoding='utf-8')
(BASE/'review_matrix.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False),encoding='utf-8')
with (BASE/'review_matrix.csv').open('w',encoding='utf-8-sig',newline='') as f:
    fields=['id','file','author','date','kind','section','line','anchor_status','highlight','change','category','next_action','score_status']
    w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore'); w.writeheader(); w.writerows(rows)
manifest={'created':stamp,'entries':len(rows),'counts':counts,'lost_current_anchors':sum(r['collapsed'] for r in rows),'historical_highlights':sum(r['historical'] is not None for r in rows),'sources':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCE.glob('*.tex')},'overleaf_mutations':False,'scored_revision':False}
(BASE/'matrix_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
assert report.read_text(encoding='utf-8').endswith(previous)
assert len(rows)==len(entries)==137
assert all(r['change'] for r in rows if r['kind']!='Comment thread')
print(json.dumps(manifest,indent=2))
