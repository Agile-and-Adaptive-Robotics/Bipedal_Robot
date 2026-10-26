from pathlib import Path
import json,shutil,re,difflib,hashlib
from collections import Counter
R=Path(r'D:\Github\Bipedal_Robot');E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E';V=E/'review/figure61_correction'
def loc(file,text):
 lines=(E/'chapters'/file).read_text(encoding='utf-8').splitlines();i=next(i for i,l in enumerate(lines,1) if text in l);return i,lines[i-1]
chapter=(E/'chapters/60-sns.tex').read_text(encoding='utf-8');app=(E/'chapters/95-AppendixD.tex').read_text(encoding='utf-8')
assert r'\subsection{Block Programming Neural Pathways}' in chapter and 'Graphical Circuit Construction' not in chapter
assert '2layerImprove.png' not in chapter and 'connor_ben_shared_rg_original.png' in chapter
assert r'\textbf{(A)}' in chapter and r'\textbf{(B)}' in chapter
assert 'Years later, I split the shared RG into one RG for each side' in chapter
assert r'\label{sec:app_conversion_repairs}' in app
audit={};diff=''
for name in ('60-sns.tex','95-AppendixD.tex'):
 b=(V/('before_'+name)).read_text(encoding='utf-8');a=(E/'chapters'/name).read_text(encoding='utf-8')
 nums=lambda t:Counter(re.findall(r'\d+(?:\.\d+)?',t))
 labs=lambda t:set(re.findall(r'\\label\{([^}]*)\}',t))
 audit[name]={'labels_retained':labs(b)<=labs(a),'added_labels':sorted(labs(a)-labs(b)),'removed_numeric_tokens':list((nums(b)-nums(a)).elements()),'added_numeric_tokens':list((nums(a)-nums(b)).elements()),'sha256':hashlib.sha256((E/'chapters'/name).read_bytes()).hexdigest()}
 assert audit[name]['labels_retained']
 diff+=''.join(difflib.unified_diff(b.splitlines(True),a.splitlines(True),fromfile='before/'+name,tofile='after/'+name))
(V/'conservation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');(V/'changes.diff').write_text(diff,encoding='utf-8')
result=json.loads((E/'review/build_result.json').read_text());assert result['errors']==0 and result['undefined_references']==0 and result['missing_figures']==0
B=Path(result['build']);pdf=E/'Revision_E_corrected.pdf';shutil.copy2(B/'main.pdf',pdf)
(E/'review/current_output.json').write_text(json.dumps({'pdf':pdf.name,'build':result['build'],'reason':'Corrected original Figure 6.1 and conversion-repair appendix. Earlier Revision_E_local.pdf is open/locked and retained.'},indent=2),encoding='utf-8')
directives=[('Editable source','60-sns.tex',r'\chapter{Synthetic','The complete editable chapter and supporting TeX project remain in Revision_E.'),('Original versus later controller','60-sns.tex','Years later, I split','Original shared RG and Ben’s later per-side RG/feedback extension are explicitly separated.'),('Figure 6.1 panels','60-sns.tex',r'\caption[Original AnimatLab','A and B labels added in TeX; B replaced with Ben’s exact attachment; old image retained but unused in Figure 6.1.'),('Section heading','60-sns.tex',r'\subsection{Block Programming','Heading changed verbatim to Block Programming Neural Pathways.'),('MyoConverter summary','60-sns.tex','My first conversion','Retained readable main-text problems and added appendix cross-reference.'),('MyoConverter details','95-AppendixD.tex',r'\section{OpenSim-to-MuJoCo','Consolidated source-supported repairs in Appendix D.1; existing Appendix A force comparisons and Appendix C robot geometry details remain.')]
rows=[]
for directive,file,anchor,treatment in directives:
 line,text=loc(file,anchor);rows.append(dict(directive=directive,file=file,line=line,text=text,treatment=treatment,acceptance='Author review pending'))
(V/'compliance.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
source=R/'Code/MuJoCo_SNS/spinal/runner.py';lines=source.read_text(encoding='utf-8').splitlines()
evidence=[]
for key in ['hip joint axes are flipped','prune the tiny','ACTIVE Fmax = 1','repair 3: patella','pathpoint marker spheres','range-vs-equality','All 14 driver hinges','19 explicit ground-vs-body']:
 line=next(i for i,l in enumerate(lines,1) if key in l);evidence.append(dict(file=str(source.relative_to(R)),line=line,code=lines[line-1]))
(V/'conversion_source_evidence.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
report=['# Figure 6.1, controller history, and conversion details','',
'Ben corrected the historical identity of Figure 6.1(B). The prior draft incorrectly used the later controller hierarchy beside a description of the original Connor–Ben model. The exact supplied attachment is now the B asset, preserved byte-for-byte. Labels (A) and (B) are added in LaTeX. No circuit connections or raster content were redrawn.',
'','The original model has a shared four-neuron RG with bilateral PF groups and MN-level afferent feedback. Years later Ben split the RG into one per side and added RG/PF sensory feedback and foot-contact feedback. The text and caption distinguish those stages explicitly. The claim that Ben constructed the right-side subsystems only during that later extension was removed.',
'','The main chapter already contained specific MyoConverter problems. Before this correction, Appendix A also discussed force comparisons and rectus-femoris routing, and Appendix C discussed robot geometry; there was no single consolidated repair account in Appendix D. New §D.1 supplies it: hip-direction corrections, retained knee convention, six disabled entries within 92, restoration of eight active-force scales, coupled-path limits versus 14 anatomical driver hinges, rectus-femoris patella routing, and the 19-to-six contact-pair restriction. These numerical additions are traceable to the implementation in conversion_source_evidence.json; no new simulation measurements were generated.',
'','The main text was also corrected to say that anatomical joint limits were re-enabled separately from coupled path-point limits. Passive torques were not presented as the sole replacement for those limits. The driver-joint and follower-joint distinction prevents an ambiguous account of “removing joint limits.”',
'','## Per-directive evidence','']
for r in rows:report.append(f"- **{r['directive']}**, `chapters/{r['file']}:{r['line']}`: {r['treatment']}")
report+=['','## Build and preservation','',f"Fresh build `{B.name}`: {result['pages'][0]} pages, zero errors, zero unresolved references/citations, and zero missing figures. The 26 document-wide overflow warnings remain. No oversized-float warning was produced. Rendered Figure 6.1 (printed p.121) and new Appendix D.1 (printed p.229 onward) were checked. Output: `{pdf.name}`. The earlier Revision_E_local.pdf was open and locked, so it was retained and the correction saved under this new filename.",
'','Exact prior sources are preserved in before_60-sns.tex and before_95-AppendixD.tex. changes.diff and conservation.json record the changes. Old figure assets and all original labels remain. Diagram density and the supplied raster resolution still limit the smallest text; the attachment itself has not been altered. Overleaf, repository model code, and remote machine files were not edited. No commit or push was made.']
(V/'CORRECTIONS.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
# Keep existing current-line ledgers usable without changing their historical comment evidence or scores.
for ledger in [E/'review/clarity_pass/comment_compliance.json',E/'review/self_scoring_matrix.json']:
 data=json.loads(ledger.read_text(encoding='utf-8'))
 for row in data:
  f=row.get('file','').replace('chapters/','');f=f or '60-sns.tex'
  if not (E/'chapters'/f).exists():continue
  a=(E/'chapters'/f).read_text(encoding='utf-8').splitlines();old=row.get('text','');matches=[(i,l) for i,l in enumerate(a,1) if old[:65] in l] if old else []
  if row.get('id')=='M-105':matches=[(i,l) for i,l in enumerate(a,1) if l.startswith('Years later, I split')]
  if matches:row['line'],row['text']=matches[0]
 ledger.write_text(json.dumps(data,indent=2),encoding='utf-8')
header='''# Latest correction — original Figure 6.1 and conversion-repair appendix

Current output: `Revision_E_corrected.pdf` (270 pages). The earlier Revision_E_local.pdf was open and locked and remains the earlier version. Figure 6.1 now has (A)/(B) labels and uses Ben’s exact supplied original Connor–Ben controller as panel B. The main text distinguishes the original shared RG from Ben’s later per-side RG extension with RG/PF and contact feedback. The heading is “Block Programming Neural Pathways.” Main-text conversion details now link to a consolidated Appendix D.1 repair account, complementing the existing Appendix A force comparisons. See `review/figure61_correction/CORRECTIONS.md` for source lines and build evidence. Earlier pass reports below are preserved as history.

---

'''
p=E/'REPORT.md';old=p.read_text(encoding='utf-8');p.write_text(header+old,encoding='utf-8')
p=R/'CHATGPT_REPORT.md';old=p.read_text(encoding='utf-8');p.write_text('# Figure 6.1 and MyoConverter detail correction — 26 October 2026\n\nRevision E now contains Ben’s supplied original Connor–Ben shared-RG diagram as labeled Figure 6.1(B), with the mechanical walker labeled (A). Later per-side RG and RG/PF/contact-feedback changes are attributed to Ben and distinguished from the original model. The heading is exactly “Block Programming Neural Pathways.” New Appendix D.1 consolidates the conversion repairs; the main summary and existing Appendix A force comparisons remain. Details and line evidence: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/review/figure61_correction/CORRECTIONS.md`. Editable chapter and appendix sources remain in Revision_E/chapters. No Overleaf edits, commits, or pushes.\n\n---\n\n'+old,encoding='utf-8')
p=R/'ZCODE_HANDOFF_20261023.md';old=p.read_text(encoding='utf-8');p.write_text('# Author correction — original controller identity\n\nFigure 6.1(B) must remain Ben’s supplied original Connor–Ben shared-RG controller (`connor_ben_shared_rg_original.png`). The per-side RG split, RG/PF feedback, and foot-contact feedback were Ben’s later extensions. Do not replace the original diagram with `2layerImprove.png` again. The subsection heading is “Block Programming Neural Pathways.” Conversion repairs are consolidated in E Appendix D.1, with a main-text cross-reference. See E/review/figure61_correction/CORRECTIONS.md.\n\n---\n\n'+old,encoding='utf-8')
print(json.dumps({'pdf':str(pdf),'build':result,'directives':rows},indent=2))
