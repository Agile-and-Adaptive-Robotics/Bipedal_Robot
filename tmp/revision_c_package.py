from pathlib import Path
import re,json,shutil,hashlib,subprocess,fitz
R=Path(r'D:\Github\Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';P=D/'Notes/advisor_revision_20261022';C=P/'Revision_C';B=Path((C/'review/latest_build.txt').read_text())
log=(B/'main.log').read_text(errors='replace')
for pat in [r'LaTeX Warning:.*undefined',r'There were undefined',r'multiply defined',r'Float too large',r'Label\(s\) may have changed',r'^!']:
    assert not re.search(pat,log,re.M),pat
assert (B/'pass_5.log').exists()
shutil.copy2(B/'main.pdf',C/'Revision_C.pdf')
shutil.copy2(B/'main.log',C/'review/build.log');shutil.copy2(B/'main.toc',C/'review/build.toc')
for p in sorted(B.glob('pass_*.log')):shutil.copy2(p,C/'review'/p.name)
counts=json.loads((C/'review/word_counts.json').read_text());integrity=json.loads((C/'review/source_integrity.json').read_text())
pdfdoc=fitz.open(C/'Revision_C.pdf');pdfpages=len(pdfdoc);bookmarks=pdfdoc.get_toc()
chapterpage=lambda prefix: next(page for level,title,page in bookmarks if level==1 and title.startswith(prefix))
offset=chapterpage('Introduction')-1
sadbpages=f"{chapterpage('The Sensory')-offset}–{chapterpage('Synthetic Nervous')-offset-1}"
snspages=f"{chapterpage('Synthetic Nervous')-offset}–{chapterpage('Discussion')-offset-1}"
def line(file,phrase):
    for n,s in enumerate((C/file).read_text(encoding='utf-8-sig').splitlines(),1):
        if phrase in s:return n
    raise ValueError((file,phrase))
def evidence(file,phrase):
    n=line(file,phrase);return f'[{Path(file).name}:{n}](<{(C/file).as_posix()}:{n}>)'
chaptertexts={p.name:p.read_text(encoding='utf-8-sig') for p in (C/'chapters').glob('*.tex')}
figureledger=[]
for before in (C/'review/before').glob('*.tex'):
    old=before.read_text(encoding='utf-8-sig')
    for m in re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}',old,re.S):
        label=re.search(r'\\label\{([^}]+)\}',m.group())
        if not label:continue
        label=label.group(1);locations=[name for name,text in chaptertexts.items() if r'\label{'+label+'}' in text]
        figureledger.append({'label':label,'before':before.name,'after':locations,'status':'retained_or_moved' if locations else 'omitted_from_compiled_copy; original asset retained'})
(C/'review/figure_dispositions.json').write_text(json.dumps(figureledger,indent=2))
checklist=[
('Separate local copy; no Overleaf writes','main.tex',r'\include{chapters/95-AppendixD}'),
('SNS begins with AnimatLab attempt','chapters/60-sns.tex','The AnimatLab Bipedal Walker'),
('SNS has tools-and-verification structure','chapters/60-sns.tex','Requirements for the New Tools'),
('Detailed SNS results moved into appendix','chapters/95-AppendixD.tex','Tuning and Supported Evaluations'),
('SADb shortened; coverage/pathways in appendix','chapters/96-AppendixE.tex','Record Sets and Coverage'),
('Pulley axes/configuration caption corrected','chapters/94-AppendixC.tex','Predicted insertion-torque magnitude'),
('Updated cam geometry shown','chapters/30-torque.tex','extensor_cam_route.pdf'),
('New interconnected bilateral RG construction credited','chapters/60-sns.tex','I subsequently expanded'),
('Final torque results sourced from analysis scripts','chapters/30-torque.tex','test-specific analysis scripts'),
('Completed-study wording per latest instruction','chapters/30-torque.tex','The redesigned configurations meet or exceed'),
('Conclusion uses that drafting assumption','chapters/80-conclusion.tex','The isometric tests demonstrate'),
('Unsupported-walking claims scoped','chapters/60-sns.tex','no harness-free walking'),
]
rows='\n'.join('| '+a+' | '+evidence(f,p)+' |' for a,f,p in checklist)
statrows=[]
for n in ['10-introduction.tex','15-background.tex','30-torque.tex','50-database.tex','60-sns.tex','70-overall.tex','80-conclusion.tex']:
    v=counts[n];old=v['before'][0];new=v['after'][0];statrows.append(f'| {n} | {old:,} | {new:,} | {100*(old-new)/old:.0f}% |')
allbefore=sum(sum(v['before']) for v in counts.values());allafter=sum(sum(v['after']) for v in counts.values())
report=f'''# Revision C

An edited, independently compilable local dissertation copy. Started from the committed round-two chapter proposals at 39e3f707 and unchanged supporting ProofFinal sources. Original proposal chapters, ProofFinal sources, bibliography bytes, and protected chapter-3 figures passed the integrity check. This pass made no Overleaf writes or Git commits/pushes. A merge to 39de1ccc was observed during the work; the dissertation, companion torque-paper, workbook, and flexor-endpoint sources used here did not differ across those commits.

## Main editorial decisions

- SNS now follows the AnimatLab construction and Ben's interconnected left/right RG extension → unresolved coordination and verification questions → integration requirements → tools → representative verification checks. The main chapter has five figures and one capability table. Larger reference sweeps, pathway removals, spiking comparisons, and supported-walker evaluations are in Appendix D. GUI inventories and debugging chronology were removed from the compiled argument; original text and figure assets remain recoverable.
- SADb leads with the resource and distinguishes the paper-only application from the five-table working base. Coverage details, linked-pathway tables, and the citation neighborhood are in Appendix E. Actual application screenshots replace retired plots; the pivot is framed through LaTeX to remove its empty viewport and enlarge the matrix excerpt.
- Torque retains final knee-design results in the main chapter. Earlier multi-site batch diagnostics, whole-body route-builder history, and transmission surfaces are in Appendix C. The newer companion-paper EPS exports were used for the cam route and two transmission figures: the similarly named old PDFs were stale. The false “two pulley positions” caption is corrected to actual joint-angle, pressure, and transmission variables. The two surfaces occupy separate pages. The second flexor path is separately solved, as confirmed in the endpoint and prediction code; the old shared-route statement is corrected.
- Introduction, selected background explanations, Discussion, and Conclusion were aligned with those decisions. The erroneous cross-product moment-arm definition was removed. Captions were shortened without adding unmeasured claims.

## Latest instruction: completed-study wording

Ben reported that the extensor data arrived last night and flexor collection is in progress, and directed writing as if the intended result were established. The abstract outcome sentence is therefore retained, and the torque summary and Conclusion use completed-result wording. This is an explicit drafting assumption. No missing flexor values, error scores, or confidence intervals were invented. Before release, the final extensor/flexor comparison must substantiate the stated human-reference outcome. The spreadsheet-derived cam table and torque range were removed after Ben clarified that Excel rows 16 and 17 are sanity checks; final torque values and presentation plots come from the test-specific analysis scripts. See `review/07_Bilateral_and_torque_source.md`. The active bipulley optimization on easteregg2 is a new campaign, not a completed result represented by the older sweep figures.

## Length

Approximate narrative source-word counts, excluding captions and tables, using the campaign helper with escaped-percent and end-of-file fixes. LaTeX labels/arguments can contribute tokens, so these are comparative counts rather than publication word counts.

| Chapter | Round two | Revision C | Reduction |
|---|---:|---:|---:|
'''+ '\n'.join(statrows)+f'''

SNS plus new Appendix D has {counts['60-sns.tex']['after'][0]+counts['95-AppendixD.tex']['after'][0]:,} narrative tokens versus {counts['60-sns.tex']['before'][0]:,} in the old SNS chapter. SADb plus new Appendix E has {counts['50-database.tex']['after'][0]+counts['96-AppendixE.tex']['after'][0]:,} versus {counts['50-database.tex']['before'][0]:,}. Across changed files, including added appendices, narrative/caption/table tokens decrease from approximately {allbefore:,} to {allafter:,}. Movement to appendices was accompanied by selection and shortening.

## Per-directive check

| Directive | Source evidence |
|---|---|
{rows}

The protected isometric-force chapter is byte-identical. All {integrity['protected_figure_hashes_verified']} entries in the protection manifest (active/retired figures and chapter source) match their recorded hashes; the copied active figures also match. The bibliography is byte-identical. See `review/source_integrity.json`.

## Verification

- Full MiKTeX build: pdflatex → BibTeX → three pdflatex passes, all successful in a fresh temporary directory. No undefined references/citations, duplicate labels, changed-label rerun warning, or oversized-float warning in the final log.
- Compiled output: `Revision_C.pdf`, {pdfpages} physical PDF pages. SADb occupies printed pages {sadbpages}; SNS {snspages}. Supporting appendices remain part of the complete document.
- Rendered and inspected the changed SADb/SNS sections, new appendices, and mechanical figures. Page renders and overview sheets are in `review/rendered`; the manifest records the exact pages.
- Existing wide equations/tables in inherited Appendices A–C still produce overfull-box warnings. These pre-existing layout issues were not declared repaired; the new moved tables were fitted to the text width. This is an editorial revision with a successful build, not a claim that all inherited typesetting is publication-ready.
- Final torque results must use the test-specific analysis-script outputs. Excel rows 16 and 17 are sanity checks, not presentation data. The spreadsheet-derived table was removed. `review/torque_sheet_evidence.json` is retained only as a superseded inspection record and is explicitly marked accordingly; no reduction script was run during this editorial correction.

## Review and integration

Open `Revision_C.pdf` to read, or `main.tex` to compile this separate copy. `chapters` contains the edited TeX; `figs` contains existing assets and the three newly converted existing torque EPS exports. Preserve that folder layout.

`review/changes.diff` compares the revised chapters against the saved round-two baseline and includes the new appendix registrations. `review/figure_dispositions.json` lists retained, moved, and omitted figure labels. `review/numeric_change_review.json` contains every removed/replaced source line containing numerals, including relocated lines: it is an inspection ledger, not a claim that each number disappeared. Detailed section notes explain scientific corrections, explicit omissions, and moves. Do not apply the entire copy over an advisor-edited Overleaf project blindly; integrate the selected chapter/section changes and required assets against the current live source.

The new data campaign, the expanded SADb app export, unsupported walking, and hardware/network deployment were not completed or evaluated by this editorial pass.
'''
(C/'REPORT.md').write_text(report,encoding='utf-8')
(C/'README.md').write_text('''# Revision C

Read **Revision_C.pdf**. The full edited LaTeX project is **main.tex**, **chapters/**, and **figs/**. Open this main.tex to compile the separate copy; existing ProofFinal and Overleaf were not edited.

**REPORT.md** gives the editorial decisions, chapter reductions, checks, and exact source evidence. **review/06_Completed_study_wording.md** records the completed-result drafting instruction while the flexor data are collected.

Use **review/changes.diff** and the section notes for selective integration. The SNS and SADb chapter rewrites are intentionally substantial; torque and other changes are more localized. New Appendix D and E require the corresponding two includes in main.tex and their existing figure files.

Canonical verification used MiKTeX pdflatex, BibTeX, and three further pdflatex passes, with a fresh output directory and BIBINPUTS/BSTINPUTS pointing to this folder. The final PDF and build log are included. Do not overwrite the actively reviewed Overleaf project wholesale.
''',encoding='utf-8')
print('Packaged',C/'Revision_C.pdf');print('Counts',allbefore,allafter);print('Integrity',integrity)
