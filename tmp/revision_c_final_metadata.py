from pathlib import Path
import re,json
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
texts={p.name:p.read_text(encoding='utf-8-sig') for p in (C/'chapters').glob('*.tex')};ledger=[];seen=set()
for before in sorted((C/'review/before').glob('*.tex')):
    old=before.read_text(encoding='utf-8-sig')
    for m in re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}',old,re.S):
        for label in re.findall(r'\\label\{([^}]+)\}',m.group()):
            locations=[name for name,text in texts.items() if r'\label{'+label+'}' in text];seen.add(label)
            ledger.append({'label':label,'before':before.name,'after':locations,'status':'retained_or_moved' if locations else 'omitted; original asset/source preserved'})
for name,text in texts.items():
    for m in re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}',text,re.S):
        for label in re.findall(r'\\label\{([^}]+)\}',m.group()):
            if label not in seen:
                ledger.append({'label':label,'after':[name],'status':'added in Revision C or supporting appendix'});seen.add(label)
(C/'review/figure_dispositions.json').write_text(json.dumps(ledger,indent=2))
for name,note in [
('03_SADb.md','Final-sweep update (2026-10-23): the actual app now embeds all 1,393 records, including non-paper records. The old paper-only distribution limitation below is historical. The bubble screenshot is in the main chapter, and the table screenshot is in Appendix E. See final_sweep_20261023/app_status_correction.md.'),
('06_Completed_study_wording.md','Final-sweep update (2026-10-23): the physical tests are completed and their committed analysis-script plots are now inserted. The drafting assumption described below is historical. Final claims distinguish measured capacity from pointwise model agreement.'),
('07_Bilateral_and_torque_source.md','Final-sweep update (2026-10-23): finished physical-test plots are included; the knowledge-gap statement now enumerates all five combinations in Ben\'s advisor reply. The monthly scraper and new bipulley run remain active. Earlier status below is historical.')]:
    p=C/'review'/name;s=p.read_text(encoding='utf-8')
    if not s.startswith('Final-sweep update'):p.write_text(note+'\n\n'+s,encoding='utf-8')
print('Updated figure moves and marked superseded review notes as historical')
