from pathlib import Path
import json,datetime,sys
root=Path(r'D:\Github\Bipedal_Robot')
w=root/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026'
c=json.loads((w/'followup_cat_transition/changes.json').read_text(encoding='utf-8'))
text=(w/'followup_cat_transition/after_15-background.tex').read_text(encoding='utf-8')
line=text[:text.index(c[0]['new'])].count('\n')+1
done='--complete' in sys.argv
status='Overleaf recompile and local full-build verification pending.'
if done:
    build=json.loads((w/'build_result.json').read_text(encoding='utf-8'))
    (w/'followup_cat_transition/build_result.json').write_text(json.dumps(build,indent=2),encoding='utf-8')
    status='Overleaf recompile: Errors 0 / Warnings 0. Local full build `build_astra_04` completed with no errors, undefined references, or missing figures; existing overflow messages remain.' if not any(build['returncodes']) and not build['undefined'] and not build['missing'] else 'See `followup_cat_transition/build_result.json` for local build results. Live compilation is recorded separately.'
begin='<!-- CAT_TRANSITION_20261026_START -->'; end='<!-- CAT_TRANSITION_20261026_END -->'
entry=f'''{begin}
# Follow-up: transition from cat schematics to experimental procedures

Updated {datetime.datetime.now().astimezone().isoformat(timespec='seconds')}. Ben asked that the paragraph transition explain how apparently normal or happy-looking cat schematics differ from the experimental preparation. Applied in live Overleaf, `15-background.tex`, Section `sec:bg_ethics`, local snapshot line {line}.

Restructured the first paragraph around the schematic's appearance, the brainstem-transection preparation, the archival figure, and the specific meaning of replacement. It describes the impression a drawing can give without claiming to know the animal's emotional state. Retained the existing source keys for the preparation and footage; removed the unsupported “almost exclusively in schematic” generalization and the vague accounting/counterweight wording.

New first paragraph:

> {c[0]['new']}

New opening of the next paragraph:

> {c[1]['new']}

The rest of the procedural paragraph, the figure, its caption, and the citations remain as they were in the freshly captured live source. Two exact text selections were replaced in Reviewing mode. Whole-file readback confirms these were the only changes. No comments resolved or tracked changes accepted/rejected. {status}

Before/after source and exact replacements: `LIVE_ASTRA_20261026/followup_cat_transition/`.
{end}

'''
p=root/'CHATGPT_REPORT.md'; old=p.read_text(encoding='utf-8-sig')
if begin in old:
    start=old.index(begin); stop=old.index(end,start)+len(end)
    old=old[:start]+old[stop:].lstrip('\r\n')
    if old.startswith('---\n\n'): old=old[5:]
p.write_text(entry+'---\n\n'+old,encoding='utf-8')
print('CHATGPT_REPORT.md updated at source line',line)
