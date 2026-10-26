from pathlib import Path
import json,sys
root=Path(r'D:\Github\Bipedal_Robot')
w=root/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026'
c=json.loads((w/'followup_ethics_closing/change.json').read_text(encoding='utf-8'))
text=(w/'followup_ethics_closing/after_15-background.tex').read_text(encoding='utf-8')
line=text[:text.index(c['new'])].count('\n')+1
complete='--complete' in sys.argv
status='Exact whole-file comparison passed. Compilation pending.'
if complete:
    build=json.loads((w/'build_result.json').read_text(encoding='utf-8'))
    (w/'followup_ethics_closing/build_result.json').write_text(json.dumps(build,indent=2),encoding='utf-8')
    status='Exact whole-file comparison passed. Overleaf Errors 0 / Warnings 0. Full local build `build_astra_05` passed with no errors, undefined references, or missing figures; existing overflow messages remain.'
begin='<!-- ETHICS_CLOSING_20261026_START -->'; end='<!-- ETHICS_CLOSING_20261026_END -->'
entry=f'''{begin}
# Follow-up: ethical-case paragraph closing

Ben requested a replacement for “The failure modes of these experiments damage hardware, not tissue.” Applied a one-sentence tracked replacement in live `15-background.tex`, Section `sec:bg_ethics`, local snapshot line {line}:

> {c['new']}

The replacement closes the paragraph around repeatability and reduced reliance on invasive animal preparations, then leads into the next paragraph's discussion of how those preparations appear in circuit schematics. It makes no claim that animal use has already been quantitatively reduced. No numbers, citations, other prose, comments, or existing tracked-change dispositions were changed. {status}

Before/after source and exact replacement: `LIVE_ASTRA_20261026/followup_ethics_closing/`.
{end}

'''
p=root/'CHATGPT_REPORT.md'; old=p.read_text(encoding='utf-8-sig')
if begin in old:
    a=old.index(begin); b=old.index(end,a)+len(end)
    old=old[:a]+old[b:].lstrip('\r\n')
    if old.startswith('---\n\n'): old=old[5:]
p.write_text(entry+'---\n\n'+old,encoding='utf-8')
print('Updated CHATGPT_REPORT.md; source line',line)
