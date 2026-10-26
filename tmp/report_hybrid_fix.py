from pathlib import Path
import json,shutil
root=Path(r'D:\Github\Bipedal_Robot')
w=root/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026'
c=json.loads((w/'followup_hybrid/change.json').read_text(encoding='utf-8'))
text=(w/'followup_hybrid/after_15-background.tex').read_text(encoding='utf-8')
line=text[:text.index(c['new'])].count('\n')+1
entry=f'''# Follow-up clarity correction: feed-forward and feedback control

Applied in live Overleaf, `15-background.tex`, Section `sec:bg_neuromech`, local snapshot line {line}. Ben identified the sentence beginning “Hybrid work makes the complementarity…” as unclear. Replaced that sentence with three explicit sentences:

> {c['new']}

The revision identifies the two signals, defines the step-down test, and states what the muscle-dependent finding means. It was verified against the [original paper]({c['source']}). The same citation key is retained; no numerical result was added or altered. Whole-file readback confirms this was the only content change in the chapter. Applied in Reviewing mode; no comments resolved or tracked changes accepted/rejected.

Verification: Overleaf Errors 0 / Warnings 0. Full local build `LIVE_ASTRA_20261026/revised_project/build_astra_03/` passed all four compilation steps with no errors, undefined references, or missing figures. Existing overflow messages remain. Before/after source, the exact replacement, and a screenshot are in `LIVE_ASTRA_20261026/followup_hybrid/`.

---

'''
p=root/'CHATGPT_REPORT.md'; p.write_bytes(entry.encode('utf-8')+p.read_bytes())
shutil.copy2(w/'build_result.json',w/'followup_hybrid/build_result.json')
print('CHATGPT_REPORT.md updated; live source line',line)
