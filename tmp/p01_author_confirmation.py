from pathlib import Path
import json,shutil
R=Path(r'D:\Github\Bipedal_Robot');E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E';V=E/'review/targeted_reconciliation';W=V/'tools_inventory'
p=V/'P01_replacement.tex';s=p.read_text(encoding='utf-8');shutil.copy2(p,W/'P01_before_author_confirmation.tex')
anchor=r'    \item Integration of SNS-Simscape controllers'
i=s.index(anchor)
s=s[:i]+r'    \item A SolidWorks-derived robot assembly in Simscape Multibody coupled to a neural controller and BPA actuators using the force model from Chapter~\ref{ch:muscles}.'+'\n'+s[i:]
s=s.replace('Table~\\ref{tab:sns_capability} summarizes representative functional-verification tests and proof-of-concept demonstrations of these tools.', 'Table~\\ref{tab:sns_capability} presents selected functional-verification tests and proof-of-concept demonstrations illustrating the capabilities of this toolset.')
p.write_text(s,encoding='utf-8')
p=V/'P01_companion_table.tex';s=p.read_text(encoding='utf-8');s=s.replace('Each result refers to the named circuit or mechanical model.', 'The examples illustrate the toolset rather than exhaustively evaluating every implementation; each result refers to the named circuit or mechanical model.');p.write_text(s,encoding='utf-8')
note='''# Author confirmation — coupled CAD assembly and representative demonstrations

Ben confirms that the SolidWorks-derived robot assembly was coupled to a neural controller and BPAs using his force model. Its saved model is on the temporarily unavailable laptop. This is accepted as author-supplied project information and included in P01's contribution list. The specific saved artifact was not independently inspected in this review; no filename, numerical result, or performance measurement is invented. Laptop availability is editorial provenance and does not belong in dissertation prose.

Rybak's implementation is Python; no Simulink Rybak implementation is claimed. Representative proof-of-concept demonstrations illustrate the toolset. P01 does not require an exhaustive result for every tool and does not make additional ZCode experiments a condition for including an author-confirmed tool. New results can be reviewed separately for committee questions. The paired table remains a selection of demonstrations and need not have a separate numerical result for the CAD assembly.

P01_replacement.tex now explicitly lists the coupled SolidWorks-derived assembly. No Overleaf or active dissertation chapter was edited by this update.
'''
(W/'AUTHOR_CONFIRMATION.md').write_text(note,encoding='utf-8')
p=W/'TOOLS_AND_DEMONSTRATIONS.md';s=p.read_text(encoding='utf-8')
s=s.replace('Claim a complete imported robot with controller and BPA force only when its specific runnable saved model/result is identified. This is an evidence lookup, not a demand for new experiments.', 'Ben confirms that the imported robot assembly is coupled to a neural controller and BPAs using his force model. Its saved model is on the temporarily unavailable laptop. Include this author-confirmed contribution without inventing a file location or numerical result.')
s='''# Latest author clarification

The coupled SolidWorks-derived robot assembly is confirmed by Ben and is now explicitly included in P01. Its saved artifact is temporarily unavailable on the laptop. The table gives representative examples for the toolset, not exhaustive validation of every tool. Rybak remains listed as Python. See AUTHOR_CONFIRMATION.md. Earlier evidence distinctions below remain useful for attributing numerical results.

---

'''+s
p.write_text(s,encoding='utf-8')
p=V/'OVERLEAF_TARGETED_PATCHES.md';s=p.read_text(encoding='utf-8');p.write_text('# P01 author confirmation\n\nP01 now includes Ben’s confirmed SolidWorks-derived robot assembly with neural control and BPAs using his force model. Its saved model is on the unavailable laptop; this does not block the tool inventory. The table supplies selected proof-of-concept examples, not exhaustive results for each tool. Rybak is Python. See tools_inventory/AUTHOR_CONFIRMATION.md. No Overleaf edits were made.\n\n---\n\n'+s,encoding='utf-8')
p=V/'patches.json';data=json.loads(p.read_text(encoding='utf-8'))
for row in data:
 if row['id']=='P01':row['new']=(V/'P01_replacement.tex').read_text(encoding='utf-8');row['status']='Author-confirmed coupled CAD assembly included. Representative demonstrations illustrate the toolset; exhaustive per-tool results are not required.'
p.write_text(json.dumps(data,indent=2),encoding='utf-8')
for p in [R/'CHATGPT_REPORT.md',R/'ZCODE_HANDOFF_20261023.md']:
 s=p.read_text(encoding='utf-8');p.write_text('# P01 author confirmation — coupled CAD assembly\n\nBen confirms the SolidWorks-derived robot assembly with a neural controller and BPAs using his force model; its saved model is on the temporarily unavailable laptop. Include it as an author-confirmed contribution. P01 uses representative proof-of-concept examples, not exhaustive results for every tool. Do not block the prose on further experiments or ask Ben again for this artifact now. Rybak is Python. Details: Revision_E/review/targeted_reconciliation/tools_inventory/AUTHOR_CONFIRMATION.md.\n\n---\n\n'+s,encoding='utf-8')
print('Updated P01 and editorial provenance; no measurements invented or Overleaf changes made.')
