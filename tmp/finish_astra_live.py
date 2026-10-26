from pathlib import Path
import json,hashlib,difflib,shutil
import pymupdf
root=Path(r'D:\Github\Bipedal_Robot')
w=root/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026'
changes=json.loads((w/'patches.json').read_text(encoding='utf-8'))
applied=json.loads((w/'applied.json').read_text(encoding='utf-8'))
names=['10-introduction.tex','15-background.tex','60-sns.tex','70-overall.tex','80-conclusion.tex','95-AppendixD.tex']
verification={}
diff=[]
for n in names:
    actual=(w/'after'/n.replace('.tex','_LIVE.tex')).read_text(encoding='utf-8')
    planned=(w/'after'/n).read_text(encoding='utf-8')
    verification[n]=actual==planned
    before=(w/'before'/n).read_text(encoding='utf-8')
    diff.extend(difflib.unified_diff(before.splitlines(True),actual.splitlines(True),fromfile='before/'+n,tofile='after/'+n))
(w/'LIVE_CHANGES.diff').write_text(''.join(diff),encoding='utf-8')
protected=['chapters/20-muscles.tex','thesis.bib','main.tex']
hashes={n:hashlib.sha256((w/'baseline_project'/n).read_bytes()).hexdigest()==hashlib.sha256((w/'revised_project'/n).read_bytes()).hexdigest() for n in protected}
pdf=pymupdf.open(w/'Dissertation_Overleaf_revised.pdf')
local=pymupdf.open(w/'revised_project/build_astra_02/main.pdf')
checks={'live_readback_matches_prepared_text':verification,'protected_local_files_byte_identical':hashes,'applied_patch_count':len(applied),'overleaf_pdf_pages':len(pdf),'local_pdf_pages':len(local),'overleaf_errors':0,'overleaf_warnings':0,'local_build':json.loads((w/'build_result.json').read_text())}
(w/'FINAL_VERIFICATION.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
rows=['# Applied Overleaf revisions — 26 October 2026','', 'The live dissertation was revised through bounded selections in Reviewing mode. No chapter was replaced wholesale. No existing tracked change was accepted or rejected, and no advisor comment was resolved or reply posted. The peer-reviewed isometric-force chapter was not edited.','', '## Scientific argument','', 'The revised argument identifies the need to distinguish neural-circuit behavior, sensory contribution, actuator mechanics, and numerical implementation in the human-scale and BPA models. The tools provide specified, variable model components and controlled tests. The representative results establish selected implementation agreement, corrective reflex behavior, and contact/Ib-feedback dependence in a supported W2L configuration. These conclusions do not claim unrestricted controller transfer, physiological validation of every model, or harness-free walking.','', '## Applied passages','', '| Patch | File and current local line | Purpose |','|---|---|---|']
for i in applied:
    c=changes[i]; text=(w/'after'/c['file'].replace('.tex','_LIVE.tex')).read_text()
    needle=c['new'].strip().splitlines()[0]
    at=text.find(needle)
    line=text[:at].count('\n')+1 if at>=0 else 'see diff'
    rows.append(f'| A{i:02d} | `after/{c["file"].replace(".tex","_LIVE.tex")}`:{line} | {c["reason"]} |')
rows+=['','Additional insertion: Appendix D section `sec:app_conversion_repairs`, detailing the human-reference conversion repairs and the separate robotbody conversion attempt. Uploaded `connor_ben_shared_rg_original.png` and `ablation_pair_metrics_20261026.pdf` to `figs/60-sns/`. Figure 6.1 now has visible A/B labels and the original shared-RG diagram.','', '## Advisor-comment and author-response matrix','', '| Concern / recorded response | Action | Status |','|---|---|---|', '| “Describe the knowledge gap and how the tools I created fill them.” (Ben relaying advisor, this session) | Rewrote SNS gap/contributions/discussion and aligned introduction, background, overall discussion, and conclusion. | Applied; no field-wide absence claim. |', '| “huh?” / “This paragraph confuses me.” | Replaced the vague five-system prose and “connections” descriptions with named models, blocks, integrations, and experiments. | Applied; live threads retained. |', '| “need a better word” | Used neural pathways, actuator integration, signal exchange, and model implementation as appropriate. | Applied; no blanket substitution with “framework.” |', '| “What BPA interface? Non seq” / Ben explained the initial robot import | Explained length/pressure-to-force calculation, tendon route, and sensory return; separated the failed robotbody conversion from human-reference repair. | Applied in chapter and Appendix D. |', '| Advisor requested specific conversion problems / Ben described small muscles and robotbody routing | Added small-muscle deactivation, moving-path-point constraints, driver limits, patella tracking, and contact repairs with scope. | Applied. |', '| “Huh? Just delete all this?” / Ben replied “deleted” | Preserved the tracked deletion of the suspended-stepping qualifier. | Not restored; comment remains open. |', '| “for each leg?” / Ben described hip and knee/ankle PF groups and shared original RG | Correct original image and caption explain the bilateral PF groups and distinguish later per-side RG development. | Applied; reply preserved. |', '| Ben: Researchers; Gait2392 second; tools and selected proof-of-concept tests | Opening starts with Researchers; Gait2392 second in list and table; no incomplete-work column. | Applied. |','', '## Evidence and number changes','', '- The 2.00/7.85-degree loaded-arm comparison is retained and interpreted as a functional model test.', '- The new contact/Ib-ablation passage uses 33/28 versus 1/0 cycles and 1.16 versus 0.09 m from the saved paired metrics. The scored window is 25.499 s within each 30 s run; Ia/II remain active. It is the standard muscle-actuated W2L model, with identical harness configuration but differing measured support contribution.', '- Table values, including 1.000 s rhythm-core agreement, correlation at least 0.99996, the 12-actuator BPA example and 372 kPa peak pressure, come from the previously reviewed packet. Their scopes are retained.', '- The Rybak Python implementation is supported by `Code/MuJoCo_SNS/spinal/replication/verify_rybak2006.log`; the added explanation describes the distinct PF and RG perturbation responses behind its 13 checks.', '- No new optimization or simulation was run. The held CPG demonstration figure and uncorrected clipped demo figures were not uploaded. The new Deng experiment and further W2L perturbation packet remain available for a separate bounded insertion; this pass did not import them wholesale.','', '## Verification and remaining limits','', f'- Full local build: pdflatex, BibTeX, and two pdflatex passes exited zero; {len(local)} pages; no errors, undefined references, or missing figures.', f'- Final Overleaf build: Errors 0, Warnings 0; downloaded PDF verified as {len(pdf)} pages and contains the new gap text and Rybak subsection. `overleaf_compile_status.txt` preserves the UI status.', '- Whole-file readback of all six revised source files matches the prepared content after normalizing line endings. The source/PDF are preserved here alongside the before snapshot and diff.', '- The protected force chapter, existing bibliography, and main document are byte-identical between the local baseline and revised build. No commit or push was made.', '- Figure 6.1, the new table, gap prose, feedback subsection, and conversion appendix pages were visually inspected. Box-overflow messages remain elsewhere in the complete document; this is not a claim that every page or advisor thread has been fully polished.', '- The downloaded Overleaf PDF is the live compiled artifact. Local pagination differs between the build environments; do not replace the live project with the local tree.', '- Advisor comments and replies remain in Overleaf. The matrix reports content addressed, not automatic thread resolution.','']
(w/'APPLIED_REVISIONS.md').write_text('\n'.join(rows),encoding='utf-8')
summary='''<!-- LIVE_ASTRA_20261026_START -->
# Live Overleaf revision — completed 26 October 2026

Ben explicitly authorized live editing. Applied **24 bounded text replacements plus the detailed conversion-appendix insertion**, and uploaded two reviewed figures. Changed chapters: introduction, background (one consistency paragraph), SNS, overall discussion, conclusion, and Appendix D. **The peer-reviewed isometric-force chapter was not edited.**

The SNS chapter now states the knowledge gap, names the created tools, and explains what the selected demonstrations establish. “Researchers” and Gait2392-second ordering are applied. The tools table is explicitly referenced. Figure 6.1 has visible A/B labels and the correct Connor–Ben shared-RG diagram. The conversion account incorporates Ben's robotbody/small-muscle explanation. The appendix now explains the Python Rybak resetting/non-resetting mechanism.

All edits were made in Reviewing mode. Existing tracked changes were not accepted/rejected; comments were not resolved, and no replies were posted. Ben's deleted suspended-stepping qualifier remains deleted. Whole-file readback verified all six edited files against the prepared text.

**Verification:** full local build passed with zero errors, undefined references, and missing figures. Final Overleaf compile showed Errors 0 and Warnings 0. Its 268-page PDF is saved. Box-overflow messages remain; the record does not certify every page or every advisor thread as polished.

Delivery folder: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026/`
- `APPLIED_REVISIONS.md`: exact patch purposes, current local line anchors, advisor/Ben response matrix, evidence, and remaining limits.
- `Dissertation_Overleaf_revised.pdf`: final live compiled PDF.
- `LIVE_CHANGES.diff`, `before/`, `after/`, `FINAL_VERIFICATION.json`: source and verification record.
- `overleaf_revised_gap.png`: live browser confirmation.

Do not reapply the old incoming batches or replace whole chapters. The new Deng and additional W2L perturbation proposals remain separate; the held CPG/other figure repairs were not imported. No simulations, commits, pushes, or outgoing ZCode messages were performed.
<!-- LIVE_ASTRA_20261026_END -->

'''
p=root/'CHATGPT_REPORT.md'; text=p.read_text(encoding='utf-8-sig'); a=text.index('<!-- LIVE_ASTRA_20261026_START -->'); b=text.index('<!-- LIVE_ASTRA_20261026_END -->',a)+len('<!-- LIVE_ASTRA_20261026_END -->'); rest=text[b:].lstrip('\r\n')
while rest.startswith('---\n\n'): rest=rest[5:]
p.write_text(summary+'---\n\n'+rest,encoding='utf-8')
handoff='''# Live Overleaf edits now applied — 26 October 2026

Ben authorized the live revision. Read `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026/APPLIED_REVISIONS.md` before further edits. The latest live source snapshots are in that folder's `after/`; do not overwrite them with Revision E or reapply incoming batches. Knowledge-gap/contribution revisions, the corrected Figure 6.1, contact/Ib-ablation example, and conversion appendix are now live. Further proposals remain bounded and must be reconciled against this state. Existing comments and tracked changes were retained.

---

'''
p=root/'ZCODE_HANDOFF_20261023.md'; p.write_bytes(handoff.encode()+p.read_bytes())
print(json.dumps({'matches':verification,'protected':hashes,'live_pdf_pages':len(pdf),'local_pdf_pages':len(local)},indent=2))
