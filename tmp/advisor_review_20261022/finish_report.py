from pathlib import Path
import json,shutil,hashlib
H=Path(__file__).resolve().parent;R=H.parent.parent;N=R/'Notes/advisor_review_20261022';p=R/'CHATGPT_REPORT.md';s=p.read_text(encoding='utf-8')
s=s.replace('Updated from the live **Bolen-Dissertation** Overleaf Review → Overview on 22 October 2026.','Updated from the live **Bolen-Dissertation** Overleaf Review → Overview, final capture **22 October 2026, 12:22:09 Pacific**.')
s=s.replace('Line numbers and context refer to the fresh pre-upload source snapshot in tmp/advisor_review_20261022/overleaf_before/chapters, not the differing local ProofFinal chapters or PDF line numbers.','Line numbers and context refer to the fresh pre-upload source snapshot in tmp/advisor_review_20261022/overleaf_before/chapters, except the later live background snapshot saved as Notes/advisor_review_20261022/background_final.tex. They do not refer to the differing local ProofFinal chapters or PDF line numbers.')
s=s.replace('A subsequent layout check led to six further exact-fragment caption/sizing changes to the two combined figures: shorter captions and 0.94-textwidth panels retain figure text above 10 pt.','Subsequent layout and technical checks led to nine further exact-fragment changes to the two combined figures: shorter captions and 0.94-textwidth panels retain figure text above 10 pt, and “fix axial displacement” precisely identifies the transverse constraint.')
s=s.replace('The new background notes target the motor drawing,','The final background batch also identifies CPG-specific material misplaced in the broader computational-neuroscience section (BG-31 and BG-43–45). Move or condense it rather than repeating it in both sections. The new background notes target the motor drawing,')
s+='''
## Final synchronization checks

- **Live Overleaf:** 295-page current preview; Errors 0, Warnings 4. Both new combined-figure size warnings were eliminated. The four remaining warning cards concern pre-existing floats and a database table outside this figure-sync scope. The advisor's live source was not replaced with the different local draft.
- **Local verification:** ProofFinal/build_advisor_figures_20261022_final/main.pdf, 298 pages; four MiKTeX passes exit 0, zero undefined references/citations, and zero missing files. This local draft differs from the live review copy, so its page count is not a shortening result and must not be presented as the live Overleaf page count.
- Final live SNS and Discussion sources match the fresh pre-upload source plus only the intended figure patches (UTF-8 BOM ignored for comparison). The earlier copied introduction matches its intended single caption patch. No general prose rewrite was applied on Overleaf.
- The final rendered pages for Figures 6.6 and 7.2 pass independent visual and technical review. Minimum editable figure text is 10.152 pt. Their lowest text lies at 672.51 and 671.01 pt on a 792 pt page, respectively. No clipping remains.
- Chapter 3's 19 protected files remain byte-identical. Reports, prompts, native editing sources, and evidence stayed local. No advisor comments were resolved, no tracked edits were dispositioned, and no commits or pushes were made by this agent.
- GitHub was fetched again at completion; HEAD and origin still match 65b05aae. The latest torque-manuscript source commit remains fd1c13f4. No remote change was merged over the dirty working tree.
- Overleaf source/PDF export attempts later in the pass were blocked by Chrome. Final verification therefore uses the visible live PDF, compile UI, and exact source-copy comparisons, plus the independent local build. A completed earlier PDF download was inspectable, but it is not claimed as the final live PDF.

### Per-directive disposition

| Ben's request | Completed action / evidence |
|---|---|
| Review advisor contents, not the whole dissertation | Review Overview capture, exact anchors, pending edits, and relevant reply threads; no new full quality review. |
| Update local CHATGPT_REPORT | Current 73-thread/51-edit inventory above; previous report preserved locally. |
| Insert changed figures | Eleven PDF assets uploaded to the existing chapter figure folders; exact-fragment integration patches and Overleaf proof image saved in Notes/advisor_review_20261022/. |
| Prepare ZCode for reduction and new SNS structure | Replaced Notes/new_session_prompt_20261022.md; explicit AnimatLab → gap → tools/results → verification sequence, limited future-work handoff, and no mandatory Methods/Results/Discussion headings. |
| Avoid more bloated prose and captions | Prompt requires measured cuts, brief captions, comment-specific repairs, and independent review; both combined-figure captions shortened in this pass. |
| Check latest torque work | Fetch and commit/source inspection recorded above, including contradictory sweep status and newer test/sweep records. No unverified numerical transplant. |
| Keep active advisor review intact | No chapter-file uploads, wholesale replacements, comment resolution, or tracked-edit acceptance/rejection; future prose changes remain targeted local proposals. |
| Never save reports on Overleaf | No report/prompt uploaded; only figure PDFs and targeted figure-related source fragments changed there. |
| Protect isometric-force/MATLAB figures | All 19 protected Chapter 3 files retain their hashes. |

Start the new ZCode session using the entire prompt file. It is a handoff for work still to do, not a claim that SADb/SNS reduction or torque revision has already been completed.
'''
p.write_text(s,encoding='utf-8')
for n in ['constraint_patch.json','final_caption_patches.json','final60.tex','final70.tex']:
    shutil.copy2(H/n,N/n)
intro=(H/'live70.txt').read_text(encoding='utf-8-sig');expected=(H/'overleaf_expected/10-introduction.tex').read_text(encoding='utf-8-sig')
print('INTRO_SINGLE_PATCH_MATCH',intro==expected)
assert intro==expected
shutil.copy2(H/'live70.txt',N/'final10.tex')
assert s.count('### SNS-')==15 and s.count('### DB-')==26 and s.count('### INTRO-')==34 and s.count('### BG-')==49
print('REPORT_SECTION_COUNTS_VERIFIED')
