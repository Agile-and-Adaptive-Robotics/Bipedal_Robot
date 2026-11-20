"""Complete interrupted delivery records without altering manuscript/build bytes."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
DISS=HERE.parent.parent
ROOT=DISS.parents[2]
resume=json.loads((HERE/'RESUME_VERIFICATION.json').read_text())
build=json.loads((HERE/'BUILD_RESULT_r7.json').read_text())
pdf=Path(build['pdf'])
assert hashlib.sha256(pdf.read_bytes()).hexdigest()==build['pdf_sha256']
assert build['status']=='compiled' and not any(build['gates'].values())
report=HERE/'REVISION_REPORT.md'
text=report.read_text(encoding='utf-8')
first=text.index('## Build and rendered-page verification')
last=text.index('## Model-change record',first)
section=r'''## Build and rendered-page verification

The complete canonical `Final/main.tex` build succeeded in the fresh `tmp\ch6_authoritative_20261116\build_final_r7` directory, outside Final. It produces **318 pages**:

`D:\GitHub\Bipedal_Robot\tmp\ch6_authoritative_20261116\build_final_r7\main.pdf`

SHA256: `03314775553c05483786f42950e991512e18e8958169d6cc6c29a920a97ce104`.

The accepted run used the installed MiKTeX executables by absolute path, absolute BIBINPUTS/BSTINPUTS, BibTeX, and four LaTeX passes. Command/verification elapsed time was **34.984 s**, excluding preparation of the private runtime and human page review. `BUILD_RESULT_r7.json` preserves the exact commands, individual times/exit codes, runtime configuration, package/font provenance, frozen source hashes, and final gates. The installed normal runtime was blocked by registry access; the accepted private runtime uses hash-verified copies of existing support/fonts and scratch-only format timestamps. No installed/user-original package, font, format, registry setting, or manuscript support file changed.

All final automated gates pass: no build errors, undefined references/citations, duplicate TeX labels, missing graphics/files, oversized floats, or rerun requests. Nonfatal notices remain explicitly recorded in `BUILD_REVIEW.md`: the unchanged Appendix C 0.9562 pt overfull box, underfull boxes, PDF inclusion/version/page-group notices, and the duplicate PDF destination `page.1` (distinct from a duplicate TeX label). No zero-warning claim is made.

All **23 current Final files** match the accepted r7 source snapshot, rechecked after the connection returned. The saved PDF hash is unchanged. No source changes or another build were needed on continuation. The original-revision conservation check still finds only the two authorized TeX files changed and the other **21 files** byte-identical.

Publication-page renders are 144 dpi. Final Chapter 6 occupies **PDF pages 153-183**, printed pages 133-163; Appendix D occupies **PDF pages 273-314**, printed pages 253-294. Complete coverage combines individually reviewed pages, exact body-pixel matches to already reviewed r6 pages, and the resumed inspection of all eleven outstanding changed Chapter 6 pages. The air and CAD reviews close all final Appendix D pages. The current summaries are `BUILD_REVIEW.md`, `COMPLETION_PROSE_REVIEW.md`, `AIR_COMPILED_REVIEW.md`, and `CAD_COMPILED_REVIEW.md`. **No material figure, caption, table, sign, clipping, or page-transition defect remains in the revised chapters.** Some normal float whitespace and paragraph continuation remain documented as minor typography, not concealed as evidence problems.

The actual compiled vector-figure span scan reports a minimum **10.0006 bp**, with zero below-floor spans. Retained raster figures were checked visually. The original walker illustration uses a disclosed TeX viewport to omit tiny internal labels while preserving the source PNG bytes. Float barriers prevent earlier W2L/stimulus figures from interrupting later PF interpretation. Chapter 6 saves/restores its incoming widow/club penalties so other chapters do not inherit the local layout settings. The CAD caption was shortened after an actual 18.81631 pt oversized-float finding, preserving every numerical threshold and the 14 cm plot width.

The execution connection failed after the accepted build and before final documentation was written. On continuation, hostname/workspace/write-read-delete access, the clean working tree, and the absence of relevant running processes were checked first. Saved artifact hashes were verified, the remaining page checks were completed, and these records were updated. SSH was not retried, at Ben's direction. `RESUME_VERIFICATION.json` records that continuation and demonstrates that completed work was reused rather than rerun.

'''
section=section.replace(chr(92)*2,chr(92))
text=text[:first]+section+text[last:]
first=text.index('## Delivery and commit scope')
text=text[:first]+r'''## Delivery and commit scope

`DELIVERY_PATHS.txt` lists the full Windows paths to both canonical TeX files, all eight new/revised vector figures, the accepted PDF, reports, ledger, provenance, and commit records. The PDF remains an ordinary ignored local build file; no ZIP or download-only substitute was created.

The original revision began with a clean working tree. After the connection returned, the tree was again clean and **HEAD `0f493de8c3bcc0986e4dbf3a83ccaf3eba406e10` already contained the earlier manuscript/figure revision**. No agent made that commit. Completion changes are reports, conservation/ledger refreshes, and final review records; canonical manuscript/figure bytes remain unchanged on continuation. `COMMIT_FILES.txt` and `.json` now distinguish exact remaining completion paths to commit from earlier revision paths already committed. Pre-existing uncommitted changes at resume: **NONE**. No agent commit, push, Overleaf upload, scientific rerun, or unrelated source edit was performed.

Remaining evidence limits are disclosed in both the relevant prose and audits: original air launch adapter/complete consumed-command history unrecovered; early hybrid-relay full waveform absent; historical CAD compiled binary/full 2 ms trajectory absent; some secondary historical compiled identities, paired Simulink arrays, and MuJoCo matched-force samples unavailable. Those secondary claims remain source-bound saved-summary/script checks, not falsely described as full raw rescoring or fresh dynamic replay. The earlier EB475WS4 inventory was unavailable; Ben waived further SSH work on continuation. There is **no remaining local build or required page-review blocker**.
'''
report.write_text(text,encoding='utf-8')

ledger=HERE/'NUMERICAL_CHANGE_LEDGER.md'
value=ledger.read_text()
value=value.replace('local paragraph widow penalties affect layout only.',
    'chapter-scoped widow/club penalties affect layout only and restore their incoming values before the next chapter.')
ledger.write_text(value,encoding='utf-8')
layout_path=HERE/'LAYOUT_CORRECTIONS.json'
layout=json.loads(layout_path.read_text())
layout['main_pf_float_barrier']={'location':'60-sns immediately before Intermediate verification exposed',
    'reason':'Keep the preceding historical W2L figure with its experiment rather than inside later PF interpretation.'}
layout['widow_club_control']={'scope':'Chapter 6 only','value':10000,
    'incoming_values_saved':True,'incoming_values_restored_before_next_chapter':True,
    'scientific_parameters_changed':False}
layout['closure']='Verified in accepted r7; no remaining material visual finding.'
layout_path.write_text(json.dumps(layout,indent=2)+'\n')

status_path=DISS/'Notes/DISS_STATUS.md'
old=status_path.read_bytes()
start=b'<!-- BEGIN AUTHORITATIVE CH6 20261116 -->'
end=b'<!-- END AUTHORITATIVE CH6 20261116 -->'
status=r'''## Authoritative Chapter 6 and Appendix D completion - 16 November 2026

Canonical `Final/main.tex` compiles to **318 pages** using BibTeX and four MiKTeX passes, with no build errors, unresolved references/citations, duplicate TeX labels, missing graphics, oversized floats, or remaining rerun request. Final target-page review is complete; compiled vector text is at least 10 pt. The accepted PDF is `D:\GitHub\Bipedal_Robot\tmp\ch6_authoritative_20261116\build_final_r7\main.pdf`. Exact hashes, source conservation, nonfatal notices, and review coverage are in `Notes/ch6_authoritative_20261116/BUILD_RESULT.json`, `BUILD_REVIEW.md`, and `REVISION_REPORT.md`.

The revision preserves the tools framing while correcting PF chronology, air-command provenance/weight denominators, stimulus windows, controller/body/actuation identities, W2L sides/support/push protocols, timing scopes, native literature-model baseline chronology, and the matched spiking-cup endpoint. All numerical changes and sources/windows are in `NUMERICAL_CHANGE_LEDGER.md`; source diff/conservation and figure moves/hashes accompany it. No model parameter, controller weight, scientific recording, bibliography, or unrelated Final file changed. The other 21 Final source/support files remain byte-identical, including Introduction, protected muscle chapter, and Chapter 7's supplementary/shared-drive distinction.

Earlier sources/figures already reside in commit `0f493de8`; continuation completes reports and review records only. The connection interruption is resolved. The exact remaining commit list is `Notes/ch6_authoritative_20261116/COMMIT_FILES.txt`; no agent commit/push/upload/ZIP or simulation/optimization run. Further SSH was waived by Ben. Earlier dated checkpoints below remain historical.
'''.strip().replace(chr(92)*2,chr(92)).encode('utf-8')
block=start+b'\n'+status+b'\n'+end+b'\n\n'
if start in old:
    i=old.index(start);j=old.index(end,i)+len(end)
    new=old[:i]+block.rstrip()+old[j:]
else:new=block+old
status_path.write_bytes(new)

figs=[('60-sns',s) for s in ['simulation_workflow_authoritative','sns_knee_reflex_authoritative',
    'sns_cup_reflex_authoritative','biw2lvar_air_aligned_authoritative']]
figs += [('95-AppendixD',s) for s in ['AIR_stimulus_neural','AIR_stimulus_mechanics',
    'PF_timing_authoritative','CAD_supported_traces']]
paths=[DISS/'Final/chapters/60-sns.tex',DISS/'Final/chapters/95-AppendixD.tex',pdf]
paths += [DISS/'Figures'/chapter/(stem+'.pdf') for chapter,stem in figs]
paths += [HERE/n for n in ['REVISION_REPORT.md','NUMERICAL_CHANGE_LEDGER.md','BUILD_RESULT.json',
    'BUILD_REVIEW.md','FIGURE_MOVES.json','COMMIT_FILES.txt','COMMIT_FILES.json']]
paths += [ROOT/'CHATGPT_REPORT.md']
(HERE/'DELIVERY_PATHS.txt').write_text('Canonical sources, accepted local PDF, active vector figures, and review records:\n\n'+'\n'.join(str(p) for p in paths)+'\n')
print('Completed revision/status/layout/delivery records; manuscript and PDF bytes untouched.')
