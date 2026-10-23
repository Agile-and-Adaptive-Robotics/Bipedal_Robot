from pathlib import Path
import json
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';S=C/'review/overleaf_upload_20261023'
actions=json.loads((S/'review_actions.json').read_text(encoding='utf-8'));assert len(actions['resolutionLog'])==40 and len(actions['remaining'])==34 and all('23 October' in t for t in actions['remaining'])
text='''# Overleaf upload and advisor review — 2026-10-23

Uploaded Revision C into the existing Bolen-Dissertation project at Ben's explicit authorization. Replaced individual chapter files and main/support files, preserving folders and review threads. Updated the roadmap, Background figures, torque figures including final physical-test exports, SADb screenshots, SNS figures, Discussion schematics, and Appendix B circuit figures. Existing unchanged force-chapter and reference figures remain. No reports or review artifacts uploaded.

Chapter 4 now contains separate complete pre-optimization, optimized pinned-flexor, and optimized pinned-extensor GoF tables (4.1–4.3), plus a two-row flexor/extensor correction table (4.4), with four parameter columns chi0–chi3 as Ben confirmed. The earlier cross-configuration fit remains in Appendix C. No GoF results were recalculated or changed. The extensor-identified wrapping coefficient is explicitly distinguished from the pinned-flexor fit.

Final Overleaf compilation: Errors 0, Warnings 0; inherited typesetting warnings remain. Overleaf PDF has 259 physical pages; local MiKTeX verification PDF has 257. Both show Tables 4.2 and 4.3 on printed page 96. The Overleaf table page was visually inspected. No source text or comments were added to the protected isometric-force chapter.

Advisor review: 40 straightforward comments resolved; 34 substantive threads remain open, each with a new response describing the revision. This does not claim advisor acceptance. Accepted all 13 remaining grouped tracked-edit markers after replacement; Overview shows no pending tracked edits. The upload regrouped the original 51 advisor tracked edits and Ben's 3 tracked deletions, so grouped-marker counts differ from the original inventory. review_actions.json retains thread text, submitted replies, resolutions, and accepted groups.

The source-ZIP download event could not be verified. A labeled pre-upload history checkpoint was created instead. User authorization covered uploads and routine advisor replies/resolutions. No commit or push was performed.
'''
(S/'REPORT.md').write_text(text,encoding='utf-8')
p=C/'REPORT.md';s=p.read_text(encoding='utf-8');s=s.replace('# Revision C — final local sweep, 2026-10-23','# Revision C — uploaded and advisor review updated, 2026-10-23')
s=s.replace('These are local editorial dispositions, not Overleaf resolutions or advisor approval. No comment was resolved, no tracked change was accepted/rejected, and no source was uploaded. The historical resolved archive remains preserved.','Updated in Overleaf at Ben\'s authorization: 40 straightforward comments resolved, 34 substantive threads left open with replies, and all remaining tracked edits accepted. This does not indicate advisor approval. The historical resolved archive remains preserved. See review/overleaf_upload_20261023/REPORT.md and review_actions.json for the actual UI actions.')
s=s.replace('Overleaf file replacement may disturb source-linked review anchors; the complete current comments and local dispositions are retained outside the upload bundle. No upload or review-state change was performed by this pass.','Upload completed into the existing project using individual files rather than folder replacement. Comments and substantive replies were preserved. Final Overleaf compilation shows Errors 0 and Warnings 0; inherited typesetting warnings remain. The compiled Overleaf PDF has 259 pages. The local verification PDF has 257 pages.')
p.write_text(s,encoding='utf-8')
p=C/'COMMENT_DISPOSITIONS.md';s=p.read_text(encoding='utf-8');s='> Overleaf review updated 2026-10-23: 40 straightforward threads resolved; 34 substantive advisor threads remain open with replies; no pending tracked edits. The table below records editorial treatment; actual UI actions are in review/overleaf_upload_20261023/review_actions.json.\n\n'+s;p.write_text(s,encoding='utf-8')
p=R/'CHATGPT_REPORT.md'
with p.open('a',encoding='utf-8') as f:f.write('\n\n## Revision C uploaded; Chapter 4 tables restored — 2026-10-23\n\n'+text.split('\n',2)[2]+'\nEvidence: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/review/overleaf_upload_20261023/`.\n')
print('Recorded upload, table restoration, build status, and advisor-review actions')

