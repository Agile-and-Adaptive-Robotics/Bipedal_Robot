---
description: Dissertation-first status check and ranked task list for Bolen_Dissertation
---

# /diss

Run the bipedal-pm **dissertation-focus** skill end to end: status scan (chapter mtimes,
placeholder grep), figure integrity (\includegraphics existence/staleness), recent-work
mapping (reports/git -> INCLUDE / ALREADY IN / SKIP / LATER), then deliver the ranked
task list with the single highest-leverage item called out. Update
`Documentation\Reports and Papers\Dissertation\Notes\DISS_STATUS.md` if rows changed —
in-place status refreshes of existing rows only. Structural changes (adding, removing,
or renaming rows, e.g., after a chapter restructure) are proposed in the report for
Ben's confirmation, not applied.

Distinct command: `/diss-status` (dissertation plugin) is the read-only summary of the
recorded state and never edits; `/diss` is this active scan. If Ben types a short
`/diss`, confirm which one he means before editing anything.

$ARGUMENTS may narrow scope, e.g. "figures only", "CPG section", "defense prep".
