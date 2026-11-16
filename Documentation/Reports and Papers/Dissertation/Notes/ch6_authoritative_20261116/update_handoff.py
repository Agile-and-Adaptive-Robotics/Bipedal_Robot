"""Append/update only this session's report block; preserve all prior bytes."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
path=ROOT/'CHATGPT_REPORT.md'
start=b'<!-- BEGIN CH6 AUTHORITATIVE 20261116 -->'
end=b'<!-- END CH6 AUTHORITATIVE 20261116 -->'
result=HERE/'BUILD_RESULT.json'
status=json.loads(result.read_text()).get('status','pending') if result.exists() else 'pending'
body=r'''
## Authoritative Chapter 6 / Appendix D revision (2026-11-16, easteregg2)

Requested local revisions completed in:

- `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Final\chapters\60-sns.tex`
- `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Final\chapters\95-AppendixD.tex`

Review, numerical/source-change ledgers, generators, and exact commit paths:
`Documentation\Reports and Papers\Dissertation\Notes\ch6_authoritative_20261116\REVISION_REPORT.md`,
`NUMERICAL_CHANGE_LEDGER.md`, `SOURCE_CHANGE_LEDGER.json`, `MANUSCRIPT.diff`, and `COMMIT_FILES.txt`.
Active figures are in sibling Figures/60-sns and Figures/95-AppendixD; superseded assets and companions are preserved with original paths and hashes in Deprecated/20261116_authoritative and FIGURE_MOVES.json.

The original syn6 0/12 timing failure remains valid; later native/adaptive timing succeeds offline but live-body acceptance remains unsuccessful. Air commands/sensory conditions, support mass denominators, stimulus phases/windows, W2L right/left knees, CAD actuation/support, and distinct timing scopes are explicit. Spiking cup 10 s endpoint corrected to 1.78 deg at 0.625 kg; former 2.02 deg retained as maximum at 9.688 s. Fine-step knee discrepancy rounds to 5.20e-5 deg. See the ledgers for every numerical addition/removal and interpretation change.

Historical checks also corrected the shared W2L95/26-core and2ms best-lag scope, the common40N forward-pulse/duplicate-arm protocol, older speed8-30s, support16.75% using411N, right-side force comparison, composite IK score, spiking asymmetry, and a stale Di Russo native-runtime blocker. The later native literature-model baseline is separate from the rejected960-trial Gait2392 adaptation and from the developed supported controller-transfer demonstrations.

Model changes: OLD = NEW for all points, transforms, stiffness ordering, bounds, flags, Xi arguments, controller weights, actuator laws, and recording arrays. None changed. No scientific simulation, optimizer, training campaign, new dynamic replay, or timing benchmark was run. Saved arrays/configurations/scorers were independently read and compared. Only two Final sources changed; all 21 other source/support files, including Introduction, protected muscle chapter, and Chapter 7, are byte-preserved.

Independent evidence, complete prose reread, and standalone figure review completed; all active vector text meets 10 pt at source-derived manuscript scaling. Build and page outcome is recorded separately in the revision report and BUILD_RESULT.json. Current build status: STATUS_PLACEHOLDER.

Unresolved evidence: EB475WS4 SSH blocked before handshake with port-22 permission denied after name-resolution and direct cached-host-key-address attempts; no remote files inspected. Air launch body_probe hash/full consumed-command traces not recovered; F/E 0 versus 0.15 are retained adapter declarations, not independently verified original input samples. Captured physical evaluator matched. Historical CAD compiled binary/full 2 ms trajectory absent. Early hybrid-relay full raw waveform absent. These are not promoted to fresh verification.

Initial working tree clean; pre-existing changes NONE. Exact session paths/statuses are in COMMIT_FILES.txt/JSON. No commit, push, Overleaf upload, ZIP, or unrelated manuscript edit. Generated builds and scratch are ignored.
'''.replace('STATUS_PLACEHOLDER',status).strip()
block=start+b'\n'+body.encode('utf-8')+b'\n'+end
old=path.read_bytes()
if start in old:
    i=old.index(start);j=old.index(end,i)+len(end)
    new=old[:i]+block+old[j:]
else:new=old+b'\n\n'+block+b'\n'
path.write_bytes(new)
print('Updated only the current-session CHATGPT_REPORT block; build='+status)
