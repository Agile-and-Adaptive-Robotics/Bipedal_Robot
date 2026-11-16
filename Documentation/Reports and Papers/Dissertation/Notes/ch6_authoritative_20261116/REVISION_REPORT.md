# Authoritative Chapter 6 and Appendix D revision

16 November 2026, easteregg2. Canonical sources are `Final/main.tex`, `Final/chapters/60-sns.tex`, and `Final/chapters/95-AppendixD.tex` under `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation`.

## Access and scope

`hostname` was the first command and returned `easteregg2`. The working directory and selected workspace were `D:\GitHub\Bipedal_Robot`. A unique text probe under repository `tmp` was created, read back exactly, and deleted before source research. Actual repository writes succeeded. Initial `git status --short` was empty.

Read the repository instructions, handoff, machine-specific dissertation instructions, both complete manuscript sources, prior laptop chapter/appendix/graphics/independent-review/synchronization/figure/build records, current dissertation status, newest available digest, and Session 1 index. Older campaign files named by standing instructions were not present in this checkout; their historical ProofFinal/ZIP workflow was not substituted for the explicit current Final instructions. The interrupted external staging directory was not an accepted revision or build source.

Only the two requested canonical TeX files changed. `CONTENT_CONSERVATION.json` verifies all other 21 Final source/support files byte-identical, including main.tex, the latest Introduction, protected published muscle chapter, Chapter 7, both bibliographies, class, and bibliography style. No source tree was replaced. Review and provenance files remain under Notes; active figures remain in the appropriate sibling chapter folders.

## Substantive result

The laptop's tools-focused chapter structure was useful and retained. Independent evidence review required corrections to experimental scope, chronology, force normalization, side identification, and timing definitions, rather than another general rewrite.

- The original syn6 0/12 PF timing failure remains. Later relay and native timing experiments are now included: fixed native timing works at nominal replay cadence; adaptive timing passes the specified replay cadences; subsequent live-body tests still fail their retained timing or physical acceptance. Offline computation, transmitted neural behavior, and coupled acceptance are distinguished.
- Air-stage sensory and sustained command differences are explicit, so unchanged learned weights do not imply a feedback-only comparison. Cycle, assistance, display, and response windows are separate.
- The adopted 99.860% air support statistic retains its legacy numerator and denominator. The fixed-ground mass explains the discrepancy with moving-model metadata. No adopted result was silently renormalized.
- Right/left W2L108 knees, pose-offset conventions, upward versus signed assistance, and scoring windows are explicit.
- CAD remains a separate supported hybrid-control demonstration with ideal reference tensions. Source hashes and saved replay arrays support the result, while the absent historical compiled binary and absent full 2 ms trajectories remain disclosed. New quantitative traces show loading, clearance, knee bounds, and contact penetration.
- Timing comparisons identify controller family, workload, baseline, repeats, and timer scope. The compact diagnostic's saved evaluator wall time is no longer called complete-process time.
- The spiking cup comparison now uses the same 10 s, 0.625 kg endpoint as the graded case: 1.78 versus 2.00 degrees, compared with 7.85 without feedback. The former 2.02 value is preserved as a different maximum statistic. The knee fine-step discrepancy is correctly rounded to 5.20e-5 degrees.
- The distinction between zero supplementary hand drive and retained shared descending drive remains intact. No harness-free walking or continuous curriculum across unrelated configurations is implied.
- Historical source checks further corrected the W2L representative-core/best-lag agreement claim, disclosed the default forward pulse and duplicate push arm, bound the older speed to its actual fitting window, corrected legacy support rounding to 16.75%, and replaced a stale Di Russo runtime blocker with the separately verified native baseline. Spiking asymmetry, force-table side/reference interval, and composite-score terminology are explicit.

`NUMERICAL_CHANGE_LEDGER.md` gives each addition, replacement, removal, and interpretation change with exact sources and scoring windows. `SOURCE_CHANGE_LEDGER.json` preserves every changed source passage, all old/new numerical tokens, and token changes. `MANUSCRIPT.diff` supplies the complete readable source diff. None of these records substitutes for the source-specific raw-data audits.

The final air review identified an additional material provenance limit: F/E commands 0 versus 0.15 are the retained adapter's declared settings, with 0.15 also in the fitted manifest. The exact original launch adapter and complete consumed-command traces were not recovered. Both manuscript files now say so explicitly. The declared stage definitions and unchanged learned weights therefore cannot establish a controlled feedback-only effect or sample-by-sample command identity. This qualification does not alter the independently recalculated motion, support, and stimulus measurements.

## Independent evidence and prose review

| Review | Record and result |
|---|---|
| Air, noise, stimulus, W2L108 | `AIR_EVIDENCE.md`, `AIR_RESULTS.json`, `AIR_CYCLE_RESULTS.json`; saved-array measurements and archived cycle detector independently recalculated. All 32 W2L108 numeric arrays match the saved independent replay. |
| Original and later PF | `PF_EVIDENCE.md`, `PF_ORIGINAL_RESULTS.json`, `PF_AUDIT_RESULTS.json`; original and later native timing rescored, source/model bindings checked, offline/live scopes retained. |
| CAD | `CAD_EVIDENCE.md`, `CAD_AUDIT_RESULTS.json`; raw counts/support, all 26 arrays in both saved replay pairs, model/source hashes, local cases, and geometry readback checked. |
| Timing, cup, knee | `TIMING_CUP_EVIDENCE.md`, `TIMING_RESULTS.json`; actual MAT samples decoded, exact timing source recovered read-only, 45- and 62-array comparisons and model identities verified. |
| Independent complete prose reread | `PROSE_FINAL_REVIEW.md`; narrow findings corrected and read back. |
| Final evidence/figure checks | `AIR_FINAL_REVIEW.md`, `PF_FINAL_REVIEW.md`, `CAD_FINAL_REVIEW.md`; material source/typography findings resolved. Compiled-page verification is separate below. |
| Historical implementation and spiking | `HISTORICAL_IMPLEMENTATION_REVIEW.md`, native and spiking reviews/checks/results; raw bridge sensor arrays, shared W2L signals, original Di Russo baseline/repeat, and spiking comparison metrics checked. Secondary Deng/export results checked against retained primary logs and scorers, with absent historical arrays/identities disclosed. |
| Historical W2L, IK, and mechanics | `HISTORICAL_W2L_REVIEW.md`, `HISTORICAL_IK_REVIEW.md`, `HISTORICAL_MECHANICAL_REVIEW.md` and reproducible calculations; windows, shared force protocol, phase/displacement/speed, support denominator, repeats, gates, IK values, force CSV, reference peaks, and repair code checked. |

No training, optimization, scientific simulation, fresh dynamic replay, or performance benchmark was run. Saved independent repeats were compared; the review does not claim they were newly simulated. A read-only MATLAB MAT-export attempt failed before script execution, first on preferences access and then licensing 63/-63.2. Base Anaconda SciPy successfully read the actual saved MAT contents instead. The initial myo SciPy loader error was resolved for needed calculations with the available base environment.

## Figures and conservation

| Active chapter | Delivered files / purpose |
|---|---|
| `Figures\60-sns` | `simulation_workflow_authoritative.pdf`, `sns_knee_reflex_authoritative.pdf`, `sns_cup_reflex_authoritative.pdf`, `biw2lvar_air_aligned_authoritative.pdf`; readable workflow and compact demonstrations, direct aligned air arrays. PNG/alt companions remain with active assets. |
| `Figures\95-AppendixD` | `AIR_stimulus_neural.pdf`, `AIR_stimulus_mechanics.pdf`, `PF_timing_authoritative.pdf`, `CAD_supported_traces_20261116.pdf`; direct raw-array stimulus/PF/CAD evidence. PNG/alt and available SVG companions remain with active assets. |

All 12 active vector figures meet the 10 pt floor under source-derived manuscript scaling. New figures use exact fixed 14 cm inclusion widths and reserve caption space. `CAD_FIGURE_QC.json` and `BUILD_FIGURE_FONT_CHECK.json` record the measurements. Standalone previews were visually inspected. Native signs, units, sample clocks, and every plotted sample are retained. Declared calculations include unit conversions, model transmission readouts, signal differences, and region-force sums. The PF intended primitive is evaluated by its original retained mapping on the calibration clock, rather than shifted to improve agreement. The inherited knee summary itself is an earlier 1 ms interpolation grid; this session does not relabel it as original variable-step data.

Four unreadable or redundant inherited inclusions were removed: normalized force-sequence grid, duplicate appendix knee figure, multi-perturbation grid, and low-resolution original shared-RG bitmap. The force table, knee protocol and main response, intervention results, original architecture explanation, source recordings, and original assets remain preserved. The shared-RG architecture is still distinguished from later bilateral RGs.

Every superseded active asset and companion is preserved in its chapter's `Deprecated\20261116_authoritative` directory, with original absolute path and SHA256 in `FIGURE_MOVES.json`. Intermediate session layout drafts are marked separately. No source recording or model asset was deleted.

## EB475WS4 follow-up requested during this session

Read-only SSH was attempted using the configured alias and fully qualified host name. Hostname resolution failed; explicit DNS lookup returned access denied. A cached address was then verified against the existing eb475ws4 host key and used for a direct connection. SSH failed before handshake with `port 22: Permission denied`. The DNS client cache query also reported unavailable CIM access. No remote file inventory or evidence was obtained, so this report makes no claim about locally saved EB475WS4 files. Machine-local connection details and keys were not copied into the public repository.

## Build and rendered-page verification

The installed normal MiKTeX startup failed on registry access. A repository-local runtime was recovered using existing readable installed assets and scratch copies of format/database configuration; installed/user originals were preserved. A minimal PDF probe succeeded. The full canonical build was then launched from Final with absolute BIBINPUTS/BSTINPUTS and outputs in a fresh directory outside Final.

Full-build outcome, exact commands, elapsed time, errors/warnings, page ranges, rendered-page review, and any remaining gap will be finalized here after the owned build finishes. Static checks and standalone typography are not described as a successful dissertation build.

## Model-change record

OLD = NEW for every model point, transform, stiffness ordering, bound, flag, controller weight, actuator law, and Xi argument. None changed. No result recording, SLX, MATLAB model, or controller source was modified. Figure scripts perform read-only extraction or geometry inspection and write their own derived artifacts.

## Delivery and commit scope

`COMMIT_FILES.txt` and `COMMIT_FILES.json` enumerate exact absolute paths, statuses, and any deletion/move, separating this session from pre-existing changes. The initial working tree was clean, so pre-existing changes are NONE. The two canonical TeX files, active and preserved figure changes, this Notes folder, and the appended root `CHATGPT_REPORT.md` are the intended session changes. Generated full-build output and scratch remain ignored and are not commit deliverables. No commit, push, ZIP, or Overleaf upload was performed.

Remaining evidence limits: remote EB475WS4 inventory blocked; original air launch-time `body_probe.py` hash and complete command telemetry not recovered (captured physical evaluator matched); early hybrid-relay report lacks full raw waveform; historical CAD compiled binary and full 2 ms trace not retained. Some secondary historical compiled identities, paired Simulink verification arrays, and MuJoCo matched-force samples are unavailable; those specific claims remain checked to saved summaries/primary scripts rather than fully rescored trajectories. These limits do not turn retained source-bound or saved-array checks into fresh dynamic verification.
