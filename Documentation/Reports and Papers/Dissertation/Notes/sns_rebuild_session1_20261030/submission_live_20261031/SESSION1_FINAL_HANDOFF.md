# Session 1: final live SNS revision

Machine: **easteregg2**. Ben explicitly authorized live Overleaf edits during the final submission pass. This supersedes the earlier editor-only manuscript handoff instruction for this pass.

## Delivery

- Live project: https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb
- Final PDF: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\sns_rebuild_session1_20261030\submission_live_20261031\Bolen_Dissertation_SESSION1.pdf`
- Complete downloaded live source archive: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Archive\session1_submission_20261031\Overleaf_source.zip`
- Exact edited chapter copies: `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Archive\session1_submission_20261031\source`. These are archival snapshots; the local working ProofFinal was not overwritten with an older source.
- Browser-download PDF: `C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation (2).pdf`.
- Final downloaded PDF: **310 pages**. Live compiler: **Errors 0, Warnings 0**. PDF text has no unresolved `??` references. No edits to the published muscle chapter or torque chapter.

## Requested changes completed

- Figure 6.1 now pairs the original walker (A) with Ben's supplied ground-test fall frame (B), in a compact float.
- The old shared-RG controller is standalone **Figure 6.2**, on its own sideways page (PDF page 154, printed page 133). It is no longer a continued panel.
- The newer supported treadmill still is separate Figure 6.4; 9 left / 8 right qualifying swings over 10 s, with 40.5% mean net body-weight assistance.
- Native W2L108 subnetwork schematic includes contact and bilateral joint flexor/extensor type-II feedback (Figure 6.6).
- Native public rat subnetwork figure includes the biarticular subsystem. Caption identifies the public precursor rather than claiming the final 2023 model was reproduced.
- Added actual Gait2392 geometry with muscle routes. Preserved the CAD/SNS-Simscape/BPA standing example and constructive future-work wording.
- Moved the matched-input engine plot and Gait2392 IK diagnostic plot to Appendix D while retaining labels and references. Existing editor D4 replacement was preserved.
- Added concise GUI/headless/API/SDK explanation and acceleration procedure: stored topology, sparse connections, preallocated arrays, compiled Numba calculations, and reference checks.
- Table 6.4 compares complete loops: Gait2392 8 s simulation, 122.120 to 13.253 s wall time (9.21x); W2L80 4 s simulation, 17.208 to 0.970 s wall time (17.75x). Different references are named. No unmatched AnimatLab/Simscape timing was invented.
- Retained compliant-tendon and damping prototype coverage; added the separate BPA pressure-dynamics and synthetic hysteresis test, with dynamic identification framed as future research.
- Additional old execution measurements were moved to Appendix D, preserving their numerical record rather than silently deleting them.

## Source line evidence

- `60-sns.tex:10` — fig:animatlab_model_architecture.
- `60-sns.tex:10` — fig:animatlab_shared_rg.
- `60-sns.tex:67` — Headless execution.
- `60-sns.tex:108` — fig:w2l108_subnetworks.
- `60-sns.tex:306` — fig:rat_native_subnetworks.
- `60-sns.tex:316` — tab:sns_runtime.
- `60-sns.tex:337` — synthetic hysteresis term.
- `60-sns.tex:335` — compliant-tendon prototype.
- `95-AppendixD.tex:51` — fig:w2l_crossengine.
- `95-AppendixD.tex:166` — fig:gait2392_ik_overlay.
- `95-AppendixD.tex:127` — Additional execution checks..

## Timing provenance

- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\backend_benchmarks.json`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\BENCHMARK_FINAL.json`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\EXPERIMENT_RECORD.md`

## Saved next modeling instruction

> Keep working on Gait2392. Fix stage 3 with the same controller, then work toward less weight support in stage 4. Don't let the knees hyperextend. Use the speed-ups and keep the separate training record updated.

No new controller training was run during this final manuscript pass. The separate experiment records remain authoritative for model performance.
