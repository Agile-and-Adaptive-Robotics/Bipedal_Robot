# Li–Deng-derived W2L108: bounded ground calibration

Machine: **easteregg2**, 4 November 2026. This is the newer 108-cell,
216-chemical-edge graded W2L model, not the biarticular rat model or BiW2Lvar+.
It retains source muscle Ia/II/Ib feedback and contact circuitry.

Rustuna completed 32 trials in 130.59 s using a persistent scientific worker
and the existing accelerated neural implementation. The selected 10 s trial
repeated all saved numeric arrays byte-exactly in a fresh process. It passes
the existing supported ground-stepping gate: 8 left and 7 right unloaded,
whole-foot-clear swings, and 100% alternating qualified transitions.
Minimum anatomical knee flexion was 1.946°/2.087°, without hyperextension.

Upward rail/guide assistance averages **29.519% BW after 1 s**. Continuous
harness force is zero; mechanical rails and 30 mm-slack guides remain.
Signed net assistance over the entire trial, including startup, is 35.069% BW.
These are different windows and metrics. Belt-relative progress is +0.908 m
in 10 s; the treadmill is moving, so this is not overground displacement.
Foot-force p95 is 1.109 BW and peak is 11.727 BW; impact transients still need
work. The selected fresh simulation ran at 2.94 simulated seconds per wall
second, excluding rendering; this is throughput, not a speedup ratio to SNS.

Shared bilateral fitted parameters:

| Parameter | Value |
|---|---:|
| PF→MN recruitment scale | 1.115222745 |
| Hip PF→MN additional scale | 1.261249111 |
| Neural time-constant scale | 1.144302445 |
| Toe→MN maximum conductance | 0.0488441541 µS |

The source home pose and explicit four-edge ankle PF phase correction were
retained. Muscle physics, damping, rail geometry, and feedback signs were not
fitted. Existing toe→MN wiring was tuned; no toe→RG pathway was added. Objective:
first require finite, nonhyperextending, alternating loaded/clear steps; then
minimize upward assistance plus p95/peak impact and low-progress penalties.
The inherited adapter's source-condition phrase `runtimeclock1` is stale for
this expanded search: the actual consumed clock scale is in the table and
MANIFEST/BEST/RESULT metadata. Parameters and actual trajectories are authoritative.

Full paths on easteregg2:

- Source: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer_Standalone.asim`.
- Source SHA256: `01ece6edcbcf17c274410b23f1ce7603eab288aadb946d6b548c5102ed5c1154`.
- Study: `D:\Users\Ben Bolen\NeuroCampaign_20261031\LiDeng_W2L_20261104\SUMMARY.json`.
- Configuration/ranges: same root, `LAUNCH.json`; trials: `TRIALS.jsonl`.
- Fresh proof: same root, `FRESH_VERIFICATION.json`.
- Exact model: same root, `fresh_repeat\W2L_VirtualWalker_Newer\repeat\model.xml`.
- Actual recording: same folder, `trajectory.npz`.
- GIF/still: same folder, `SUPPORTED_STEPPING_REPLAY.gif` and `SUPPORTED_STEPPING_STILL.png`.

Portable GIF/still are copied into
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104`
as LiDeng_W2L108_ground.gif/.png. Rendering uses the recorded 5–8.2 s segment,
not a newly integrated or corrected motion. This calibration is not a complete
Hunt/Bolen curriculum or an independent-balance result.
