# easteregg2 prompt — feedback ablation + Deng/Li original-experiment reproduction (Ben, 2026-10-26)

Paste-ready. Two orders, Ben-approved. Pre-register predictions BEFORE
any run (the 2026-10-26 post-hoc correction applies): write expected
outcomes to the campaign folder first, then run.

## Order A — the feedback-ablation pair (MuJoCo, this machine)

The chapter needs ONE saved fixed-tuning comparison: feedback intact vs
fully ablated on the walking BPA-routed model.
1. Freeze the `reports_20261025_persyn` w25 winner (conductances, drive,
   wsup=0.25 — no retuning in either arm).
2. Arm 1 INTACT: fresh-process 30 s eval. Arm 2 ABLATED: same config
   with ALL afferent + contact feedback gains zeroed (Ia/II/Ib/heel/toe
   classes → 0; efferent/central conductances untouched). Save both runs
   (npz + csv) with a comparison json: cycles/leg, inter-leg phase,
   forward progress, pelvis height, fall yes/no per arm.
3. Plot both traces side by side into the campaign folder.
4. Label arms Deng-2023-style (active / feedback-disabled). Pre-register:
   predict ablated-arm failure mode (fall time, phase collapse, or
   cycling-in-place) BEFORE running.

## Order B — original-experiment reproduction (split with EB475WS4)

B1 (yours, if AnimatSimulator runs there; else mine): headless-run the
ORIGINAL Deng 2019 .aproj (ModelDB 264844 copy at EB475WS4
`Code\Matlab\SNS_Simscape\dev\deng2019_modeldb\K_A wirh RG  march 21\
K_A synergy.aproj`) with the Deng 2023 perturbation protocol: 2 nA to
the RG flexor neuron, 2 s to ~2.x s window; second arm: 2 nA to the knee
PF flexor. Capture RG/PF/MN traces both arms.
B2 (EB475WS4, queued): apply the identical stimulus to the verified
Simulink Deng twin (period 1.938 s).
B3: compare BOTH to the paper's reported behavior — perturbation
sensitivity via the weak cross-excitation is the headline claim; targets
= poster p1 (`_zcode_backups\lit_20261026\poster_p1.png`) + figure
outcomes in kaiyu2023_fulltext.txt.
Li contact feedback, same recipe once B lands: original
`walk tester rearranged.aproj` run with its contact script; same timing
into the Simulink Li twin; compare stance durations + inter-leg phase.
Assets/mine: `deng2019_mine_20261025.json`, `DENG2023_SOURCES_20261025.md`
(both on EB475WS4 under SNS_Simscape\dev).

Rules: every number traceable to a saved run; honest verdicts ("compare
periods, not samples" applies); no silent simplifications; report gaps
rather than claiming reproduction.
