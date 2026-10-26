# easteregg2 prompt — Order B1: run the ORIGINAL Deng 2019 AnimatLab project headless + the 2 nA perturbation (Ben, 2026-10-26)

Paste into a fresh session on easteregg2. Machine guard first: `hostname`
must be easteregg2; read only that machine's AGENTS.md bullet. Ben
authorized reverse SSH to `eb475ws4` (read-only: scp files FROM it; never
write outside its `tmp\`).

## Goal

Complete the first true original-experiment reproduction in the program:
run the ORIGINAL Deng 2019 project (ModelDB 264844) headless in
AnimatSimulator, apply the Deng 2023 perturbation protocol, and compare
three things: the ORIGINAL's response, our Simulink twin's response
(already banked, below), and the paper's reported behavior.

## Assets

- ORIGINAL project (on eb475ws4; scp it over):
  `Code\Matlab\SNS_Simscape\dev\deng2019_modeldb\K_A wirh RG  march 21\
  K_A synergy.aproj` (plus the Back-up and K&A-synergy-1 variants).
  Census: 90 non-spiking neurons, 12 muscles, 196 synapses, 58 named
  synapse types (per-target conductances; per-instance values are the
  untouched 0.5 uS GUI default — type conductance is the effective
  strength; mine: `dev\deng2019_mine_20261025.json` on eb475ws4, plus
  `dev\check_deng2019_instances_20261025.py`).
- Deng 2023 protocol + targets (already local to you):
  `D:\GitHub\_zcode_backups\lit_20261026\kaiyu2023_fulltext.txt` and
  `poster_p1.png`. Protocol: external 2 nA current to the RG FLEXOR
  neuron, window t = 2.0 to ~2.5 s (mine quotes "from 2 s to 2.x");
  second condition: 2 nA to the knee PF FLEXOR. Node names in the
  project (from the mine): "LH_RG F"/"RH_RG F" (RG flexors),
  knee PF flexors under the K&A layer — resolve exact GUIDs from the
  mine json before editing.
- Headless run + XML-surgery patterns: YOUR OWN
  `spinal\reports_20261020\w2l_animatlab\` machinery (five headless
  30.1 s Vortex runs, Root telemetry, RG burst counting; stimulus
  injection = lockstep GUID edits adding TonicCurrent stimuli,
  deterministic prefix — see build_tuned_variant_20261020.py and
  edits_log_20261020.json). AnimatSimulator.exe is on eb475ws4
  (`D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\`);
  if your box lacks it, run on eb475ws4 via SSH instead and say so.

## Our twin's banked result (the comparison target, from eb475WS4)

Simulink Deng twin, same protocol (2 nA -> RG flexor, 2.0-2.5 s):
extensor bursts fully SUPPRESSED during the window (0 bursts in-window);
post-window rhythm resumes at period 1.9377 s vs the 1.938 s reference
(results\deng_perturbation_twin_20261026.mat/.png on eb475ws4). The
question for the original: does it show the same family — perturbation
sensitivity (burst interruption or timing shift) WITH recovery?

## Procedure (in order)

1. PRE-REGISTER predictions to
   `D:\GitHub\_zcode_backups\diss_campaign_20261024\PREREG_B1_deng_original_20261026.md`
   BEFORE any run (program rule since the 2026-10-26 post-hoc
   correction). Include: baseline rhythm expectation, in-window
   behavior, recovery.
2. BASELINE first: run the unmodified original headless (working COPY —
   never edit the original file), capture RG E/F + knee MN traces. If it
   does not run unmodified, that is a finding — report it, do not force
   it.
3. PERTURBATION: XML surgery on the copy adding the 2 nA stimulus on
   the RG-flexor GUID for the 2.0-2.5 s window (your drop-test
   TonicCurrent pattern). Run both conditions (RG flexor; knee PF
   flexor).
4. COMPARE: original vs twin (banked numbers above) vs the paper's
   figures/poster (weak cross-excitation = perturbation sensitivity is
   the headline claim; "compare periods, not samples").
5. Save everything under
   `spinal\reports_20261026_dengoriginal\` (in-repo) with an honest
   verdict. Follow-on if budget remains: Li contact-feedback recipe —
   run `Neuromechanical_Models\...\walk tester rearranged.aproj` (the
   `li` mine's source) headless with its contact script, then the
   Simulink Li twin (eb475ws4) with the same timing; compare stance
   durations + inter-leg phase.

Rules: working copies only; no commits/pushes (Ben commits); every
number traceable to a saved run; gaps reported, reproduction never
claimed beyond evidence; terminology — "time steps" for integrator
steps, "stride/foot steps" for gait, experimental conditions are
"conditions" (the model's limbs are legs).
