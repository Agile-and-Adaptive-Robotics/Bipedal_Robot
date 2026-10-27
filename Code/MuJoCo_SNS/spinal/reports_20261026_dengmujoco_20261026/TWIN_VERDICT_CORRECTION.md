# CORRECTION — the Simulink twin "suppression" was an analysis artifact (2026-10-26)

Root cause, proven: the ODE baseline's extensor peaks fall at 1.954,
3.893, 5.831, ... s — ZERO peaks occur in [2.0, 2.6] even unperturbed
(ode_pert_groundtruth.py + the baseline run above). The twin runner
counted prominence peaks in that window (2 pre / 0 during / 4 post) and
I read "0 during" as suppression without a baseline peak count. The twin
was correct all along: its peaks (kickoff + 1.95, then 3.89.. 9.71) are
all present, period 1.9377 vs 1.938 s reference.

## What this means

- SNS-Simscape is NOT broken on this test. The neuron/synapse models are
  independently verified (units test 6.07e-04 mV; tutorial twins <=1.1e-05
  mV; 410-neuron network 4.2e-06 mV vs the toolbox backend), and the
  Deng twin's perturbation response MATCHES the ODE and MuJoCo once
  analyzed correctly: perturbation, no suppression, full recovery.
- Four-way agreement now holds for the RG-flexor condition across
  implementations of the Nourse-template lineage (ODE / Simulink /
  MuJoCo), and with the AnimatLab original and paper qualitatively
  (perturbation without destabilization). Knee-PF null: paper, original,
  MuJoCo (exact), twin untested for that condition.
- The earlier statement "twin's 0 in-window bursts" in
  deng_perturbation_twin_20261026 outputs and the MuJoCo REPORT verdict
  is superseded by this note: the artifact was peak-slot placement, not
  window length, and not any model.

## "Everything from AnimatLab verified in MuJoCo" — current map

| Published model | AnimatLab original | MuJoCo verification |
|---|---|---|
| Deng 2019 predecessor circuit (Nourse tables) | run by egg2 (period 0.221 s, its own params) | DONE: this folder (predecessor protocol, 3/4-way agreement) |
| Deng 2023 biarticular model (BFP/RF/GA) | project on hand | NOT BUILT (egg2 H1; needs biarticular actuators) |
| Li contact-driven CPG | aproj on hand | INDIRECT: the W2L MuJoCo lineage is Li-derived and runs; a DIRECT li-template MuJoCo run (test_li_stepping.py) is a bounded follow-on |
| Shevtsova / Rybak | (no AnimatLab originals; standalone python) | out of MuJoCo scope by construction |
