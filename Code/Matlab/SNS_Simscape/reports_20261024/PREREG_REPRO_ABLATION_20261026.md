# Pre-registered predictions — Ben-approved orders A + B, 2026-10-26 (ZCode, EB475WS4)
Written BEFORE any run, per the 2026-10-26 post-hoc correction rule.

## Order B2 — Deng perturbation, Simulink twin (2 nA -> RG flexor, 2.0-2.5 s)
Baseline: self-sustained rhythm, period 1.938 s from the 10 nA start kick.
PREDICTIONS:
- P1: the rhythm PERSISTS through and after the 2 nA window (no arrest);
  bursts continue (paper's model tolerates the perturbation).
- P2: the perturbation shifts burst timing within the window (period
  measured across the window deviates >10% from 1.938 s), i.e.
  perturbation sensitivity WITHOUT destabilization — consistent with the
  weak cross-excitation claim (removing it would make the CPG MORE
  sensitive; we do not ablate it here).
- P3: post-window period returns toward 1.938 s within 2 cycles.

## Order A — feedback-ablation pair, walking W2L model (MuJoCo, desktop harness)
Config: Ben's literal aproj defaults (te 3/tf 4, documented contact2-era
tuning), 30 s fresh-process evals, frozen in both arms.
Arms: INTACT vs CONTACT-FEEDBACK-ABLATED (ch_hc, ch_ank, ch_contra = 0;
plus any afferent-chain knobs present in the harness — if the harness
exposes none beyond contact families, the ablation is labeled
contact-feedback-only, partial).
PREDICTIONS:
- P4: INTACT survives (documented: obj 90.0, 19/17 cycles, phase 0.498).
- P5: ABLATED either falls (pelvis collapse) or loses inter-leg phase
  lock while cycling; my prediction: phase lock degrades first (phase
  exits [0.35, 0.65]) because the crossed heel->contra-RG-flx edge is the
  documented antiphase mechanism — NOT immediate falling, since drives
  and commissural coupling remain.
- P6: forward progress strictly decreases vs intact.

## Order B1 — original Deng 2019 .aproj headless (attempt)
Goal this session: BASELINE run only (does the original project execute
headless and produce its rhythm/stepping record?). The 2 nA stimulus
injection into the original is XML surgery — attempted only if the
baseline run works and budget remains. No prediction beyond: the original
project starts its rhythm from its own kickoff stimuli (per its
architecture) — if it does not run headless unmodified, that is a finding,
not a failure to hide.
