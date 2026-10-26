# PREREG — Deng 2023 perturbation replication in MuJoCo (EB475WS4, 2026-10-26)
Written BEFORE any coupled run, per program rule. Template gate PASSED
first (verify_deng.py: self-sustained, 1937.8 ms, anticorrelated).

## Banked comparisons (the targets)

- Paper (Deng 2023, biarticular model): rapid short phase step; 15.44%/15%
  phase-delay class of effects; no RG effect from knee-PF stimulus in the
  antagonist-only predecessor architecture.
- AnimatLab ORIGINAL (easteregg2, 2026-10-26, ModelDB 264844): 2 nA on RG
  flexor t=2.0-2.1 s -> EXTRA SHORT FLEXOR BURST at 2.000 s (baseline next
  onset 2.433 s), NO extensor suppression, post-period 0.224 vs 0.221 s
  (phase-reset-with-recovery). Knee-PF stimulus -> no RG effect
  (-2.3 +/- 2.4%). NOTE: original's period ~0.22 s (its own parameter
  set), NOT the Nourse-template 1.94 s.
- Simulink twin (this machine, 2026-10-26): 2 nA RG flexor t=2.0-2.5 s ->
  0 extensor bursts IN-WINDOW, post-period 1.9377 vs 1.938 s. Window is
  0.5 s = 26% of the twin's period (vs 4.5% of the original's) - the
  suppression reading may be a WINDOW ARTIFACT.

## Predictions (MuJoCo coupled, suspended hold, Nourse-template circuit)

- P-BASE: coupled air-stepping alternation sustains; coupled period
  differs from the 1.94 s neural-only value (body loading shifts it) -
  reported, not tuned. Range guess: 1.6-2.4 s.
- P1 (2 nA RG flexor, 2.0-2.1 s): phase-reset-with-recovery like the
  ORIGINAL - the stimulus lands inside one cycle (10.2 s window share),
  advancing the next flexor onset; extensor rhythm NOT suppressed;
  post-period returns to baseline within 2 cycles (±3%).
- P2 (2 nA knee PF flexor): small or no RG effect (<5% period/phase
  change) - matches the original and the paper's predecessor result.
- P3 (2.0-2.5 s window variant): expect the same qualitative outcome as
  P1 with a longer tail; this arm exists to test the twin's window
  artifact hypothesis directly.
- DECISION RULE: if MuJoCo shows phase-reset-with-recovery on P1/P3 while
  the twin showed in-window suppression, declare the twin reading a
  window artifact and quantify (expected flexor onset advance vs window).

Violations will be reported, not retrofitted.
