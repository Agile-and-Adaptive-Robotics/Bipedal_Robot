# Deng 2023 perturbation in MuJoCo — EB475WS4 replication (2026-10-26)

Prereg: PREREG_DENG_MUJOCO_20261026.md (written first). Gate:
verify_deng.py PASS (1937.8 ms). Runner: deng_mujoco_runner.py (working
copy; w2l_air.xml welded root = suspended hold; LEFT leg only; NO
contact feedback; feedforward coupling — the template carries no Ia/II/Ib
populations; knee/ankle driven from the PF_ka half-centers via the
Fig-6B map, documented substitution). Two build bugs found and fixed en
route, both recorded in the runner header/probe: float-compare kickoff
failure, and a seconds-vs-milliseconds dt mismatch (substeps 1000x too
small). Neural 2 ms cadence with 0.1 ms Euler substeps; physics 1 ms;
12 s runs; data: deng_mujoco_{base,p1_rg,p2_pf,p3_rg_long}.npz.

## Results (burst onsets, −45 mV rising, t>1 s)

- base: ext onsets 1.942…11.630, period LOCKED 1.938 s; flx interleaves
  (2.910, 4.848, …). Coupled period = neural-only (1.938) — expected:
  the coupling is feedforward, so the body answers the circuit, not the
  reverse (prereg P-BASE: correct).
- p1_rg (2 nA -> RG flexor, 2.0-2.1 s): trajectory perturbed (max dV
  2.2 mV) but onset grid moves ≤2 ms; NO extra flexor burst; NO
  suppression; recovery immediate. Prereg P1 phase-reset prediction
  WRONG at this stimulus strength/period ratio — honest scorecard.
- p2_pf (2 nA -> knee PF flexor): trace BIT-IDENTICAL to baseline
  (max dV_ext = dV_flx = 0.000e+00 at the RG). NO RG effect. Prereg P2
  correct.
- p3_rg_long (2.0-2.5 s window): large divergence (max dV 24 mV — the
  flexor burst pattern is materially altered in-window) yet every ext
  burst survives (onsets 3.878, 5.816 inside/post window), no deletions,
  period returns 1.938. Prereg P3 family correct.

## Four-way table (paper / AnimatLab original / Simulink twin / MuJoCo)

| Condition | Paper (biarticular) | AnimatLab original (2026-10-26, egg2) | Simulink twin | MuJoCo (this) |
|---|---|---|---|---|
| 2 nA -> RG flexor | rapid short phase step; rhythm persists | extra SHORT flx burst at 2.000 s (baseline next 2.433); NO ext suppression; post 0.224 vs 0.221 s | 0 ext bursts in 2.0-2.5 s window; post 1.9377 s | no suppression at EITHER window (0.1/0.5 s); dV 2-24 mV; onsets ≤2 ms; post 1.938 s |
| 2 nA -> knee PF flexor | no RG effect (predecessor arch) | −2.3±2.4% (none) | not run | EXACTLY none (0.0 mV) |

## Verdict

1. Knee-PF null: THREE-WAY AGREEMENT (paper predecessor claim, AnimatLab
   original, MuJoCo) — the strongest cross-platform result.
2. RG-flexor stimulus: original + MuJoCo + paper agree on "perturbation
   WITHOUT rhythm suppression, with recovery." The Simulink twin's
   in-window suppression is the OUTLIER; per the prereg decision rule,
   the window-artifact hypothesis is SUPPORTED (MuJoCo at the twin's own
   0.5 s window shows zero burst deletion). The twin result should be
   restated/verified before any prose cites it as suppression.
3. No extra short flexor burst here vs the original's — explained by the
   period ratio (original cycle 0.22 s: a 0.1 s stimulus spans half a
   cycle; template cycle 1.94 s: the same stimulus is ~20x weaker
   proportionally). Reported, not tuned (prereg rule).
4. Scope: this is the Nourse-template antagonist-only lineage coupled
   feedforward to a suspended body — a replication of the PROTOCOL and
   the predecessor architecture's responses, not of the paper's
   biarticular model (which requires BFP/RF, easteregg2's H1).

Files: npz data (4), runner, prereg, probe/diff/check/onset scripts, this
report. No existing files touched; nothing committed.
