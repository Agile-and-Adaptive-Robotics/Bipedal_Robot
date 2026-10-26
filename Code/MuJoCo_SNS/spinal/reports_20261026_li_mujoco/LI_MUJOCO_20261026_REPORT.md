# Li follow-on in MuJoCo — lever-1 descent fix (2026-10-26, egg)

Ben's redirect of the Order-B1 follow-on (2026-10-26): do the Li contact
feedback recipe in MuJoCo on this box. Preregistered before any run:
`D:\GitHub\_zcode_backups\diss_campaign_20261024\PREREG_LI_MUJOCO_20261026.md`.

## VERDICT (honest)

**The 2026-10-01 never-contacts dead end is FIXED by lever 1, and the Li
CPG's engaged side steps at Li's own reference period — but only ONE leg
engages, and the left-leg blocker is NEURAL, not geometric.**

1. **Lever 1 works.** Lowering the descent/retained-support PD target below
   the keyframe contact plane (`--drop=0.05 --sink=0.01`, the W2L-proven
   recipe; working copy `test_li_stepping_zdrop.py`, Ben's runner untouched)
   takes the pelvis from the 2026-10-01 never-descends state (0 stance
   episodes in all seven prior runs) to ground contact and R-foot stepping
   in every configuration tested.
2. **The R-foot step period sits on Li's AnimatLab reference.** At support
   0.5 and 0.4: 1.294 s and 1.295 s (Li ref 1.30 s; duty 0.07). The
   0.658–0.724 s readings at support 0.7 are the same rhythm sampled at
   half-period by the episode splitter, not a different rhythm.
3. **The left foot never loads — and that is an activation problem, not a
   geometry problem.** Eight active configurations (drop 0.03–0.10, support
   0.4–1.0, kick on R, on L, or on both) all give 0 L-foot episodes: the
   active L leg holds its knee at −66…−26° while the R knee extends to
   +7.6°. The passive probe (`--passive=1`, muscles zeroed) drops the same
   body onto the same platform with the L knee relaxing to −10…−4° and the
   L foot resting in contact for 15 s: the leg REACHES fine. The loop is:
   the L stance CPG needs heel input to phase into extension, and the L
   heel only fires if the knee extends. The R side broke this loop through
   its swing trajectory; the kick side is NOT the determinant (kicking L
   still steps only R).
4. **Forward speed stays ~0.04 m/s vs Li's 0.638** — expected and predicted:
   his original runs on a virtual walker; this protocol retains the
   harness, so the claim is rhythm, stance, and phase, not distance.

**Prereg scorecard:** P2 HELD (R period within ±30% of 1.30 s — measured
within 0.5%); P1 HALF-HELD (R ≥3 episodes in every run; L 0 everywhere);
P3 UNTESTABLE (antiphase r needs both feet); P4 as predicted; the
"contact without rhythm engagement" falsifier applies to the left leg
specifically, with the passive-probe diagnosis above.

## Comparison of record

| metric | Li AnimatLab original (his asim, DataTool_7) | Simulink twin (2026-09-30) | this MuJoCo run (best: drop .05, support .5) |
|---|---|---|---|
| step period L / R | 1.299±0.024 / 1.333±0.080 s | 1.000 s (contact-locked) | L: none / R: 1.294 s |
| stance duration | 0.61–0.67 s | (period claim only) | 0.053 s (toe-tap; not loaded stance) |
| duty | 0.443 / 0.555 | — | 0.07 (R) |
| L/R phase | antiphase | r ≥ 0.9999 with contacts | untestable (L idle) |
| forward speed | 0.637 m/s | — | +0.043 m/s (harness-retained) |

The R-period agreement (1.294 vs 1.299±0.024 s) is the strongest
cross-engine number this stack has produced for the Li circuit: the
MuJoCo-coupled contact-driven CPG, once it engages, runs at the original's
period. Stance duty is NOT yet comparable — the R touches are brief heel
taps, not loaded stances, so the contact encoder fires but the body rides
the harness.

## Runs (all logs in this folder)

| log | knobs | R episodes / period | L |
|---|---|---|---|
| li_zdrop_d05.log | drop .05, support .7 | 14 / 0.658 s | 0 |
| li_zdrop_drop003.log | drop .03 | 5 / 3.530 s | 0 |
| li_zdrop_drop007.log | drop .07 | 5 / 2.916 s | 0 |
| li_zdrop_drop010.log | drop .10 | 14 / 0.658 s | 0 |
| li_zdrop_d05_s05.log | drop .05, support .5 | 11 / **1.294 s** | 0 |
| li_zdrop_d05_s04.log | drop .05, support .4 | 8 / **1.295 s** | 0 |
| li_zdrop_d05_s05_kL.log | kick on L | 8 / **1.310 s** | 0 |
| li_zdrop_d05_s05_kB.log | kick both | 12 / 0.724 s | 0 |
| li_zdrop_d05_passive.log | muscles zeroed (probe) | 1 touch / — | **15 s resting contact** |

Runner: `w2l_mujoco\test_li_stepping_zdrop.py` (working copy; knobs
`--drop --sink --kickside --passive` added, all documented in-file; Ben's
`test_li_stepping.py` untouched).

## Next lever (needs Ben's choice, not silently done)

The single missing event is the FIRST L heel strike. Candidates:

1. **Brief bilateral extensor bias during the descent phase only** (land,
   both knees extend momentarily, both heels fire, then the CPG owns it).
   This is the cousin of the settle-assist co-contraction that Ben's
   2026-09-30 protocol ruling REMOVED as "an invention, not the protocol" —
   so it needs his explicit OK, scoped to the descent phase only.
2. **Lever 2 (body fidelity)**: stiff solimp on the 8 hinges + Kpe as
   tendon stiffness (the 2026-09-30 probe recipe) so a stance leg bears
   load; with support reduced further, both feet must load. Bigger job,
   no protocol questions.
3. **Contact-gated descent**: keep descending until BOTH heel encoders have
   fired once, then hold. Protocol-clean (a descent rule, not a neural
   invention), but the pelvis bottoms on the R leg before the L heel
   reaches under activation, so it likely needs 1 or 2 with it.

## Not claimed

- No bilateral stepping, no antiphase number, no stance-duration comparison
  beyond the table above; the AnimatLab original's 0.637 m/s is not claimed
  comparable to the harness-retained protocol.
- The 2026-10-01 air-gate breakage (test_w2l_air FAIL airborne after the
  freejoint rebuild) is unrelated and untouched.

— Li follow-on executor, 2026-10-26
