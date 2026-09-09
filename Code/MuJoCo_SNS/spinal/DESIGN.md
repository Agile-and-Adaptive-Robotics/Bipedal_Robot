# Spinal network for Gait2392 in MuJoCo (SNS-Toolbox) — design & status

Built 2026-09-09. Files in `Code\MuJoCo_SNS\spinal\`.

## Goal

Two-level spinal cord network (McCrea–Rybak RG + PF) with proprioceptors for
each of Gait2392's 92 muscles, driving the MyoConverter MJCF of
`gait2392_simbody` in MuJoCo: standing → walking → standing, speed via
descending drive AND reflex-gain modulation, with PF→MN weights back-solvable
from OpenSim IK/SO activation patterns.

## Architecture (per side; ~406 neurons total, 368 inputs)

| layer | cells | notes |
|---|---|---|
| descending | DRIVE, POSTURE, BAL_PF, BAL_DF | MLR/postural surrogates, external inputs |
| RG | RG-E, RG-F + ADAP-E, ADAP-F | half-center: mutual inhib. + slow (0.9 s) self-adaptation; frequency rises with DRIVE; stance-biased drive split (E stronger) sets duty |
| PF | PF-E1, PF-E2, PF-F1, PF-F2 + PFA-* | forced mode (no own rhythm); (tau, adapt) shape multipliers stagger their windows; reciprocal inhib. on conflicting pairs |
| MN+afferents | per muscle: MN, Ia, II, Ib | Ia ~ tendon velocity (pure-signal: no resting tone — baseline drove constant reciprocal inhibition and crushed flexors); II ~ length (small baseline); Ib ~ force |
| Ib load sharing | IBEXC per stance group | Ib afferents → group IN, gated by RG-E → excites that group's MNs (stance reflex reversal / load sharing) |
| balance | BAL_PF/BAL_DF | pelvis-COM PD → ankle PF/DF + hip flex/ext MNs (ankle+hip strategy), fades with DRIVE |

Standing uses a **solved posture pattern**: static optimization at the
keyframe (torque rows = gravity − contact − passive, plus co-contraction
preload; activations regularized toward a physiological prior) injected as
per-MN bias currents (`POST_<muscle>`). This is the same machinery as the
planned IK/SO back-solve, applied to standing.

Ia/II/Ib gains are speed-modulated presynaptically (`params.MOD`, scaled by
DRIVE) and phase-gated from RG potentials — the "reflex modulation = speed
control" knob.

## Literature grounding (verified; mostly in Ben's Zotero)

- **RG+PF two-level**: Rybak et al. 2006 J Physiol (PF = multiple
  populations exciting synergist MN pools, forced by RG, low persistent Na);
  McCrea & Rybak 2008; Rybak et al. 2015 eNeuro review. The 2006 model
  explicitly anticipates multiple phase-shifted PF units for bifunctional
  muscles → Ben's "similar but phase-shifted PF" observation.
- **Speed**: increased drive shortens stance, swing ~constant (Danner 2017,
  Rybak 2024 eLife RP98841 — both RG-only, no PF). **Reflex gain modulation
  as speed control**: Bunz, Ijspeert, Schmitt 2026 Sci Rep
  10.1038/s41598-026-48509-z (offline+online, walk↔run); Ross/Rybak-lab mouse
  spindle-gain work (2018, in Zotero).
- **Precedent for exactly this stack**: SNS-Toolbox paper (Nourse et al.
  2023, Biomimetic 8:247) MuJoCo case study = rat hindlimbs, two-layer RG+PF
  CPG (Deng-style), PF shapes joint commands, MNs adjusted by Ia/Ib, speed
  set by RG inhibition level. Our network is the human-92-muscle analogue.
- **Synergies**: Ivanenko et al. 2004 (5 modules, speed-invariant);
  Kibushi et al. 2018 (speed carried by timing/amplitude); Avaltroni 2024
  spinal maps. → back-solve design in `fit_synapses.py`.
- **Closest competitor**: Di Russo, Ijspeert, Bouri 2023 J Neural Eng
  (CPG+reflex human neuromechanical model) — read before claiming novelty.
- Also in Zotero Sensory Afferent Database: Markin et al. 2016
  (large-scale PF→MN neuromechanical cat model), load-receptor review 2000,
  stick-insect multisensory control, air-stepping by tonic drive (2009),
  human balance NN papers (2017/2000).

## Verified status (this session)

- `check_rhythm.py`: half-centers oscillate; DRIVE=4 nA → period 0.885 s,
  E duty 0.36, left–right antiphase corr −0.88; PF E/F windows alternate;
  within-arcade staggering small (v1).
- Network builds/compiles/steps at dt=5 ms in ~3 ms/step (numpy backend).
- Runner closes the loop on the full 92-muscle model with logging/summary.

## Open problems (in priority order)

1. **Joint-axis sign audit** — `rect_fem` (hip flexor) emerges as the
   solver's hip-*extension* provider and drives at ctrl=1.0 in every standing
   solve: strong evidence the MJCF hinge axis for hip_flexion (and possibly
   others) is inverted vs the OpenSim convention assumed in
   `muscle_map.GROUPS`/`W_PF_MN`. Audit every hinge's axis sign against
   OpenSim coordinate conventions before any further behavioral tuning.
2. **Leg-DoF explosions / NaN** even with the pelvis rig: ankle angles spin
   to 1e4–1e5 deg. `balanceinertia` applied; next suspects: contact
   solref/solimp at dt=5 ms, ankle/subtalar damping,ctrl slew limits.
3. **Pathpoint weld** — OpenSim conditional pathpoints arrived as massless
   slide bodies (singular mass matrix) driven by equality polycoefs; we weld
   them at keyframe geometry (`apply_harness`). Proper fix: per-pose spline
   pathpoints or re-conversion with MyoSuite's own handling (check how
   `myosuite_gait2392_simbody` deals with them before re-deriving).
4. **Standing without the rig** — needs (1) then an actual balance
   controller (ankle+hip strategy exists; COM-y/vel and foot-placement
   absent). Swing-bench / rigid-rig first, free balance later.
5. Stance duty 0.36 (human ~0.6): load-receptor prolongation of stance
   should improve this once 1–3 stop the explosions; can also bias ADAP-E.
6. Fit `W_PF_MN` from OpenSim SO activations (`fit_synapses.py`), then
   reflex gains per speed from multiple-speed trials (Bunz-2026-style
   modulation study).
