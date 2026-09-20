# Gait-variants plan: uneven terrain, walking lunges

## Roadmap (Ben, 2026-09-18): balance posture -> walking -> running ->
## walking -> standing again

The long-range target is a single network that CHANGES GAIT on command.
The joint-layer PF structure (G["joint_pf"], wired 2026-09-18) is the
mechanism Ben expects to use: layers triggered per phase of the
sequence.  How that composes:
  * DRIVE schedule changes the rhythm regime (the network already
    self-sustains 0.4-0.8 Hz; running needs ~2.5-3 Hz = flight phase,
    a genuinely new regime requiring E-burst shortening and likely a
    running-specific ankle push) — this is the hardest transition;
  * PF layer GAIN scheduling (which layers get rg_to_pf drive and
    afferent-loop gain) sculpts the pattern per gait: standing =
    POSTURE-dominant with layers near zero; walking = the T1-fitted
    layer pattern; running = ankle-E + hip-E dominant with KNEE-F
    re-timed;
  * transitions are phase-triggered (heel/toe + hip sensors already in
    the net) rather than clock-triggered.
Prerequisite for ALL of it: the current model dialed in (curriculum
rerun launched 2026-09-18) and, ideally, the scaled-subject conversion
so tuning happens on the right anthropometry.

Ben, 2026-09-18: "I'm interested in different types of gait. Walking
lunges. Walking over uneven terrain. Come up with a plan to implement
that."  Order chosen: terrain FIRST (same gait, tests robustness of the
existing circuit), lunges SECOND (new motion class, needs reference data
and a slower rhythm).

## Prerequisites (gating, in order)

1. A valid ground-walking winner: curriculum stages 2-3 currently hold
   -100 sentinels; re-run them (sentinel purge decision first).  Terrain
   work without a ground walker has nothing to stress-test.
2. Anthropometry fix (recommended, from audit_ik_ground.py 2026-09-18):
   convert the SCALED subject01_simbody.osim with MyoConverter (easteregg2,
   ~90 min) and redo the backsolve chain, OR consciously accept the
   ~9 cm offset.  Reference-activation work for lunges needs the fixed
   model; terrain work does not (it stresses the controller, not the
   reference).

## Stage T1-T3: uneven terrain

* T1 — terrain generation: MuJoCo heightfield XML (seeded RNG): (a) gentle
  undulation (+-1 cm, wavelengths 10-30 cm), (b) discrete steps/blocks
  (0.5-2 cm), (c) ramps.  Three difficulty tiers; fixed seeds for
  reproducibility.  Same rig/pelvis support machinery as ground walking.
* T2 — sensing: the heel/toe mechanosensors are the primary terrain
  sensors; verify they fire on the heightfield (contact normal force vs
  height change).  Consider adding a per-foot earliest-contact probe if
  midfoot contacts dominate.  Ankle PlantarForce/knee compliance: the
  existing joint limits cover it; no model change expected.
* T3 — controller stresses: keep the circuit identical, vary terrain
  tier.  Measure: stays-up fraction over N>=10 seeds x 20 s, pelvis
  height variance, per-step duty variance, trip-recovery count.  If the
  flat-ground winner fails tier (a), the mediolateral/frontal-plane
  pathway (BAL_LAT/LBIN + any future frontal PF bias) is the first
  suspect, not the sagittal CPG.
* T4 — curriculum hook: a `--terrain tier` flag in runner.py; optional
  stage-4 study "terrain_1" seeded from the ground winner with terrain
  randomized per trial (domain randomization).

## Stage L1-L4: walking lunges  (Ben 2026-09-18: CONTINUOUS cyclic gait —
## NOT held poses; corrected from the earlier quasi-static framing)

Walking lunges are a slow rhythmic gait: one lunge cycle ~2.5-3.5 s
(0.3-0.4 Hz) with a deep front-knee bend, a long double-support weight
transfer, and continuous flow into the next lunge.  This is IN the
demonstrated RG regime — the NaP air smoke run already self-runs 2.1 s
periods (0.47 Hz) and Ivanenko air-stepping sits at 0.3 Hz — so no
architecture surgery is expected for the cadence.  The hard parts are
amplitude and balance, not timing.

* L1 — reference motion (no repo data exists): generate a parametric
  CYCLIC lunge-walk trajectory as joint-angle targets on the scaled
  model: alternate lead leg each cycle; front knee to -90..-100 deg,
  rear knee toward -90..-110 deg (model limit -120), pelvis drop
  ~0.2-0.3 m, trunk upright, ~0.3 Hz continuous; cross-check amplitudes
  against a lunge biomechanics paper in Ben's Zotero before finalizing
  (provenance-tag it).
* L2 — activation reference: track the trajectory in MuJoCo
  (equality/actuator tracking as in bsolve's IK replay), run the
  ridge/NNLS backsolve for activations; this gives the lunge synergy
  target on the FIXED anthropometry.  Expect large ankle-PF (push-off)
  and knee-extensor eccentric (controlled descent) features.
* L3 — controller work: slow the RG to lunge cadence (DRIVE down; the
  period knob) with duty ~0.65, larger PF amplitudes for the deep RoM.
  The stress points, in order: (1) the long double-support transfer —
  mediolateral + A/P balance with the CoM far outside the base of
  support (this is where the frontal-plane CPG-vs-reflex question gets
  answered empirically); (2) rear-leg swing clearance past the front
  leg; (3) eccentric knee-extensor control during the descent.
  Fallback (only if the cyclic drive fails): CPG-triggered posture
  sequencing — demoted to fallback, not the plan.
* L4 — success metrics: N continuous lunges without fall, pelvis drop
  reached each cycle, front-knee angle within +-10 deg of target,
  frontal stability during transfers; compare kinematics against the
  L1 reference with the kine-style score.

## Interactions with the joint-layer PF work

The T1 layer fit is WALK-specific by construction.  Lunge/terrain trials
should NOT retune the layer weights; they test generalization.  If
terrain demands a second ankle-E burst (foot placement), that is
evidence for per-layer afferent loop gains (aff_e_pf per layer), not for
restructuring.
