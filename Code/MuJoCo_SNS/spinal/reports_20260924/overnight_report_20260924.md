# Overnight report — 2026-09-24 (EB475WS4)

Compiles five finished work packages from `reports_20260924\`: `scone_hyfydy_hands_on.md`,
`animatlab_w2l_airwalk.md`, `gait2392_muscle_function_pf_assignment.md`,
`balance_data_catalog.md`, and the new `spinal\w2l_cpg\` package. A gatekeeper audit re-ran
the headline checks (fresh w2l_cpg smoke reproduced `W2L_SMOKE PASS` bit-exactly; all quoted
numbers re-derived from the delivered logs and chart files) and shipped everything.

## 1. TL;DR — your two questions, answered per model

**"Did the models stand with proprioceptive + vestibular feedback?" — No, not anywhere
tonight, and honestly so.**

- **SCONE Tutorial 3a (reflex + vestibular, OpenSim):** at default (unoptimized) parameters
  it falls at t = 1.05 s — COM drops to 59.9% of initial, tripping `termination_height = 0.6`.
  This is the tutorial working as designed: the vestibular gains start at 0
  (`$KP`/`$KV = 0~0.1`) and CMA-ES is supposed to find them. Nobody has run the optimization
  yet.
- **AnimatLab W2L / BilateralRG:** no balance layer at all — the body is suspended/frozen;
  feedback untuned (your own words, confirmed).
- **MuJoCo s3k walker:** not re-run tonight. Its known state stands: real walking only inside
  the rig (weaning boundary S ≈ 0.8–1.0); the pelvis-balance piece is still the missing
  architecture.
- **w2l_cpg (new port):** network only — no body, nothing to balance.

**"Walking in air?" — Yes, four times over.**

| Model | Air-walk? | Numbers |
|---|---|---|
| AnimatLab W2L 2023 asim | YES, cleanest | 8 bursts in 10 s, ~0.77 Hz, L/R antiphase (+0.457 cycle), stays suspended at 0.96 m |
| AnimatLab W2L modern asim | YES, sub-threshold | coherent 2.222 Hz alternation, hip 38° / knee 61° / ankle 16° swings — but peaks −55.8/−56.4 mV, just under the −55 mV burst threshold |
| AnimatLab BilateralRG modern | YES, above threshold | RG crosses −55 mV: 10–11 bursts, 0.464 s (2.16 Hz), antiphase, joints step at the same period |
| w2l_cpg (SNS-Toolbox port) | rhythm only | PASS: period 1.000 s, antiphase r = −0.532, 22+22 RG-E bursts, locked to alternating 1 Hz heel-contact trains |

## 2. SCONE + Hyfydy (`scone_hyfydy_hands_on.md`, logs in `scone_logs\`)

Four headless `sconecmd` runs (SCONE 2.4.4.3333) on your Tutorials3 copies; sconestudio.exe
untouched, no keys reproduced anywhere.

- **Balance–OpenSim at defaults: falls** at 1.05 s (result 102.92 = BalanceMeasure 96.5;
  COM 0.965 → 0.578 m, drifts 0.76 m backward, pelvis pitch −10° → +54°). BalanceMeasure
  scores 0 whenever the model is still up when `max_duration` ends (proven by the A/B run).
- **Gait–OpenSim at defaults: STEPS** — 7 steps, 2.50 m in 2.78 s (0.523 m/s vs the 1.0 m/s
  target), then falls (COM 84.9% < 0.85; hips slam to −1.55/−1.31 rad; knee limit-torque
  penalties 7.42/5.65 > 1). Runs 1.17× real time.
- **Controller reading (banked for our goal-2 VEST stage):** proprioception = MuscleReflex
  arcs U = C0 + KF·[F−F0]⁺ + KL·[L−L0]⁺ + KV·[V−V0]⁺ with tiered delays (hip 10 / hamstrings
  15 / knee 20 / ankle 35 ms, vestibular 100 ms), negative-gain cross-muscle arcs as
  antagonist inhibition, phase gating by a 5-state FSM on leg-load + foot position; balance =
  PD on a torso body point broadcast to all 9 muscles, or phase-separated DofReflex on
  pelvis_tilt around P0 = −0.105 rad. Their `termination_height` 0.6 (stand) / 0.85 (gait)
  are ready-made acceptance metrics.
- **Hyfydy blocker + the 1-minute fix:** `Tutorial 4a - Gait - Hyfydy.scone` fails before any
  sim step — "no active license key was found" (exact text in `gait_hyfydy.log`; it's a trial
  license that isn't activated, not a SCONE bug). Fix: open SCONE Studio → **Tools →
  Preferences → Hyfydy → paste your license key → Enable** (or run sconecmd with
  `--hyfydy <key>`). The Hyfydy scenarios are testable the moment you do that.
- **Gotcha banked:** scenario overrides need the `CmaOptimizer.` root
  (`CmaOptimizer.SimulationObjective.max_duration=N`); the short path is silently ignored
  with an "unused properties" warning — proven by A/B runs (1.000 s, BalanceMeasure 0 vs the
  warning and a 1.05 s fall).

## 3. AnimatLab results (`animatlab_w2l_airwalk.md`)

`AnimatSimulator.exe` IS installed here (registry-found; it takes no flags). Four headless
runs, all exiting cleanly. **Your claim is CONFIRMED — walking in air works — with one
nuance:** the modern W2L export's flexor half-centers never cross −55 mV because that export
carries zero tonic drive (its only stimulus is a 10 ms 1e-8 current tick at t=0, vs the 2023
model's 14 neurons held at 5–6 nA). The rhythm is real but sub-threshold: ~5–6 mV swings
peaking at −55.8/−56.4 mV, coherent across all charted neurons at 2.222 Hz, half-centers
antiphase (r = −0.956), L/R hip flexors antiphase (r = −0.994; mid-swing lag −0.492 cycle),
joints stepping (hip −15..+23°, knee to +60°, ankle −20..−4°), neural leading joint by
57–89 ms (r = 0.90–0.92). The 2023 asim is the reference air-walk (above). BilateralRG's RG
does cross −55 mV and its ground variant shows real intermittent foot contact (toe_R mean
0.49 / max 2). Honesty notes: apply the −55 mV criterion in chart units (−0.055 — charts
store volts, rest −60 mV); the modern W2L RG chart logs the LEFT side only. No model file was
modified — only chart .txt byproducts were rewritten (backups kept next to the report); your
open GUI and sconestudio were left alone.

## 4. New: w2l_cpg — Ben's AnimatLab CPG ported to SNS-Toolbox (`spinal\w2l_cpg\`)

A self-contained package (`build_w2l_net.py`, `smoke_w2l.py`, `README.md`) mechanically
transcribing your `bilateralrg` connectome template (105 nodes / 226 edges, SESSION_NOTES
build chain) — it designs no topology; every synapse keeps its template tag. Builds 91
neurons / 214 synapses / 10 inputs / 12 muscle output ports; NaP half-centers use the
toolbox class with the fixed-τh backend (τ 0.25 s) since stock tau_h(V) quenches them;
AnimatLab synapse gains re-expressed on the 0–5 mV scale per pathway (table in the README),
all overridable. Fresh smoke (deterministic, gatekeeper-reproduced):
`W2L_SMOKE PASS period=1.000 antiphase_r=-0.532 rg_e_bursts_L=22 rg_e_bursts_R=22`.
One deviation, documented: the template's R side is missing the 2 direct half-center
excitation edges the L side has (a `make_editor_templates.py` generation gap; SESSION_NOTES
says R mirrors L) — added with L's own gain (asserted 214 = 212 + 2). **Deliberately NOT
done:** no hookup to the 92-muscle MuJoCo model yet (muscle nodes are bare output ports),
afferent channels silent in the smoke, no tuning campaign — existing tuned winners untouched,
nothing outside `w2l_cpg\` modified.

## 5. Muscle-function design doc (`gait2392_muscle_function_pf_assignment.md`)

All three source tracks fetched and read (OpenSim Confluence Gait2392 page + its two PDFs
incl. per-muscle isometric forces; John et al. 2012 mediolateral-GRF paper full text; repo
cross-check of muscle_map.py / build_network.py / runner.py / params.py with line cites; 92
actuators verified in the cvt3 MJCF), every claim tagged [FETCHED]/[REPO]/[LIT-STD], with a
per-group table for both legs. **The adductor/abductor PF-pair scheme in 5 sentences:** add
per-side PF_ADD_/PF_ABD_ half-centers, both driven from RG_E as a staggered stance pair,
because John 2012 finds abductors the largest medial contributor in all three ML-GRF windows
and adductors lateral, peaking late stance. Cross-inhibition is laminated through new
PF_IN_ADD_/PF_IN_ABD_ interneurons — no direct PF-PF synapses, matching the existing PF
build pattern. MN membership is 6 adductor pools (agonist) vs 6 gluteus medius/minimus pools
(antagonist), with add_mag3/sartorius/TFL ride-alongs excluded. Per-pool Renshaw is already
generic, and the whole thing sits behind `G['front_pf']` default-0 conditional topology plus
JSON-RULE loader branches, so defaults stay bit-identical. Phase-mode and joint_pf-mode
variants are specified alongside. Delivered as DESIGN — no code touched.

## 6. 3-DoF foot sensor (same doc, DESIGN)

Per contact, take MuJoCo's `mj_contactForce` 6-vector and rotate it to world with the
contact frame rows (`c.frame`): normal = buf[0], tangents = buf[1:3] → signed
mediolateral / fore-aft / vertical channels per foot. This extends the existing
HEEL_c/TOE_c/LOAD_c ports without touching their schema (the runner currently reads the
normal component only), and each channel has a natural consumer: vertical → gait phase/load,
ML → the new frontal-plane PF pair, FA → shear. The doc includes a MuJoCo-side sign audit
and a capture-completeness recipe before anything is trusted.

## 7. Balance-data shortlist (`balance_data_catalog.md`; licenses verified live via APIs)

1. **Wang & van den Bogert 2020** standing-balance random pulses — CC BY 4.0 (Zenodo
   3819630): 8 subj × 80 min perturbed, markers + 6-DOF GRF + 9 EMG, and it SHIPS processed
   angles/torques (no IK needed) — flagship for balance training refs, feeds --stand-eval/VEST.
2. **GaitRec** — CC0 (figshare c.4788012): the only dataset with explicitly named
   mediolateral GRF items; 2,084 patients + 211 controls / 75,732 trials — flagship for
   frontal-plane validation.
3. **Schreiber & Moissenet 2019** — CC BY 4.0 (383 MB): 50 adults × 5 speeds with kinematics
   + 3D GRF&M + EMG together — best experimental upgrade to our predicted speed sweep.
4. **Moore et al. 2015** perturbed walking — CC0 (Zenodo 13030): 120 s normal + 480 s
   perturbed per trial at 3 speeds; AP-only perturbations, no EMG.
5. **Gutenberg Gait Database** — CC0 (~100 MB): 350 healthy subjects, vertical + ML GRF +
   COP — cheapest normative frontal-plane win.

Verified local gaps: zero perturbation references and zero mediolateral GRF (the loader
reads vertical columns only; refs store no GRF waveforms) — M/L validation needs a loader +
npz-schema + sim-side contact-wrench extension (notes in the report). Caveats: Wang's
perturbation direction isn't stated in the record; figshare listings are asymmetric on side
— check both at download.

## 8. Open items for Ben

1. **Hyfydy:** paste your license key (Tools → Preferences → Hyfydy → Enable) — one minute,
   then the Hyfydy tutorials run.
2. **w2l_cpg:** OK the R-side mirror fix (and whether to patch `make_editor_templates.py` so
   the template itself carries the 2 missing edges); say when to wire the 12 muscle outputs
   to the gait2392 actuator map.
3. **PF pair + 3-DoF foot sensor:** DESIGN awaiting your go (alongside the goal-2 VEST
   connectome-spec confirmation from 09-23).
4. **Balance data:** pick downloads from the shortlist (Wang + GaitRec recommended first);
   integration needs the loader/npz/contact-wrench extension.
5. **Standing status:** nobody stands tonight without a rig — SCONE needs its CMA-ES run to
   find the vestibular gains, and our stack still needs the pelvis-balance piece; SCONE's
   delay/gain priors (sec. 2) are banked for that design.
