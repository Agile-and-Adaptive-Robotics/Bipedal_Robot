# Executive summary — seven goal tracks, 2026-09-23 (for Ben)

All seven tracks completed their deliverables. Goal 1 landed rubber-band multi-select +
group move in `connectome_block_editor.html` v2.2 (not the task-named rules/gains editor)
and the verification probe also found and fixed genuinely broken v2.1 synapse gestures.
Goal 2 adds a default-off standing-balance stage (VEST pathway + II length loop) with its
supervisor block resolved (stage-4 contact-symmetry bias → the abs() form, re-verified on
disk this session); no optuna study was launched. Goals 3–5 are analysis: MuJoCo 2.3.7
cannot natively do elastic tendons, exact Hunt-Crossley contact, or adaptive integration
(two of three admitted in MuJoCo's own docs), and any model change voids every recorded
winner — keep 2.3.7/2 ms, run an offline convergence ladder, revisit SCONE+Hyfydy
post-deadline; the optuna review measured 22–70% sentinel waste and answers your ensemble
question "yes for seeds + small plant jitter, no for synthetic IK perturbation" (we hold
exactly one reference cycle per leg); SimTK downloads are now login-gated, so nothing was
fetched. Goal 6 produced first-pass Rybak/Shevtsova/Shinohara circuit drafts with every
edge cited, both supervisor blocks resolved, and a short list of genuinely-NEW edge
classes. Goal 7 found — did not recompute — the six synergies of record, replayed them at
pooled VAF 0.945, with a clean 6 s suspended open-loop demo and honestly poor tracking.
Nothing touched the dissertation tex, `spinal_run.npz`, `optuna_walk.db`, or your
connectome spec.

## 1. Connectome editor group-select — `connectome_block_editor.html` v2.2

**What changed.** Rubber-band multi-select on empty-canvas drag; group drag preserves
relative offsets and mutates each node's `n.x`/`n.y` exactly like single drags, so
Export/undo/localStorage all persist group moves. SHIFT+click toggles selection (the v2.1
two-click synapse arming moved to ALT+click); Esc clears. Implemented in the **block**
editor because the task-named `connectome_editor.html` has no canvas, palette, or
bLoad/bClear — every concrete anchor in the task maps only to the block editor. The rules
editor (whose export `runner.py` consumes), `CONNECTOME.md`, and `connectome_gains.json`
were left untouched (mtimes 09/22, re-confirmed on disk this session).

**Key finding.** The baseline probe showed the on-disk v2.1 synapse gestures were broken:
one-gesture SHIFT-drag created **no edge** (listener never attached), and the documented
SHIFT-click-A-then-click-B flow created a **duplicate edge** (2× A→B). v2.2 makes both
gestures self-contained — one drag = one edge, two-click = exactly one.

**Verification.** `node --check` exit 0; 37/37 behavioral suite against a DOM stub (band
hit-test, group drag, snapshot/persist, undo, palette/bLoad/bClear, templates, tab
persist); full 255-line diff re-read — zero hunks touch the palette/load/save/clear/
template/serialize paths.

**For Ben to decide.** (a) Confirm the file choice was right. (b) The pre-existing
localStorage lag (active tab's `data` only refreshed on tab switch) has a 2-line fix —
left alone without your OK. Known non-goals: no group delete; band replaces selection.

## 2. Balance curriculum stage — default-off stage 4 + VEST pathway

**What changed.** `_curriculum.py` stage 4 ("balance") plus conditional circuitry:
VEST_r/l cells (tau 0.1 s stands in for the SCONE vestibular delay) → ipsilateral
extensor MNs (exc, `vest_ext`, 42 edges = 21 muscles/side) and flexor MNs (inh,
`vest_flex_inh`, 26 edges = 13 muscles/side — the counts differ only because the
extensor MN pools are larger); a stance-gated II length-loop presynaptic boost
(`vest_prop`, the SCONE KL analog — Ia deliberately untouched, flagged);
`runner --stand-eval` with new `bal_*` metrics (sway, tilt, contact symmetry, COM
height, fall). Modeled on SCONE Tutorial 3a + Di Russo 2023 eq. 5.

**Defaults-off proof + smoke.** Params diff = exactly three new G keys at 0.0; topology
fingerprint identical (410/376/1186/Σg 668.17); vest-ON build arithmetic matches (412
neurons, +68 edges → 1254 synapses). A 2 s smoke through the exact stage-4 path: sway
0.1205 m, tilt 26.94°, contact_sym 0.271 (the normal.mot asymmetric start pose), no fall,
no NaN; output isolated to `goal2_smoke_balance.npz`.

**Supervisor block (round 2) resolved.** The stage-4 symmetry term was `−40·(0.5 − sym)` —
asymmetric: right-side overload earned a score **bonus** up to +20. Now
`−40·|0.5 − sym|` (either side overloading costs the same, 20 max). Re-verified this
session: `findstr` shows the abs() form at `_curriculum.py:279`. Defaults-off unaffected;
no study launched (`optuna_walk.db` untouched, 09/22 5:15 PM).

**For Ben to decide.** The six-item connectome-spec confirmation list in goal2 §6
(VEST cells/edges/signal, vest_prop + the Ia-untouched call, stage-4 weights + searching
`rig_scale`), and when to run `_curriculum.py 4 30` (not launched).

## 3. MuJoCo muscle fidelity — analysis only, nothing edited

**Findings (all measured this session on cvt3/simbridge2, mujoco 2.3.7).** All three
Hyfydy differentiators are real gaps, two admitted verbatim in MuJoCo's docs ("We assume
inelastic tendons…"; fixed-step-only). Muscles carry exactly 1 state (actnum=1); all 92
tendons have stiffness=damping=0; only 4 fixed-step integrators exist (mjmodel.h:127-132).

- **SEE:** impossible natively. A demonstrated extra-DOF workaround (slide joint + spring
  tendon; 1000 N muscle step) is numerically stable at the 2 ms timestep only with a
  **1 kg** virtual tendon mass — and even then the load-force transient peaks at
  16,300 N ≈ 16× the step input (ringing). A pathpoint-like **0.1 kg** mass diverges at
  2 ms (4.6×10⁸ N) and only becomes stable again at 0.5 ms (10,760 N transient).
  Upgrading MuJoCo does not buy SEE and would break the Simulink bridge.
- **Contact:** depth-dependent stiffness **and** damping are native via `solimp`;
  negative-`solref` gives direct k,b. Exact Hunt-Crossley and viscous friction are
  impossible (no custom contact callback). **Adaptive integration:** does not exist;
  hand-rolling it dissolves the 2 ms lockstep and determinism.

**The standing warning.** ANY model change (muscle/contact params, timestep, integrator)
invalidates every recorded winner, optuna study, basin-gate dump, Simulink E1/E2 parity,
and the regenerated dissertation figures — a re-tune-everything cost.

**Recommendation.** Keep 2.3.7 + implicitfast + 2 ms for all dissertation work; buy
confidence with an offline 1/0.5 ms convergence ladder (ensemble statistics, not
pointwise); revisit SCONE+Hyfydy post-deadline as the standing-balance engine (installed
here, gated only on a free Hyfydy key; port = SconePy co-sim rewrite + `.hfd` conversion +
full retune).

**For Ben to decide.** Whether to schedule the convergence ladder; whether to request the
Hyfydy key after the deadline.

## 4. Optuna + SCONE review — read-only; ensemble scaffold delivered

**Measured (on a %TEMP% copy of optuna_walk.db; 19 studies — v8b/v9/v10 are NOT in this
db, and the raw-sqlite contiguity check on the copy (study_ids 1–19, v6 flowing directly
into curr_s1) shows they were never in the file at all — not lost in copying; per
AGENTS.md that work ran on another machine, so the main db shares the same gap).**
22–70% of every stage-3 study's budget sat on the −320 sentinel (s3c 70%, s3b 68%);
the space grew to 32 params at 1.3–2.5 trials/param. fANOVA (curr_s3k): `desc_e` .187,
`ib_rge` .153, `ib_e_central` .131, `desc_f` .101, `contact_onset` .049 — 27 of 32 params
≤ .035. Median eval cost 1.0–1.1 min/trial. SCONE practice we lack: robustness inside the
objective (NoiseController/PerturbationController/sampled initial states), covariance
adaptation, multi-start pools, tiered thresholds.

**Answer to your ensemble question.** **Yes** for seeds + small plant jitter (±5–10%
Fmax/damping; Koelewijn 2022 shows training under uncertainty moves the optimum). **No**
for synthetic IK perturbation: the repo holds exactly ONE reference cycle per leg
(walk1_ik.mot spans 0.5–2.5 s; `subject01_walk1.mot` is all zeros — a placeholder), so
jittering it would train robustness to a biased target. Get more reference data first
(goal 5).

**Delivered.** `ensemble_objective_scaffold.py` (score = mean − 0.25·std over nominal +
plant-jitter + reference/pose variants; 2-stage finalist rescoring; ~3 min/trial) —
smoke-tested (−169.543, finalists [0,1,3]) but **not wired in**; plant jitter needs a
runner-side Fmax hook that does not exist yet.

**For Ben to decide.** (a) Adopt the ensemble (N=3) for the next study? (b) Switch to
`TPESampler(multivariate=True, constraints_func=feasibility)` + a ~10-param space + two
enqueued seeds? (c) Approve the default-off Fmax jitter hook?

## 5. Gait-library catalog — six SimTK projects, downloads now login-gated

**What was done.** All six projects catalogued (contents, formats, licenses, sizes) from
their home pages and download listings, with per-project integration tasks grounded in
`kine_ref.py` / `bsolve_ik.py` / `fit_synapses.py` at line citations. **Key finding:**
SimTK file downloads now require a (free) login — verified on 3 of 6 projects via
login.php redirect stubs; inferred on the other three (simtk.org rate-limited the client
before those probes). Nothing was downloaded; exact zip contents are labeled expectations
pending your manual download.

**Priority shortlist.** (1) `muscfib_walkrun` Walk zips — speed-binned walking references
for the duty/cadence gaps + F-L-V validation of the converted muscles (**CC BY-NC 3.0 —
flag before any commercial reuse**; pilot = README + Subject 01/02 Walk ≈ 132 MB).
(2) `predictivesim` — 17 MB for a reflex+contact-state walking-controller reference
implementation + a second experimental dataset. Then nmbl_running EMG (211 MB selective),
assistloadwalk (load-perturbation eval), runningsim (moment-arm cross-check); crouchgait
deferred per your note.

**For Ben to decide.** Log in once and pull the download table (goal5 §"manual download
list"); pick a landing folder OUTSIDE the repo — only derived tables get committed.

## 6. Circuit drafts — Rybak / Shevtsova / Shinohara (replication/)

**What was produced.** First-pass drafts + rules JSONs for Rybak 2006a/2015/2024,
Shevtsova 2026 (eLife RP107480), and Shinohara 2025 (bioRxiv 687930; both full texts
local) — every edge cited to a paper table/figure with an EXISTS-vs-NEW flag
cross-checked against the compiled net (`_net_edges.py`: 82-edge default + corrected
170-edge all-conditionals probe).

**Both supervisor blocks resolved.** (1) IaIN↔IaIN mutual inhibition reclassified
EXISTS-conditional (g 0.5 under `ia_in`+`full_rules`, build_network.py:795-800; Rybak
2006a Table 2 −0.1 is a second source — merge provenance, do not duplicate); the dump also
shows IBIN↔IBIN g 0.5, so both ask-named gaps already exist. (2) The misfetched
`lit_rybak2006b_fulltext.txt` (an unrelated HIF-1α paper from a guessed PMCID) is deleted
with a warning in rybak_draft.md — no Rybak-2006b text exists locally (paywalled).

**Genuinely-NEW edge classes for you:** RG-E↔RG-F/self recurrent excitation + opposite-RG→PF
gating + MLR→PF (audit row #23's literature basis was incomplete); V3-F crossed flexor
excitation + the V2a/V0V chain (Shevtsova Table 1); per-muscle flexor-length→RG-F/IN-F/PF-F
fan-in with the 0.9·lmax + v^0.6 encoding (Shinohara eqs 10–11 — the most concrete spec for
our duty/frozen-leg gap). Note the lineage allocates the rhythm to the FLEXOR center; our
DRIVE biases E.

**For Ben to decide.** rybak §8's five items (provenance merge, trial recurrent excitation
at 0.0125-scale default-0, opposite-RG→PF gate, DRIVE→PF reinstatement, crossed F↔F route),
plus biped-reduction choices (lumbar-only Shevtsova; cat 7-muscle → 92-muscle grouping for
Shinohara). Ben edits before any build.

## 7. Six-synergy neuromechanical model — basis found, replay + demo run

**What was produced.** `synergy_model.py` + `synergy_basis.npz` + demo npz. The six
synergies of record were **found, not recomputed**: the per-leg rank-6 NMF basis from
`fsa_backsolve.py` on `bsolve_ik.py`'s back-solved activations, exported verbatim
(W [43×6], H [6×121] per leg; 86/92 actuators — the 6 outside are exactly the prunes).

**Key numbers.** Six-channel replay (fixed W, per-timestep NNLS) reconstructs the
reference activations at pooled VAF **0.945 / R² 0.920** (0.900/0.860 vs raw acts); the
NNLS coefficients coincide with the NMF H to ~1e-6 of peak — the converged NMF was already
per-timestep optimal. Weakest groups: hip_abd/hip_add/hip_ext. The suspended no-ground
open-loop demo ran a clean 6.00 s (3000 steps, leg damping 0.5, knees −109..+11°);
**undamped it fails** at t=1.786 s (rect_fem_l patella-follower singularity, DOF 66).
Open-loop tracking is honestly poor (RMSE 28–58°, negative r) — the activations were
solved for loaded stance; the neuromechanical loop is doing real work this pass omits.
`runner.py`/`params.py`/`build_network.py`/`_curriculum.py` untouched.

**For Ben to decide.** (a) 1:1 synergy↔PF-population mapping vs shared half-centers;
(b) S5/S6 are not bilaterally phase-stable — keep or drop per side; (c) phase-windowed H
vs measured-H replay; (d) retune on longer records before porting into SNS topology
(wiring would go through your connectome-spec route).

---

## Decisions for Ben

1. **Goal 1:** confirm the feature belongs in the block editor; optional 2-line
   localStorage autosave fix.
2. **Goal 2:** sign off (or amend) the six VEST/balance items in goal2 §6; when to launch
   `_curriculum.py 4 30` (NOT launched).
3. **Goal 3:** approve/schedule the offline 1/0.5 ms convergence ladder; decide on the
   post-deadline Hyfydy key. No model change before the deadline without accepting the
   void-everything cost.
4. **Goal 4:** adopt ensemble N=3 + multivariate TPE with `constraints_func` + ~10-param
   space for the next study? Approve building the default-off runner Fmax-jitter hook
   (it does not exist yet)?
5. **Goal 5:** manual SimTK download of the shortlist (login required); pick a landing
   folder outside the repo; note the muscfib CC BY-NC 3.0 license.
6. **Goal 6:** the five Rybak §8 rulings + biped-reduction choices; provenance merge for
   IaIN↔IaIN (no duplicate edge); whether to reinstate DRIVE→PF (audit #23 amendment).
7. **Goal 7:** the four port decisions above (mapping, S5/S6, H form, retune-first).

## Verified vs not verified

**Verified by me this session (the checks I ran):** all seven goal reports + the three
replication drafts + README + `ensemble_objective_scaffold.py` read in full;
`goal4_db_analysis.json` structure confirmed (per-study n/best/fANOVA/sentinel fields);
file existence + mtimes for every deliverable (`dir`); the goal-2 abs() fix on disk
(`findstr` → `_curriculum.py:279`); v2.2 markers in `connectome_block_editor.html`;
`connectome_editor.html` (09/22 12:38 AM) and `CONNECTOME.md` (09/22 09:18 AM) unchanged;
`synergy_basis.npz`, `synergy_replay_air_20260923_130556.npz`, `goal2_smoke_balance.npz`,
`lit_rybak2006_fulltext.txt`, and the three replication rules JSONs all present.

**Carried from the goal agents' reports (commands + outputs recorded there, not re-run
here):** every behavioral/numerical check — the 37/37 editor suite and node --check; goal
2's py_compile, defaults-off diffs, topology fingerprints, 2 s smoke; goal 3's MuJoCo
inspections + SEE toy; goal 4's db analysis + scaffold smoke; goal 5's curl/login probes;
goal 6's `_net_edges.py` dumps; goal 7's replay/demo runs.

**Standing gaps (stated in the reports, none filled):** goal 1 not run in a real browser;
goal 2's stage-4 optimization never exercised (no study, by instruction; 8 s eval not run)
and the vestibulospinal-extensor citation is textbook-level, not full-text verified; goal 3
ran no Hyfydy and no convergence ladder; goal 4 analyzed a db that never contained
v8b/v9/v10 (that work ran on another machine) and executed no ensemble study; goal 5 downloaded nothing (login gate) and zip contents are
expectations; goal 6 built/ran no simulations and four paper details stay unverified
(Rybak 2006b text, 2006a appendix params, Shevtsova kτ, Shinohara cost terms); goal 7 made
no ground-contact attempt and no SNS wiring.

*Sources: the deliverable files read in full 2026-09-23, plus the on-disk spot-checks listed above.*

## Post-run addenda (2026-09-23, same session)

- Goal 3: the recommended convergence ladder was EXECUTED
  (spinal/goal3_timestep_ladder.py -> goal3_convergence_ladder.md): same
  20-muscle recorded control replayed into the harnessed plant at 2/1/0.5
  ms. Ground: 2 ms vs 0.5 ms = 0.44 deg RMS (max 11.9 deg, knee/subtalar
  transients); air: 0.19 deg RMS; error roughly halves with the step -
  first-order converged, so the 2 ms contract is now quantified rather
  than asserted. The ladder costs 1-4 s per variant (cheap enough to run
  before any recorded-winner claim). Caveat: replay drives the 20 key
  muscles from a recorded run, not all 92 - a representative loading
  scenario, not a bit-exact rerun.
- Goal 5: Falisse 2022 predictsim_mtp staged openly at
  D:\temp\gait_lib_staging (46 MB) as a substitute candidate for the
  login-gated predictivesim; generalized loader gait_lib_loader.py built
  and regression-verified vs kine_ref (< 1e-9); Falisse reference cycle
  extracted; our six synergies do NOT transfer to it (VAF ~ 0, R2-vs-0
  0.60) - argues for multi-reference training. muscfib_walkrun still
  needs Ben's SimTK login.

## Second post-run addendum (2026-09-23, goals 3 and 5 pushed further)

- Goal 3, IMPLEMENTED at wrapper level (goal3_fidelity_variants.py -> goal3_fidelity_variants.md):
  (A) ERROR-CONTROLLED INTEGRATION now exists - a step-doubling error
  estimator inside the 2 ms lockstep (coarse 2 ms vs two 1 ms probes;
  refine to 4 x 0.5 ms when local error > tol). Measured: only 5.0% of
  samples exceed tol=0.01 deg; tol behaves as designed (loose tol
  degenerates to the 1 ms probe path, RMS 0.251 deg = the ladder's 1 ms
  row). Wall 3.2-3.6 s vs 1.7 s (fixed 2 ms) / 4.4 s (fixed 0.5 ms).
  Honest limit: local error is controlled; the chaotic global divergence
  (phase) is not, and cannot be, by any local controller.
  (B) CONTACT DAMPING variants demonstrated stably in the impedance
  framework: tc 0.01/zeta 0.7 triples peak impact force (3237 N vs
  1001 N baseline); solimp mid 0.90/width 0.003 softens it (882 N). The
  naive direct-(k,b) negative-solref route was tried and is UNSTABLE at
  2 ms (NaN by t=0.012 s) - concrete confirmation of the report's
  'partially native' verdict. Nothing production-facing changed.
- Goal 5, FIRST USE of the staged Falisse data (gait_lib_score.py ->
  goal5_falisse_scoring.md): the s3k winner scores kine -189.5 against
  subject01 vs -200.5 against Falisse's independent predicted gait -
  generalizes within ~11 points on a ~190 scale. Multi-reference scoring
  is now a one-liner against any loader reference.
- Still Ben's: SimTK login downloads (muscfib Walk zips + predictivesim)
  and the goal-3 production ruling (keep 2.3.7/implicitfast/2 ms through
  the dissertation; SCONE+Hyfydy for balance after).

## Third post-run addendum (2026-09-23: first training-side use of library data)

- Reference-tension pilot (gait_lib_pilot.py v1 -> gait_lib_pilot2.py ->
  goal5_pilot2.md): the PRODUCTION s3k winner was run through the exact
  stage-3 eval machinery against BOTH references (REF_CACHE swap, no core
  edits). Result: the s3k gait matches Falisse's predicted walking BETTER
  than subject01 at every drive tested (kine -228.0/-222.5/-224.9 vs
  -241.6/-240.9/-239.4 at drives 1.77/2.08/2.39) - the subject01 misfit
  concentrates in the known knee-flexion/hip-range gaps. Drive response is
  flat at +-15%; the tension lives in which reference is fit, not the
  drive derivative.
- REGRESSION FLAG for Ben: the v10-era winner `runner --fitted --best10`
  NO LONGER WALKS under current physics (pilot v1: 100% double-support,
  zero cycles at any drive 2.6-3.7; npz diagnostics in %TEMP%). All
  v8b/v9/v10 recorded scores predate the RoM-limit/contact-surgery era -
  treat them as historical until re-verified. Only s3* winners are
  current.
- Still Ben's: SimTK login downloads (muscfib Walk zips + predictivesim)
  and the goal-3 production ruling.

## Fourth post-run addendum (2026-09-23: goal-3 features IN the production runner; goal-5 intake ready)

- runner.py now carries the two goal-3 features as DEFAULT-OFF flags (the
  goal-2 conditional-topology pattern): `--adaptive-tol X` (error-controlled
  integration - any 2 ms step whose step-doubling error estimate exceeds X
  deg is refined to 4 x 0.5 ms; 0 = the exact previous single-step path) and
  `--contact-damp lessviscous|nonlinear` (the two stable impedance-framework
  contact variants; off = production solref/solimp untouched). Verified:
  py_compile exit 0; three 20 s smokes all exit 0, stayed up, no NaN
  (defaults / --adaptive-tol 0.05 / --contact-damp nonlinear - the contact
  variant visibly changes behavior: COM ends 0.85 m vs 0.80, subtalar range
  25 deg vs 35); AARL_NPZ honored (smoke npz in %TEMP%, spinal_run.npz
  mtime untouched). Tendon elasticity remains the one Ben-gated piece.
- NEW gait_lib_intake.py: drop Ben's SimTK zips into
  D:\temp\gait_lib_staging\downloads\ and run one command (myo env) -
  it extracts, classifies every .mot by its columns, pairs IK+GRF, builds
  kine_ref-schema references via gait_lib_loader, saves them to
  spinal/gait_refs/*.npz, and reports unpaired files. Verified on the
  staged Falisse data (1 reference saved: Case_40_motion.npz). muscfib
  fiber-length files will list UNPAIRED - they feed the F-L-V validation
  route, not the reference route.
- Ben's two items unchanged: SimTK downloads (now a one-command intake) and
  the goal-3 production ruling (SEE conversion per goal3_mujoco_fidelity.md
  3.2 if approved).
