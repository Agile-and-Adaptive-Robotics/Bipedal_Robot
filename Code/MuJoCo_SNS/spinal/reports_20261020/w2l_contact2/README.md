# w2l_contact2 — Ben's corrected W2L connectome in MuJoCo + SNS (2026-10-20)

Ben's ask (2026-10-20): "Now that the w2l model has been corrected in
animatlab, I want you to run the updated connectome in
Mujoco_SNS_toolbox and tune it. /supervisor-gate the results. Reflect
these contact pathways in the other Mujoco models (syn6, s3k, ...) as
well."

## Source of truth

`Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG.aproj`
(Ben, 2026-10-18 21:56, commits 2955e4e9 "physical contact sensing closed
(Li-style chain)" + 404bd246 "model updates"). Mined by
`spinal\w2l_cpg\make_bilateralrg2.py` (GUID-level; mine record:
`spinal\w2l_aproj_mine_20261018.json`) into the `bilateralrg2` template
key of `spinal\connectome_templates.json`. Census **96 neurons / 204
synapses / 10 base inputs / 12 muscle outputs**. The 2026-10-03 canonical
export (`ben_connectome_canonical.json`, 110n/216e) PREDATES this aproj
and is superseded for runnable purposes.

## The wiring delta vs the runnable baseline (bilateralrg, 2026-09-26)

1. **Toe cascade REMOVED** (10 edges): toe no longer drives ipsi RG ext /
   Hip-Knee PF ext / Hip-Knee MN ext.
2. **Heel cascade re-typed WEAK** (10 edges): dedicated types
   `Heel_Ipsi_*_Excite` gmax 1e-8 S = the build_aff.pl "Afferent HC
   Excite SynAmp 0.01" class → tag `contact_hc` (gain family
   aff_HC_excite).
3. **NEW crossed flexor edge** (2): `Heel_Contra_RGflx_Excite` — heel →
   CONTRAlateral RG flx, weak class → tag `contact_contra_flex`. THE
   antiphase mechanism DESIGN.md 2026-09-22 named as missing.
4. **NEW ankle contact edges** (4): `Heel_Ipsi_AnkMNext_Excite` +
   `Toe_Ipsi_AnkMNext_Excite` at GENERIC 5e-7 S (the only strong contact
   edges left) → tag `contact_ank`.
5. **Li-style per-joint afferent chain confirmed** (Ia/II flexor → PF flx
   + RG flx; Ib extensor → PF ext + RG ext; all dedicated weak types).
6. **Knee-extensor Renshaw LATERALIZED**: the aproj's two
   "L Knee MN ext RE" nodes (GUIDs a3828afc / efa2eb22) traced at GUID
   level — a3828afc serves the R side (R Knee MN ext → a3828afc → R Knee
   MN ext / R Knee ext Ia / R Knee MN flx RE). It is Ben's R-knee-ext RE
   with a stale L label (his 2026-10-03 export renamed it RC_8971);
   relabeled **"R Knee MN ext RE"** here (LABEL fix only — every edge
   keeps its GUID-resolved endpoint).
7. Commissural labels normalized to the aproj's ("L c1"/"L V3" instead
   of the old template's "c1_L"/"V3_L") — same edges.

### FLAGGED for Ben (not blockers, transcribed with his 10-03 export's ruling)

- **F1**: the edge "R Hip MN flx RE → R Hip MN ext RE" carries the
  V3-Commissural-EXCITE synapse type (GUI clone-type artifact family,
  documented in w2l_cpg README). Ben's 2026-10-03 export rules it
  `inh/rc_mutual`; transcribed rc_mutual INH. If he INTENDED excitatory,
  flip the sign in the template edge.
- **F2**: effective-strength doctrine = the synapse TYPE's conductance
  (the .pl build chain's mechanism; per-instance values are GUI prefill).
  The new dedicated types put the heel cascade + crossed edge at the weak
  0.01 µS class — if Ben meant these at generic strength, change
  `contact_hc`/`contact_contra_flex` starting gains (they are tuning
  knobs regardless).

## Where things run

- Pure-net smoke: `spinal\w2l_cpg\smoke_w2l.py --template bilateralrg2`
  → **PASS** period 2.997 s, antiphase r −0.782 (old template regression
  re-run: PASS, bit-identical 1.000 s / −0.525 / 22-22).
- MuJoCo ground walk: `spinal\w2l_mujoco\test_w2l_ground.py
  --net=bilateralrg2 ...` (default `--net=aff` unchanged = the M5
  lineage of record).
- Tuning: this folder, `_w2l_contact2_campaign.py` (study
  `curr_w2l_contact2_r100` + wean ladder in `optuna_w2l_contact2.db`;
  metric + survival gate VERBATIM from reports_20261016\w2l_wean).
- **First datum (fresh process, Ben's literal untuned aproj defaults:
  te 3/tf 4, ch 0.08/1.0/0.08, ctn 4/50, ibnA 1, comm 1): SURVIVES —
  obj 90.0, 19/17 cycles, inter-leg phase 0.498, fwd +0.545 m, tilt
  15.8°.** The aff lineage's untuned defaults at weaned drive were DEAD
  (obj 7).

## Reflected into the other families (default-0, bit-identical off)

`params.G` keys `contra_flex` (heel IN side A → RG-F side B, EXC) and
`contact_ank_mn` (heel+toe INs → that side's ankle_pf MN pools, EXC) in
build_network.py (s3k lineage), build_network_spiking.py,
build_network_w2lvar.py, build_network_syn6.py. Gates 2026-10-20:
defaults census (410,376,1186) + s3k anchor −161.56754173676563
BIT-EXACT; w2lvar (888,382,7866); syn6 (794,382,3033); ON-probe delta
exactly +2+2×18=38 synapses in w2lvar/syn6 (standard/spiking deltas
include the pre-existing stance_fb machinery their heel_rge=0.5 also
activates); one-step finite everywhere
(`reports_20261020\w2l_contact2\probe_contact_knobs.py`).

## Basin + robustness outcomes (2026-10-20, post-M1-gate)

STRICT basin gate (±1% joint perturbation of all 9 knobs, 12
fresh-process trials, pass bar 11/12): **NO configuration passes; NO
basin pass is claimed.**

| config | basin |
|---|---|
| tuned winner te=5.27 (round_best_r030) | 7/12 |
| r100 winner te=4.35 | 1/12 (one-leg marginal) |
| Ben's literal defaults | 6/12 (cycles robust 26-35 both legs; PHASE drifts out of [0.35,0.65]) |
| defaults + comm ×2 | 4/12 (over-coupling drags phase toward sync) |
| defaults + comm ×4 | 1/12 (phase 0.076 ≈ synchronized) |
| defaults + ch_contra 0.5 | **9/12 (best)** — the crossed heel→contra-RG-flx edge RAISED from Ben's 0.08 hardens the phase lock (7→9); remaining failures are phase excursions |

The robust-objective study (`curr_w2l_contact2_robust`, 34 trials,
objective = mean over base + 2 perturbations) converged back to Ben's
literal defaults as the most robust point (obj_mean 84.03,
`robust_best.json`). Failure mode across configs: stepping is robust,
the inter-leg phase lock is the fragile piece. Only the LAST basin run
persists in `basin_contact2.json` (the ch_contra-0.5 one, survived=9);
all six counts are reconstructable from `verify_trials.csv` rows
10-93 (12-row blocks).

## Wean ladder — DEGENERATE (design flaw; no weaning claim)

`dscale` was BOOKKEEPING ONLY in this campaign's ladder: te/tf were
free search dims at every rung (`_w2l_contact2_campaign.py knobs_of`
uses p["te"] directly), so the r085–r030 "rungs" re-ran effectively the
same optimization and returned the same winner at the same effective
drive (the six ladder verifies of 04:27–04:40 are byte-identical,
obj 87.17257 each). **No weaning result is claimed for this campaign**;
a real wean would scale the winner's te/tf per rung and re-tune
feedback only (the reports_20261016\w2l_wean design).

## CSV note

`_robust_round.py` inherited `_w2l_wean_eval`'s default CSV path, so
its 100 `rob`-tagged rows (2026-10-20 11:01–11:11) appended to the OLD
campaign's `reports_20261016\w2l_wean\wean_trials.csv` — append-only,
old rows untouched. The main campaign's trials are in
`contact2_trials.csv`; winner/fresh verifications in
`verify_trials.csv`.

### ch_contra bracket (2026-10-20, supervisor pass 2)

Defaults with ch_contra swept in isolation: 0.08 → 6/12, 0.3 → 6/12,
0.5 → 9/12, 0.8 → 9/12 (the 1.10 datum at 7/12 had all other knobs
TPE-changed, not comparable in isolation). READ: the crossed-flexor
edge buys a robust +3/12 plateau at 0.5–0.8 and saturates there — the
remaining failures are phase excursions with cycles intact, so the
residual marginality is structural (soft inter-leg coupling through
the yield/rig harness), not a gain problem. Decision whether to accept
9/12 or restructure is Ben's (DESIGN.md open issue #4).
