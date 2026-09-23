# Goal 5 — SimTK gait-library catalog + integration plan for the spinal-circuit pipeline

**Date:** 2026-09-23 · **Scout:** gait-library scout (GLM subagent)
**For:** Ben (PI). **Deliverable of GOAL 5** — analysis only, no code changed.

## Read this first — what you must decide

1. **All SimTK file downloads now require a logged-in (free) SimTK account.** Verified on 3 of
   the 6 projects (runningsim, muscfib_walkrun, predictivesim): every file's
   `download_confirm.php` page redirects to `/account/login.php?triggered=1`. The other three
   use the identical download mechanism (same URL pattern on the same FusionForge install), so
   the same gate is inferred there (not directly re-tested — see *Verified / not verified*).
   Project pages and download listings stay anonymously readable. **Action:** log in once at
   simtk.org, then the "Ben downloads manually" list below is a few clicks per project.
2. **Recommended first download (~460 MB total):** `predictivesim` (17 MB) + `muscfib_walkrun`
   README + 2 Walk zips (~132 MB, scale up to ~340 MB for all five subjects). These two attack
   three of our four named gaps (duty, cadence targets, muscle F-L-V fidelity) and give a
   second walking dataset + a reflex-controller reference implementation. Details below.
3. **License watch:** `muscfib_walkrun` is **CC BY-NC 3.0** (non-commercial). Fine for the
   dissertation and lab research; flag it before any commercial/patent-adjacent reuse.
   Everything else catalogued is CC BY (3.0/4.0) or BSD-style.

---

## Our pipeline hooks this plan is grounded in

| Hook | Where | What it does today |
|---|---|---|
| Reference gait cycles | `Code/MuJoCo_SNS/spinal/kine_ref.py:40-43` | Hardwired to ONE trial: `subject01_walk1_ik.mot` + `subject01_walk1_grf.mot` |
| Reference loader | `kine_ref.py:72-124` (`load_reference`) | Per-side mean cycles phased by GRF loading onsets; duty, period, interleg lag, double-support |
| Contact-free cycle fallback | `kine_ref.py:166-172` | RG-E burst onsets when contact channels absent — reusable for datasets without GRF |
| Frozen-leg guard | `kine_ref.py:297-312` | Bilateral bonus / frozen-leg cost (the s3k fix) |
| IK/SO back-solve + muscle validation | `Code/MuJoCo_SNS/spinal/bsolve_ik.py:1-33` | OpenSim IK → MuJoCo muscle-length validation + NNLS activations vs OpenSim SO; `opensim-cmd` at `bsolve_ik.py:57` |
| Muscle↔group mapping | `bsolve_ik.py:46` (`muscle_map.classify`, `GROUPS`) | Cross-model muscle-name reconciliation point |
| NMF synergy scaffold | `Code/MuJoCo_SNS/spinal/fit_synapses.py:10-13, 33-60` | k=5 NMF (Ivanenko 2004; Kibushi 2018 speed-invariance claim) + per-phase NNLS → W_PF_MN; CSV in (t + muscle columns) |
| Known gaps (AGENTS.md, Sept 2026) | duty (E-duty 0.16–0.34 vs ref 0.61), cadence (1.18–1.4 Hz vs 0.81), hip phase inverted vs RG anchor, frozen-left-leg history; MuJoCo muscles = rigid-tendon simplified Thelen (no series elasticity) — open fidelity question | — |

Current reference stats (single trial): 1.23 s / 0.81 Hz, duty 0.61, knee −69.7° (AGENTS.md).
**Every dataset below's main job is to replace "one subject, one speed" with "subjects × speeds"
references, validation data, and perturbation scenarios.** None of them fixes our architecture
gaps by themselves (hip phase / E-duty ceiling are ours to build) — they fix *targets and
evaluation*.

---

## 1. muscfib_walkrun — Arnold et al. 2013 — **relevance: NOW (top priority)**

**Page:** simtk.org/projects/muscfib_walkrun (fetched 2026-09-23). Group 769.

**What it is:** "Simulated muscle fiber lengths and velocities during walking and running" —
simulation results behind Arnold, Hamner, Seth, Millard & Delp, *J Exp Biol* 216:2150-2160
(2013), doi:10.1242/jeb.075697. Five subjects, **walking 1.0–1.75 m/s (4 speeds)** and running
2.0–5.0 m/s (4 speeds), EMG of 11 lower-limb muscles; muscle-tendon dynamics simulations
yielding fiber lengths, fiber velocities, forces, and force-generation ability per speed
(abstract: soleus force-generation ability falls with walking speed; walk→run transition
raises it by cutting fiber velocities).

**Contents (from the download listing, sizes as listed):**
- `README.txt` (6 KB) — "Details about using these simulation results and models"
- `Unscaled Model.osim` (398 KB) — "based on the model described by Arnold et al. 2010 but
  uses the **Schutte1993Muscle** model formulation"
- Per subject (01, 02, 04, 18, 19): `Subject NN Models.zip` (~44–57 KB, scaled .osim + README),
  `Subject NN Walk.zip` (60/72/62/79/67 MB), `Subject NN Run.zip` (98/114/143/136/122 MB)
- Totals: **≈ 953 MB everything; ≈ 340 MB for the five Walk zips only.**

**Formats:** expected OpenSim workflow files (.osim models, .sto/.mot results) — the listing
calls them "simulation results" without enumerating; exact per-file inventory is login-gated,
so treat contents as *expected, verify against README.txt after download*. Whether GRF
channels ship in the Walk zips is unknown; if absent, phase cycles kinematically and reuse the
contact-free fallback pattern at `kine_ref.py:166-172`.

**License:** CC BY-NC 3.0 (NonCommercial) — Stanford 2013 (license text on the downloads
page). **Size:** ~953 MB total / ~340 MB walk-only. **Anonymous download:** NO (login gate
verified on this project's `README.txt` confirm link).

**Why it matters to us (maps to 3 of 4 gaps):**
1. It is the only catalogued project with a **walking speed series** — the raw material for
   speed-binned reference cycles → real human duty-vs-speed and cadence-vs-speed targets,
   replacing the single 0.81 Hz point our DRIVE→RG-period calibration currently aims at.
2. It is the reference for the **muscle F-L-V fidelity question** AGENTS.md flags on our
   converted MuJoCo muscles (rigid-tendon simplified Thelen, no series elasticity).
3. Its per-speed EMG-driven activations are exactly the input class `fit_synapses.py` was
   built for, letting us *test* the Ivanenko speed-invariance / timing-shift claim on human
   data before leaning on it in the dissertation.

**Integration tasks (concrete):**
1. Ben downloads `README.txt` first; verify inventory; then `Subject 01 Walk.zip` +
   `Subject 02 Walk.zip` (~132 MB) as the pilot.
2. **Extend `kine_ref.load_reference`** (`kine_ref.py:72`) to accept indexed sources and
   return a per-speed reference family {1.0, 1.25, 1.5, 1.75 m/s}: per side — T, duty, interleg
   lag, ds fraction, mean/range/knee_min (same keys `kine_ref.py:89-124` already produces);
   if GRF is absent, phase by kinematic heel-strike events (fallback pattern
   `kine_ref.py:166-172`).
3. **Cadence/DRIVE calibration table:** fit human cadence-vs-speed (and duty-vs-speed) from
   that family; use it as the target curve for DRIVE sweeps in `optuna_walk` studies (today's
   cadence targets are one point, 0.81 Hz).
4. **F-L-V validation:** run our converted model along their walking kinematics via the
   `bsolve_ik.py` comparison machinery (`bsolve_ik.py:25-31` outputs; current MuJoCo↔OpenSim
   length agreement median r 0.89 per AGENTS.md) and compare our muscle force-generation
   ability envelopes against their published per-speed fiber length/velocity results;
   account for the **Schutte1993Muscle vs Thelen** formulation difference (this comparison
   quantifies whether the rigid-tendon Thelen approximation needs the tendon-stiffness /
   Millard-plugin upgrade AGENTS.md lists as open).
5. **Synergy check:** convert their per-speed activations to the `fit_synapses.py` CSV format
   (t + muscle columns, `fit_synapses.py:33-44`; names reconciled via `muscle_map`) and run
   NMF k=5 per speed — evidence for/against speed-invariant synergy structure +
   timing/amplitude shift (`fit_synapses.py:10-13`) from *human walking* data.

**Papers (cited on the project page):** Arnold EM, Hamner SR, Seth A, Millard M, Delp SL
(2013) *J Exp Biol* 216:2150-2160.

## 2. predictivesim — Dorn et al. 2015 — **relevance: NOW (cheapest, most architectural)**

**Page:** simtk.org/projects/predictivesim (fetched 2026-09-23). Group 728.

**What it is:** "Predictive Simulation of Loaded and Inclined Walking" — Dorn, Wang, Hicks &
Delp, *PLoS ONE* 10(4):e0121407 (2015). Predictive (non-tracking) muscle-driven walking from
energy minimization: "includes controllers based on **muscle force and stretch reflexes and
contact state of the legs**" (abstract). Validation: 92% of joint-angle and 78% of
joint-torque trajectories within 1 SD of experimental data; predicts adaptations to carried
loads and inclines.

**Contents (from the download listing):**
- `ExperimentalData.zip` (17 MB) — "Raw C3D files for each subject (Static standing trials and
  dynamic walking trials)"
- `pdsim-source.zip` (108 KB) — "Source code for control simulation and optimization. Tested
  with OpenSim 3.2 and Simbody 3.3.1 on OSX 10.10.3"

**Formats:** C3D raw mocap; C++ source (old OpenSim 3.2-era API — porting to 4.x is *my
inference*, not verified). **License:** Stanford 2015 BSD-style + "All data may be used
freely for academic research purposes". **Size:** ~17 MB. **Anonymous download:** NO (login
gate verified on this project's confirm link).

**Why it matters to us:** this is the closest published cousin to our program — a walking
controller built from reflexes + contact state, not tracked data, evaluated against
experimental envelopes. It is a reference implementation for exactly the sensory-phase /
contact-reset machinery our RG lacks (AGENTS.md diagnosis: E-duty ceiling needs sensory phase
reset). Cross-reference: our local SCONE 2.4.4 install (scone skill, added 2026-09-23) ships
runnable Geyer-Herr reflex controllers of the same lineage — a live comparison target without
porting anything.

**Integration tasks (concrete):**
1. Ben downloads both files (tiny).
2. **Read `pdsim-source` controller structure:** reflex gain organization, force-feedback vs
   stretch-feedback pathways, contact-state phase switching. Map each onto our
   equivalents — heel/toe mechanosensors → RG-E (heel_rge/toe_rge), stance-gated Ib (ib_rge),
   PRESET transient reset, KINH swing gating — and list anything theirs has that ours lacks,
   as a supervisor-reviewable table (paper + code line on their side; our `params.G` knob on
   ours; gaps labeled as gaps).
3. **Second walking dataset for `kine_ref`:** run the C3D trials through the OpenSim
   Scale/IK workflow we already invoke in `bsolve_ik.py` (opensim-cmd, `bsolve_ik.py:57`),
   producing cross-subject reference cycles — move our acceptance framing from "one subject's
   cycle" to "experimental mean ± SD envelope" as in their 92%-within-1-SD methodology.
4. Their **loaded/inclined predictive gaits** become ready-made robustness references when we
   add perturbation scenarios (pairs with assistloadwalk, below).

**Papers:** Dorn TW, Wang JM, Hicks JL, Delp SL (2015) *PLoS ONE* 10(4):e0121407.

## 3. nmbl_running — Hamner & Delp 2013 — **relevance: LATER**

**Page:** simtk.org/projects/nmbl_running (fetched 2026-09-23). Group 603.

**What it is:** 10 male experienced runners (≥30 mi/week, Stanford Human Performance Lab
instrumented treadmill) at **2/3/4/5 m/s**: motion capture, EMG, GRFs, subject-specific
scaled models, and full muscle-driven simulations. Workflow shown on page: Scale → IK → RRA →
CMC → IAA; "Simulations were generated using OpenSim 2.4". Paper: Hamner & Delp, *J Biomech*
46(4):780-787 (2013) — soleus dominates upward + forward acceleration at all speeds.

**Contents (from the download listing):** `RAW_EMG_DATA.zip` (8 MB, "raw (i.e., unprocessed)
EMG recordings") + `subject01…subject20.zip`, ten files of **192–224 MB each** (203+224+198+
212+192+205+193+224+206+194 MB) ≈ **2.06 GB total**. Per-subject zips carry mocap/EMG/GRF +
scaled model + simulation results (page description; exact inventory login-gated — verify on
download).

**License:** Stanford 2013/2014 BSD-style (per-package license popovers). **Anonymous
download:** NO — same mechanism; the gate on this specific project is *inferred* (see
*Verified / not verified*).

**Why it matters to us / tasks (when we get to it):**
1. Only catalogued dataset with a **4-speed series of raw EMG** → process to envelopes and
   use as human afferent-side reference data for our Ia/II/Ib conversion checks (the
   (L−Lmid)/Lhalf, L̇/0.6, F/Fmax → current maps in our pipeline; AGENTS.md conversion-note).
   Standard filtering/rectify/envelope preprocessing is *my inference* of what's needed.
2. Extends the cadence-vs-speed calibration into the running regime if the lab ever targets
   faster gaits (RG speed modulation beyond walking).
3. 10-subject spread → reference variance bands for `kine_ref` (cross-subject envelopes).
4. Selective download advice: `RAW_EMG_DATA.zip` + one subject zip (~211 MB) covers items 1–2
   without the full 2 GB.

**Papers:** Hamner SR, Delp SL (2013) *J Biomech* 46(4):780-787.

## 4. assistloadwalk — Dembia et al. 2017 — **relevance: LATER**

**Page:** simtk.org/projects/assistloadwalk (fetched 2026-09-23). Group 1082.

**What it is:** 7 male subjects × 4 conditions — (a) free-speed unloaded, (b) 80% of free
speed, (c) **+38 kg torso load** at a new free speed, (d) +38 kg at speed (a) — plus OpenSim
simulations of **7 ideal massless single-DOF assistive devices** (hip abd/flex/ext, knee
flex/ext, ankle PF/DF) and their effect on metabolic cost and muscle activity. Paper: Dembia,
Silder, Uchida, Hicks, Delp, *PLoS ONE* 12(7):e0180320 (2017).

**Contents (from the download listing):** `…simulations_of_experiments.zip` (646 MB, DOI
10.18735/S5S69G), `…assistance_subjects_05_07_09.zip` (528 MB, DOI 10.18735/S51Q3J),
`…assistance_subjects_10_11_12_14.zip` (674 MB, DOI 10.18735/S5X12B), `assistloadwalk_scripts.zip`
(131 MB). Total **≈ 2.0 GB**. Formats: OpenSim simulation + experimental data zips; scripts =
source code (listing marks it "Source code"; language not stated on page — verify on
download).

**License:** BSD-style, "Copyright (c) 2016, Christopher Lee Dembia". **Anonymous download:**
NO (inferred — same mechanism; see *Verified / not verified*).

**Why it matters to us / tasks (when we get to it):**
1. **Load-perturbation evaluation scenario:** replay condition (c)/(d) as an *evaluation* of
   the tuned network — our stance-gated Ib load sharing (LBIN/ib_rge pathway) exists but has
   only ever been exercised by self-generated loads; their loaded-walking kinematics/GRF/EMG
   give the reference for "gait should degrade *this* gracefully". Concretely: a
   `kine_ref`-style compare of our sim vs their loaded references, or an optuna eval variant
   that adds torso mass in the MuJoCo model (patch_xml-level change).
2. **Metabolic-cost references:** if `optuna_walk` ever gains a metabolic term, their
   per-device/per-condition metabolic results are the sanity anchor.
3. Their speed-matched condition pair (a vs d) isolates *load* from *speed* — a clean
   factorial for the DRIVE-vs-load response of our RG.

**Papers:** Dembia CL, Silder A, Uchida TK, Hicks JL, Delp SL (2017) *PLoS ONE* 12(7):e0180320.

## 5. runningsim — Hamner et al. 2010 — **relevance: LATER**

**Page:** simtk.org/projects/runningsim (fetched 2026-09-23). Group 516.

**What it is:** the original single-subject, single-speed (3.96 m/s) muscle-actuated running
simulation — "92 musculotendon actuators representing 76 muscles of the lower extremities and
torso", arms included. Paper: Hamner, Seth, Delp, *J Biomech* (2010), doi:
10.1016/j.jbiomech.2010.06.025. Workflow Scale/IK/RRA/CMC, "after testing in OpenSim 2.0".
Note: our robotbody stack is also a 92-muscle model (AGENTS.md) — the matching count suggests
a shared gait2392-lineage lower-limb muscle set, which would make name-mapping via
`muscle_map` trivial; *that lineage link is my inference* (page does not state gait2392
provenance).

**Contents (from the download listing):** package "01 Simulation of Human Running":
`RunningSimulation_simTK.zip` (34 MB), `Hamner2010_SuppMaterial.zip` (75 MB),
`README_RunningSimulationInfo.pdf`; package "02 Full Body Model":
`FullBodyModel_SimpleArms_Hamner2010_Markers_v2_0.osim` (661 KB) + `hat_ribs_scap.vtp`
(168 KB; torso geometry). Package 01 totals 110 MB per the page's own metadata. **≈ 110 MB.**

**License:** CC BY 3.0 (Stanford 2012, both packages' license popovers). **Anonymous
download:** NO (verified on this project's confirm link). Stats: 29,756 downloads, last
updated Jun 30, 2025 (page stats box).

**Why it matters to us / tasks (when we get to it):**
1. **Cross-model moment-arm check:** the 92-actuator `.osim` gives an independent OpenSim
   moment-arm reference for the same muscle families; compare against our
   `bsolve_ik.py`-style fd_moments along matched poses (a robustness check on the converter's
   pathpoint handling, independent of subject01).
2. **Second reference gait for `fit_synapses`:** CMC activations from a running trial → NMF
   synergy comparison walk-vs-run on real data (they ship the CMC outputs in package 01 —
   verify in the README).
3. Running reference cycle for `kine_ref` if the lab ever targets running (lowest urgency).
4. Format caveat: OpenSim 2.0-era `.osim` — loading under our OpenSim 4.3/4.6 may need the
   GUI's model upgrade path (*inference*, not tested here).

**Papers:** Hamner SR, Seth A, Delp SL (2010) *J Biomech*, doi:10.1016/j.jbiomech.2010.06.025.

## 6. crouchgait — Steele et al. 2010-2017 — **relevance: DEFERRED (Ben's call)**

**Page:** simtk.org/projects/crouchgait (fetched 2026-09-23). Group 509.

**Ben's note (from the goal brief):** cerebral-palsy gait — "I don't think this will be
helpful yet"; catalogued and parked.

**What it is:** OpenSim simulations of children with CP walking in crouch gait: (1) mass-center
acceleration contributions, (2) muscle-weakness impacts, (3) tibiofemoral contact forces,
(4) ankle-foot-orthosis evaluations. Five publications: Steele 2010 *J Biomech* 43(11):2099-2105;
Steele 2012 *Gait Posture* 35(4):556-560; Steele 2012 *J Biomech* 45(15):2564-2569; Steele 2013
*Gait Posture* 38(1):86-91; Rosenberg & Steele 2017 *PLOS ONE*.

**Contents (from the download listing):** package "Crouch Gait Simulations": C1–C10 zips,
14–24 MB each (Jun 2010, 10 subjects, single-limb stance; 167 MB); "Crouch Severity
Simulations": MI01-03 / MO02-04 / SE01/02/05 zips, 23–41 MB + `Subjects.xlsx` (Dec 2012;
254 MB); "Simulated Ankle Foot Orthoses": GIL01/03/04 + MI/MO/SE re-runs, 4–6 MB each
(Jul 2017; 64 MB). **≈ 485 MB total.**

**License:** CC BY 4.0 (Stanford 2016). **Anonymous download:** NO (inferred — same
mechanism; see below). Stats: 17,095 downloads; last updated Dec 3, 2020.

**Deferred — but the two later-use hooks worth remembering (my inferences, not plans):**
1. Steele 2012 "How much muscle strength is required to walk in a crouch gait?" is a template
   for *gait degradation under force-budget shortfall* — directly analogous to our BPA
   torque-budget question (ankle margin 1.4×, AGENTS.md) if we ever characterize how the
   SNS walker degrades when BPA Fmax is cut.
2. The AFO simulations (Rosenberg & Steele 2017) are the published precedent for external
   assistance layered on a muscle-driven walker — relevant if the BPA legs are ever evaluated
   as assistance devices.

---

## Ben's manual download list (SimTK login required — free account)

| Priority | Project page → package | Files | ~Size |
|---|---|---|---|
| 1 | simtk.org/projects/predictivesim → "Supplemental Data" | `ExperimentalData.zip` + `pdsim-source.zip` | 17 MB |
| 2 | simtk.org/projects/muscfib_walkrun → "MuscleFiber Simulation Results" | `README.txt` first, then `Subject 01 Walk.zip`, `Subject 02 Walk.zip` (pilot; all five Walk zips = 340 MB) | 6 KB + 132 MB |
| 3 | simtk.org/projects/nmbl_running → "Raw EMG Data" + "Subject 01" | `RAW_EMG_DATA.zip` + `subject01.zip` | 211 MB |
| 4 | simtk.org/projects/runningsim → both packages | `RunningSimulation_simTK.zip`, `02 Full Body Model` | 110 MB |
| 5 | simtk.org/projects/assistloadwalk | start with `…simulations_of_experiments.zip` (646 MB) when the load-perturbation eval is scheduled | 646 MB |
| — | simtk.org/projects/crouchgait | none (deferred per Ben) | — |

Suggested landing folder (not created — Ben decides): keep zips OUTSIDE the repo (size-reduction
standing objective, AGENTS.md), e.g. `D:\SimTK_GaitLib\<project>\`, and commit only derived
references (speed-binned cycle tables) into the repo.

## Prioritized shortlist (what each buys, ranked)

1. **muscfib_walkrun Walk zips** — speed-binned reference cycles (duty/cadence-vs-speed
   targets for DRIVE calibration: our duty + cadence gaps), F-L-V validation of our converted
   muscles (the AGENTS fidelity open question), human synergy speed-invariance evidence for
   `fit_synapses`. Cheapest big win; note CC BY-NC 3.0.
2. **predictivesim (both files)** — 17 MB for a reflex+contact-state walking-controller
   reference implementation + a second experimental walking dataset; validation-methodology
   framing (mean±SD envelopes) for the dissertation.
3. **nmbl_running (selective)** — raw EMG speed series for afferent-side validation +
   running-regime cadence curve + 10-subject variance bands; 2.06 GB full, 211 MB selective.
4. **assistloadwalk** — load-perturbation evaluation scenario (stance-gated Ib pathway's
   first external test) + metabolic references; ~2 GB.
5. **runningsim** — independent 92-actuator moment-arm cross-check + running CMC reference;
   110 MB.
6. **crouchgait** — deferred per Ben; revisit only for weakness→gait-degradation template or
   AFO-style assistance precedent.

## Verified / not verified

**Verified in this session (2026-09-23):**
- All six project home pages fetched and quoted (webReader tool) — descriptions, citations,
  download counts, last-updated dates.
- All six download listings fetched via `curl https://simtk.org/frs/?group_id={516,509,603,769,1082,728}`
  → HTTP 200; file names/dates/sizes extracted by a parser script from the saved HTML
  (`%TEMP%\parse_simtk.py` run with the myo env python).
- Licenses read from the listing pages' license popovers/JSON-LD: runningsim CC BY 3.0;
  crouchgait CC BY 4.0; nmbl BSD-style (Stanford 2013/2014); muscfib **CC BY-NC 3.0**;
  assistloadwalk BSD (Dembia 2016); predictivesim BSD (Stanford 2015).
- **Login gate:** `curl` on `download_confirm.php` for muscfib README, runningsim README PDF,
  and predictivesim ExperimentalData.zip → each returned a 200 stub whose body is
  `window.top.location='/account/login.php?triggered=1&…'`. Direct `download.php` URL → 404.
- Local pipeline hooks read: `kine_ref.py` (full file), `bsolve_ik.py` (lines 1-120),
  `fit_synapses.py` (full file) — line citations above are from these reads.

**Not verified / stated as inference:**
- **Zip contents beyond the listing descriptions** (exact file formats inside every package):
  login-gated; every "expected contents" statement above is labeled and must be confirmed
  against each project's README after Ben downloads.
- **Login gate on crouchgait, nmbl_running, assistloadwalk specifically:** simtk.org began
  refusing connections from this client (timeouts, HTTP 000) before I could test them —
  almost certainly rate-limiting after my fetch burst. The identical mechanism on their
  listings makes the same gate highly likely, but treat those three as *inferred*, and note
  the project pages themselves fetched fine earlier.
- The runningsim↔gait2392 lineage link, OpenSim-2.0→4.x model-upgrade need, C3D→IK
  preprocessing details, EMG processing steps, and Schutte-vs-Thelen comparison protocol are
  my inferences, each labeled in-line.
- No literature PDFs were fetched — paper claims rest on the abstracts quoted on the SimTK
  pages. Per the supervisor standard, any wiring claim derived from `pdsim-source` must wait
  until the source is actually read (task 2.2).
- Nothing was downloaded and no code was modified (`mkdir reports_20260923` + this file are
  the only filesystem writes in the repo).

## Addendum 2026-09-23 (post-run integration)

- No open mirror exists for the Arnold muscfib data (SimTK-only) - Ben's
  login remains required for items 1-2 of the manual download list.
- OPEN ALTERNATIVE STAGED: Falisse/Afschrift/De Groote 2022 predictsim_mtp
  (github.com/antoinefalisse/predictsim_mtp, 46 MB) downloaded anonymously
  to D:\temp\gait_lib_staging\falisse\ (outside the repo). Contains full
  predicted walking trajectories WITH 92 per-muscle activations
  (Results/Case_40/motion.mot + GRF.mot, one periodic cycle), scaled OpenSim
  models, and an averaged experimental walking IK template. Not one of the
  six SimTK-listed libraries - Ben's call whether to adopt it as the
  predictive-sim reference.
- NEW loader spinal/gait_lib_loader.py: any (IK .mot + vertical-GRF .mot)
  pair -> kine_ref reference schema; GRF column auto-detect (subject01 vs
  r_/l_ prefixes) + periodic-wrap mode for single-cycle predictive data.
  VERIFIED: subject01 through the loader reproduces
  kine_ref.load_reference() to < 1e-9 (REGRESSION PASS). Falisse Case_40
  reference: T 1.113 s, duty 0.57/0.57, lag_rl 0.50, ds 0.15, knee_min
  -59.2 deg, hip range 51.6 deg (vs subject01 1.233 s / 0.61 / -69.7 deg).
- Synergy transfer check (goal5 x goal7): our six subject01-derived
  synergies (synergy_basis.npz W) explain 0.596/0.595 R2-vs-zero of the
  Falisse per-leg activations but VAF (mean-removed) is about -0.02: the
  bases do NOT transfer across datasets. Implication: train against
  multiple reference sources; do not fix one synergy basis. Smoke log:
  goal5_loader_smoke.txt.
