# Balance / gait data catalog — external candidates vs our local library (2026-09-24)

Session: EB475WS4, 2026-09-24. Source of external candidates:
https://github.com/modenaxe/awesome-biomechanics (raw README at
`master/README.md`, 189 KB, downloaded and read directly — Balance
section line 452, Walking line 479, Running line 530, EMG line 655).
Licenses were verified LIVE this session against Zenodo/figshare APIs
(commands in Provenance section) — not taken from the README, which
states no licenses.

Deliverable question: which open datasets could serve
**(a)** balance-perturbation training references, or
**(b)** frontal-plane (medial/lateral GRF) validation for the s3k walker.

---

## 1. What we already have locally (verified this session)

`Code\MuJoCo_SNS\spinal\gait_refs\` — **43 npz**, verified by `dir` and by
an inspection script run under `C:\Users\Ben Bolen\.conda\envs\myo\python.exe`:

| family | count | source (npz `grf_style` tag) | gait type |
|---|---|---|---|
| `Case_40_motion.npz` | 1 | `falisse` (Falisse 2022 predicted, predictsim_mtp) | WALKING |
| `ong_speed_050..200` + `ong_selfsel_Init050/125/200` | 10 | `ong_states` (Ong predictive walking states) | WALKING (speed sweep) |
| `subject01/02/04/08/10/11/17/20_Run_*` | 32 | `arnold` (Arnold 2010-style ExperimentalData, 8 subj × 4 speeds) | RUNNING (all duty 0.31–0.46) |

Every npz (verified structure, one per family) holds only: 6 cycle-normalized
100-point joint-angle traces (`{r,l}_{hip,knee,ankle}`) + per-leg scalars
(`duty_*`, `T_*`, `hip_range_*`, `mean_*`, `knee_range_*`, `ankle_range_*`,
`knee_min_*`, `lag_rl`, `ds`) + `grf_style`/`periodic`. **No GRF waveforms and
no EMG are stored in any ref.**

The subject01 *walking* "of record" reference is NOT an npz — it is
`kine_ref`'s built-in cycle (subject01_walk1_ik.mot), scored in-memory in
`gait_lib_score_all.py:48-54`.

**Two capability gaps for this ask, both verified:**

1. **No perturbation references at all.** Nothing local contains a
   perturbation response.
2. **No medial/lateral GRF anywhere in the pipeline.**
   `gait_lib_loader.py:27-30` — `GRF_STYLES` uses **vertical** columns only
   (`ground_force_vy`, `r_/l_ground_force_vy`); `_detect_grf`
   (`gait_lib_loader.py:55-63`) additionally accepts the Arnold duplicated-`vy`
   convention. M/L columns are never read; and refs store no GRF anyway.

Best current s3k walking scores (from `reports_20260923\goal5_allrefs_scoring.md`,
this library): subject01 of-record −189.5, Falisse −200.5, Ong −203…−224;
running refs −188…−332 (monotone in speed). The scoring is joint-angle kine
only — GRF-based scoring does not exist yet.

---

## 2. Balance section of awesome-biomechanics (all 3 entries)

| dataset | contents | n | license (verified) | (a) perturbation | (b) M/L GRF |
|---|---|---|---|---|---|
| **Wang & van den Bogert 2020** — Standing Balance with Long Duration Random Pulses (Zenodo 3819630) | 32 markers (27 body + 5 treadmill frame), **6-DOF GRF** (3 F + 3 M), 9 EMG right leg; **plus processed joint angles & torques via inverse dynamics**; 5-min trials | 8 subj × 4 trials (2 quiet + 2 perturbed) = 80 min quiet + 80 min perturbed | **CC BY 4.0** (Zenodo API: `cc-by-4.0`, open; Raw 2.9 GB + Processed 264 MB + MATLAB processing code + experiment PDF) | **YES — flagship.** Random-pulse standing perturbation reactions, explicitly recorded to identify postural feedback controllers | Yes (6-DOF), standing |
| **BDS** — dos Santos & Duarte 2016 (figshare 3394432) | force-platform F/M/COP posturography + Mini-BEST scores, 4 conditions (eyes × surface) × 3 reps | 163 subj | **CC BY 4.0** (figshare API) | No (no perturbations, no kinematics) | No waveforms needed — COP-based norms only |
| **PDS** — dos Santos et al. 2017 (figshare 4525082) | 42-marker whole-body 3D kinematics + **73 joint angles** + dual force plates (F, M, COP, resultant), 60-s standing under manipulated vision & surface; 12 trials/subj | 49 subj (27 young + 22 old) | **CC0** (figshare API; 7.4 GB) | No pulses — static/sensory-condition balance only | Yes, standing (dual plates) |

Wang's Zenodo description (fetched this session): *"fundamental for
identifying postural feedback controllers … long duration balance data (under
random external perturbations) … joint angles and torques were calculated
using a human body model and inverse dynamics."* — i.e. **v1 use needs no IK
pipeline**: processed angles/torques ship with the data.

Perturbation direction for Wang: the perturbing device is the treadmill the
subject stands on (5 frame markers); the record description does not state the
direction explicitly — treat as ground-level translation at the feet, direction
to be confirmed from the raw files at download (not verified this session).

## 3. Walking / gait section (promising entries)

| dataset | contents | n | license (verified) | (a) | (b) |
|---|---|---|---|---|---|
| **Moore et al. 2015** — elaborate gait + mechanical perturbations (Zenodo 13030) | instrumented treadmill with movable base, full marker set (human + treadmill), belt speeds/accelerations, **forces & moments from dual plates**; 120 s normal + 480 s perturbed per trial at 3 speeds; longitudinal (belt) perturbations | 15 subj × 3 speeds | **CC0** (Zenodo API; 2 tar.gz ≈ 2.4 GB) | **YES — walking perturbations**, sagittal/longitudinal only; no EMG in the record's file description | Yes (3D F&M), but perturbations are AP |
| **Schreiber & Moissenet 2019** (figshare 7734767) | 52-marker trajectories + **3D GRF & moments + EMG recorded simultaneously**, 5 speeds | 50 healthy adults | **CC BY 4.0** (figshare API; 383 MB) | No | **YES — best experimental M/L GRF walking set** + EMG + moments; natural companion to our Ong *predicted* speed sweep (adds experimental + neuromuscular dimensions) |
| **Fukuchi et al. 2018 WBDS** (figshare 5722711) | overground + treadmill walking kinematics + kinetics across speeds | 42 (24 young, 18 older) | **CC BY 4.0** (figshare API; 2.2 GB) | No | Yes (3D GRF incl. M/L); adds AGE spectrum |
| **GaitRec** — Horsak et al. 2020 (figshare collection 4788012) | **GRF only**, bilateral walking trials, shipped as separate RAW + PRO (processed) items per component: `GRF_F_V_*` (vertical), `GRF_F_AP_*` (anterior-posterior), **`GRF_F_ML_*` (medio-lateral)** | 2,084 patients + 211 controls = 75,732 trials | **CC0 on all 10 collection items** (verified individually via figshare API; data items ≈ 137–261 MB each) | No | **YES — largest M/L GRF corpus in existence**, healthy controls for norms + patient spectrum (impaired M/L patterns) for robustness targets |
| **Gutenberg Gait Database** — Horst et al. 2021 (figshare collection 5311538) | GRF + COP of two consecutive overground steps, self-selected speed, two embedded plates; components shipped as V / ML items + metadata | 350 healthy | **CC0 on all 10 collection items** (verified; data items ≈ 16–25 MB) | No | **YES — cheapest normative M/L GRF** (~100 MB total) |
| Lencioni 2019 (level/toe/heel/stairs, kin+kin+EMG) | multi-task walking | 50 subj, ages 6–72 | not checked this session | No | unverified |
| Luo 2020 uneven-surface IMU; Miraldo/GEDS 2020 IMU+contact+EMG (9,661 strides, CC0, 4.3 GB); Winter 2-D; Kirtley norms; Liu 2008 SimTK | wearables / legacy | — | GEDS CC0 (API); others not checked / SimTK registration | No | No (no lab-grade M/L GRF) |

Caveat on both figshare collections: the API currently lists **10 items
each**, and the returned listing is asymmetric (e.g. GaitRec PRO appears only
for `_left`; Gutenberg shows right-side V RAW+PRO + ML RAW only). Verify
per-side/per-component availability at download time.

## 4. EMG / neuromuscular section (Santuz family)

| dataset | contents | n | license (verified) | use |
|---|---|---|---|---|
| Challenging settings 2020 (Zenodo 3785065) | raw + filtered EMG, timings, **NMF synergies**, sMLE/HFD code | 476 trials / 86 participants | **CC BY 4.0** | synergy TARGETS under challenging locomotion — complements `synergy_basis.npz` (6-synergy fsa_backsolve model); note our six-synergy basis did NOT transfer to Falisse (VAF≈0) — train multi-reference, don't fix one basis |
| Sex-specific tuning 2022 (5171823); high-speed 2020 (3785077); treadmill-vs-overground 2020 (3932768) | EMG + synergies | 215/30/30 participants | **CC BY 4.0** each | EMG-level validation only |
| Running 2018 (3785076) | EMG + synergies | 135 adults | **CC BY-SA 4.0** | running only |

No kinematics/GRF in these — EMG/synergy comparison dimension only.

---

## 5. TOP-5 shortlist and why

1. **Wang & van den Bogert 2020 standing balance random pulses** (CC BY 4.0,
   Zenodo 3819630) — the only open set with *simultaneous* markers + 6-DOF
   GRF + EMG during random-pulse perturbations, recorded precisely to identify
   postural feedback controllers, and it ships **processed joint angles and
   torques** (no IK needed for v1). Direct feed for the goal-2 balance stage
   (`--stand-eval` / VEST work, `reports_20260923\goal2_balance_stage.md`):
   pulse-response references for CoP/CoM recovery, ankle/hip-strategy EMG
   targets. Fills our #1 gap: zero perturbation references locally.
2. **GaitRec** (CC0, figshare c.4788012) — the only dataset shipping
   explicitly named **medio-lateral GRF components** (RAW + PRO), at a scale
   nothing else approaches (75,732 trials; 211 healthy + 2,084 patients).
   GRF-only format sidesteps IK entirely; PRO items are processed bilateral
   traces ready for a normative M/L scoring dimension. Fills our #2 gap.
3. **Schreiber & Moissenet 2019** (CC BY 4.0, 383 MB) — 50 adults × 5 speeds
   with 52-marker kinematics, 3D GRF&M **and EMG simultaneously**: the best
   single experimental upgrade to our speed-sweep references (Ong refs are
   predicted; subject01 is one experimental cycle), and it adds both the M/L
   GRF and the EMG dimensions we lack.
4. **Moore et al. 2015 perturbed walking** (CC0, Zenodo 13030) — 480 s of
   perturbed walking per trial gives the statistics for recovery-step training
   references during gait (Wang is standing). Dual force plates give full
   F&M. Caveats: perturbations are longitudinal only (no M/L pushes), no EMG,
   and joint-level refs need marker IK (the `bsolve_ik.py` chain exists).
5. **Gutenberg Gait Database** (CC0, ~100 MB of data items) — cheapest, fastest
   full-library win: 350 healthy subjects' V + M/L GRF + COP of two-step
   trials. Enough for a normative M/L symmetry/frontal-plane metric for the
   walker without any kinematics processing.

**Honorable mentions:** PDS (CC0 — standing balance across vision/surface
conditions with 73 joint angles + dual plates: natural reference set for a
standing-balance evaluation stage, not perturbation training); Santuz
"challenging settings" (CC BY 4.0 — synergy targets under demanding
conditions); WBDS (age spectrum); BDS (CoP-only clinical norms).

**Explicit gap:** nothing in these sections provides **medio-lateral push /
slip perturbations during walking** — Moore is belt (AP), Wang is standing.
If frontal-plane *perturbation training* (not just M/L validation) becomes the
goal, it needs a source outside this README.

---

## 6. Integration notes for our pipeline

- **Refs are angle-only today.** M/L-GRF scoring requires (i) loader changes —
  extend `GRF_STYLES` (`gait_lib_loader.py:27`) with per-dataset column
  conventions and stop discarding non-vy columns; (ii) npz schema extension
  (add normalized Fx/Fz-Fy-side traces alongside the 6 angle arrays);
  (iii) sim-side extraction of per-foot M/L contact forces from MuJoCo contact
  wrenches (mind the predefined-contact-pair surgery noted in AGENTS.md).
- **Format fit.** Wang (txt mocap/analog + processed tables, MATLAB code
  included) and Moore (tar.gz, text-based) fit the `gait_lib_intake.py`
  staging pattern (`D:\temp\gait_lib_staging\downloads` → loader). Schreiber /
  WBDS / GaitRec / Gutenberg ship their own formats (not verified per-file
  this session) — each needs a small reader before the `.mot`-pair pattern
  applies.
- **Curriculum hooks already present:** `--stand-eval`, `--vest` family,
  `contact_onset` knob (2026-09-23 build), and the multi-reference training
  pattern (`gait_lib_pilot2.py` REF_CACHE swap; s3k vs library scoring in
  `gait_lib_score_all.py`).

## 7. Provenance / what was and wasn't verified

Verified LIVE this session (queries run via `urllib` scripts under the myo
env, `C:\Users\Ben Bolen\.conda\envs\myo\python.exe`; raw README downloaded
with `curl` and read directly):

- Local library: `dir` of `gait_refs` (43 npz) + npz structure dump (one per
  family) + `grf_style` tags.
- Loader GRF handling: read `gait_lib_loader.py:20-69`.
- Licenses: Wang `cc-by-4.0` + file list (Zenodo API 3819630); Moore
  `cc-zero` + description (13030); PDS `CC0` 7.4 GB (figshare API 4525082);
  BDS `CC BY 4.0` 376 MB (3394432); WBDS `CC BY 4.0` 2.2 GB (5722711);
  Schreiber `CC BY 4.0` 383 MB (7734767); GaitRec all 10 items `CC0`
  (collection API 4788012); Gutenberg all 10 items `CC0` (5311538); GEDS
  `CC0` 4.3 GB (7778255); Santuz 5171823/3785077/3785065/3932768 `cc-by-4.0`,
  3785076 `cc-by-sa-4.0` (Zenodo API).
- PDS contents quoted from its abstract (Europe PMC core search).

Not verified this session (stated as such above): file formats inside
Schreiber/WBDS/GaitRec/Gutenberg; Wang perturbation direction; Lencioni /
Luo / Winter / Kirtley / Liu / Fukuchi-running licenses and contents
(README descriptions only); SimTK items (registration required, not
checked).
