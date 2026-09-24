# Thumb-drive biomechanics campaign (2026-09-24, overnight)

Ben supplied `F:\Biomechanics data` (19 items, 3.03 GB) — the SimTK
login-gated haul from goal 5. The drive was READ ONLY throughout.

## What arrived

| item | what it is |
|---|---|
| SubjectNN-latest.zip x8 (01,02,04,08,10,11,17,20) | Arnold-cohort RUNNING trials: full toolchain per trial (IK in `ik/Results_191/`, GRF exports in `ExportedData/` with duplicate plate column names, per-cycle RRA/CMC states) — 4 speeds each |
| RawEMGData-latest.zip | RAW EMG .sto, 10 subjects x 4 running trials (matched cohort + subjects 03/19) |
| results-speeds.zip + results-deficits-{contracture,weakness}.zip | Ong et al. SCONE+CMA-ES predictive WALKING: 0.50–2.00 m/s x7 + 3 self-selected + gas/pf/sol x mild/moderate/severe deficits; gait9dof18musc model; 1 ms states incl. fiber lengths (states.sto in meters) |
| scone-setup-files-v2.zip | the SCONE setups (.scone/.par) that generated the above — also the reference for the goal-2 standing-balance stage |
| assistloadwalk_simulations_of_experiments.zip | 647 MB, 5270 .sto — assisted-loading walking experiments (inventory; separate states route) |
| RunningSimulation_simTK.zip + Hamner2010 supp + FullBodyModel osim/vtp | Hamner 2010 full-body running model + subject02 running IK/GRF + mass properties |
| ModelWithSampleSimulations-4.0.zip | Rajagopal 2015 model + CMC sample sims |

Staged: `D:\temp\gait_lib_staging\downloads\` (zips + extracted trees,
17.6 GB total, D: has 505 GB free). Repo untouched by the data.

## Intake (gait_lib_intake.py v2 — upgraded tonight)

- recursive zip-in-zip extraction (SubjectNN-latest -> subjectNN),
  `__MACOSX` skipped
- NEW GRF style "arnold" in gait_lib_loader: OpenSim plate-2 columns
  share names positionally; sides assigned by earlier first loading
  onset (fix for "no known GRF vy column pair")
- cross-folder trial pairing (ik/Results_191/Run_XXXXX <-> ExportedData/
  "Run_X XX_newCOP3") via trial tokens; RRA/CMC derivative states
  de-duplicated out (they self-pair angles+GRF)

**33 references integrated** into `spinal\gait_refs\*.npz`:
Falisse Case_40 predicted WALKING + 8 subjects x 4 RUNNING speeds.
All 32 Arnold trials are RUNNING (duty 0.31–0.46, T 0.56–0.80 s,
knee_min −80…−140 deg — a clean speed sweep, not walking).

## s3k production walker vs the whole library

`gait_lib_score_all.py` -> goal5_allrefs_results.csv / _scoring.md.

- of-record subject01 WALKING: **−189.5** (best)
- Falisse predicted WALKING: −200.5
- slowest running (Run_2xx): −188…−208 (subject08 best at −188.2)
- fastest running (Run_5xx): −277…−332 — monotone degradation with speed

Reading: the walker (duty 0.25, T 0.90 s, knee −49 deg) matches walking
references best and degrades smoothly with running speed — expected for
a walker-shaped gait, and now QUANTIFIED across 33 human references.
(cycles_l=0 is the known eval-npz quirk, identical to the 09-23 run.)

## Integration routes still open (not blocking)

- Ong/SCONE predictive walking states (0.5–2.0 m/s) = single-cycle,
  GRF-less -> needs the states-based loader extension (phase by pelvis
  or hip sign + periodic-wrap), then they become WALKING references
  spanning the speed range between subject01 and running.
- RAW EMG (10 subj x 4 trials) -> activation-level comparison with our
  synergies / back-solve pipeline (pair with the kinematic refs by
  trial token).
- RunningSimulation subject02 pair (IK/subject02_running_arms_ik.mot +
  RRA grf) didn't auto-pair (no Run_ token, not under a subjectNN
  folder) — one manual ref if wanted.

## Verdict

The multi-reference training corpus goal 5 needed is NOW LOCAL:
2 walking + 32 running refs live, scorer verified, REF_CACHE swap
(gait_lib_pilot2.py pattern) turns any of them into a training
objective. Next study (awaiting Ben): multi-reference objective that
samples refs per trial — walking refs for the current walker, running
refs only if/when he wants a running gait.


## Addendum (same night): Ong/SCONE predictive WALKING refs integrated too

The states.sto files turned out to carry per-leg vertical GRF IN the
states (Leg1_r.grf_y / Leg0_l.grf_y, BW-normalized) plus gait2392-named
angles in radians - no states-route workaround needed. New
`gait_lib_ong_refs.py` (own radians-tolerant .sto reader; GRF scaled
x750 N so the fixed 50 N onset detector applies; case-insensitive
column match for the SelfSelected vintage) built **10 more walking
references**: ong_speed_050..200 (T 2.10->1.02 s, duty 0.62->0.53,
knee ~-70 deg) + 3 self-selected. Caveat: ong_selfsel_Init200 has a
degenerate left cycle (duty_l 0.00) - drop it if it ever wins a
sampling draw.

**Library total: 43 references** (1 subject01 of record + Falisse
predicted walking + 8 subjects x 4 running speeds + 10 Ong predicted
walking speeds). s3k walker vs the Ong family: -203 (0.5 m/s) to -224
(2.0 m/s), monotone - the walker sits in the slow-gait cluster of the
manifold (subject01 -189.5 best, running sprint -332 worst).
Deficit/contracture runs intentionally NOT integrated (Ben: not
helpful yet); they use the same script pattern when wanted.
