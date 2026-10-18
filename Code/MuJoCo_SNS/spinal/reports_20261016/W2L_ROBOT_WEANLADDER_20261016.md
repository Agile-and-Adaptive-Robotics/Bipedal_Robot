# W2L ROBOT FEEDBACK OPTIMIZER + DRIVE-WEAN LADDER — stage report (2026-10-16, egg)

Overnight-campaign stage: **"feedback optimizer (robot)"** (phase 8 of the
11-phase plan in `CAMPAIGN_HANDOFF_20261016.md`). All work ran ON
easteregg2 (repo `D:\GitHub\Bipedal_Robot`, branch KneeTestSetup_BenBo_stw)
via ssh; python `D:\Anaconda\envs\myo\python.exe`. NO commits, NO pushes.
All outputs: `reports_20261016\w2l_wean\` (studies in `optuna_w2l_robot.db`,
per-trial diagnostics `wean_trials.csv` — 402 rows).

The ask: *same wean-ladder recipe as the human side — optuna study (new
name `curr_s6_drive_robot*`), tune feedback gains (heel reset, toe→TOEDF,
Ib autogenic, commissurals) at high L/R drive, then force the drive scale
down across rounds and report how far it weaned while alternating
progression survived (>= 6 cycles); kine = the progression/antiphase
metric this stack already reports; basin-gate any finalist; anti-lottery
per-trial diagnostics + finalist gain-sensitivity probe + wean ladder as
the structural test.*

## 1. VERDICT

**The W2L robot weans from 100% to 60% of its bootstrap drive
(te=4→2.4, tf=5→3.0 nA) with alternating supported progression INTACT,
robustly (fresh-process reproduction + 12/12 basin), and the weaned
walker is BETTER than the drive-only bootstrap on every reported axis
except lean.** Below 60% the gait dies: 0.50 breaks antiphase
(inter-leg phase 0.322), <= 0.40 produces zero heel strikes.

| metric (30 s walk window, CLI-canonical) | bootstrap (dscale 1.0, Ben's gains) | **WINNER (dscale 0.60, tuned)** |
|---|---|---|
| primary stride cycles L/R | 22 / 21 | **25 / 23** |
| inter-leg phase median | 0.553 | **0.487** (near-perfect antiphase) |
| RG-E L/R correlation | −0.732 | −0.766 |
| sustained speed (fit, 8–30 s) | 0.0609 m/s | **0.0662 m/s** (+9%) |
| forward displacement | +1.644 m | **+1.718 m** |
| combined foot duty L/R | 0.061 / 0.158 | **0.166 / 0.200** |
| tilt max | 12.0 deg | 16.3 deg |
| harness weight carried | 28% | 16% |
| objective (this study) | 91.4 | **94.9** |

Winner knobs (full file `w2l_wean\w2l_wean_best_params.json`):
te=2.4, tf=3.0 (dscale 0.6), g_heel=0.588, g_ib=0.112, ibnA=1.321,
ctn_amp=0.717, cref=46.4, comm=0.413, toe_m=2.220 (g_toe_in/df/df_mn =
2.220/11.099/6.102), legcoll=0, sink=0.06, rig=0.25, walker_yield=1.5.

Reproduce (canonical, 5/5 identical runs including the supervisor check):

```
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_wean && D:\Anaconda\envs\myo\python.exe _w2l_wean_eval.py --id=reproduce_finalist --dur=30 --dscale=0.6 --g-heel=0.58765 --g-ib=0.11164 --ibnA=1.32101 --ctn-amp=0.71708 --cref=46.37426 --comm=0.41304 --toe-m=2.21976
```
prints `RESULT {json}`: cycL 25, cycR 23, phase_med 0.487, slope 0.06617,
fwd +1.718, survive true, obj 94.86191 (deterministic in this context).

The raw test-CLI command (same knobs, full precision — in the param file)
lands a SIBLING healthy walk in its own process context (heel strikes
30/30, duty 0.25/0.27, fwd +2.196 m, tilt 13.4, deterministic 2/2): both
contexts pass the criterion, and the config difference was proven ZERO
(`_kv_diff_probe.py`: the only differing KV key is `phase`, which
`run_phase` never reads) — the split is the stack's documented
BLAS-summation-order attractor selection, not a config error.

## 2. THE LADDER (optuna curr_s6_drive_robot_r100…r010, 26–36 TPE trials
per rung at dur 20 s, seeded with the previous rung's winner; every rung
winner then re-run at 30 s; **CLI-reproduction audit** below)

Search space = feedback only (protocol/body pinned at the bootstrap
winner): g_heel [0.05,2.5], g_ib [0,2.5], ibnA [0,3], ctn_amp [0.5,6],
cref [20,150], comm [0.2,2], toe_m [0,3] (toe chain scaled together,
Ben's 1:5:2.749 ratios preserved). Drive scale s: te=4s, tf=5s (ratio
kept; the 10 nA/10 ms antiphase kickoff is part of the model, unscaled).

| dscale | in-process 30 s verify | **CLI fresh-process audit** | verdict |
|---|---|---|---|
| 1.00 | survive (20/20, 0.075 m/s) | survive (20/20, 0.075) | robust |
| 0.85 | survive (31/27) | survive (31/31, 0.067) | robust |
| 0.70 | survive (21/22) | survive (18/11, 0.036) | robust (thinner R) |
| **0.60** | survive (20/20, 0.085*) | **survive (25/23, 0.066)** | **WINNER — deepest robust rung** |
| 0.50 | survive (31/18) | **FAIL** (phase 0.322, cycR 9) | antiphase breaks |
| 0.40 | survive (11/10, 0.009) | FAIL (0/0 strikes) | dead |
| 0.30 | survive (28/29, tilt 43) | FAIL (0/0) | dead |
| 0.20 | survive (34/33, 0.065) | **FAIL (0/0)** — see §4 | context artifact |
| 0.175 | fail | FAIL | dead |
| 0.15 | fail | FAIL | dead |
| 0.10 | fail | FAIL | dead |

\* in-process numbers carry +Δobj context offset (§4).

**Death scale (structural test the ask requested): alternating progression
dies between dscale 0.60 and 0.50 — at 0.50 by antiphase loss, at <= 0.40
by total gait loss (R foot planted, 0 heel strikes).**

## 3. BASIN GATE + GAIN-SENSITIVITY (the anti-lottery requirements)

- **r060 winner: BASIN PASS 12/12** — all seven tuned scalars jointly
  perturbed ±1% (seeded 12 draws), full 30 s re-runs; min cycles 20,
  min fwd +2.145 m (`basin_results_r060.json`).
- **Sensitivity, one pathway removed at a time (r060, 30 s)**: heel reset
  →0: obj 93; Ib→0: 101; Ib current→0: 101; contact encoders→0: 86;
  comm→0.001: 90; toe chain→0: 94 (anchor 101). No single pathway is
  solely load-bearing at 0.6 drive — **but Ben's untuned default gains at
  the same drive are DEAD (obj 7, cyc 10/0)**: the tuned SET is what
  sustains the weaned gait (distributed, real tuning).
- r100 (full drive) for comparison: basin PASS (min cyc 19); comm
  silencing hurts most there (obj 40); defaults survive (94 = the
  bootstrap itself).
- r020 (the context artifact): every single pathway IS critical there
  (each removal → obj 5–7) AND the basin FAILS (5/12 die) — the classic
  marginal-winner signature; kept as the cautionary datum.
- Per-trial diagnostics: every one of the 402 runs logged antiphase /
  progression / periodicity (cycles, periods, phase median, rg_r, slope,
  duties, tilt, harness) to `wean_trials.csv` AND as optuna user_attrs.

## 4. THE CONTEXT AUDIT (why the CLI column exists — a caught trap)

The campaign driver runs in one python process; its 30 s verifications
"survived" down to dscale 0.20. Reproducing the r020 winner in a FRESH
CLI process gave cyc 0/0 (R foot planted, rg_r +0.30 IN-PHASE) — and its
basin had already failed 5/12. Diagnosis trail: (1) CLI repeats are
deterministic among themselves (4/4 identical, obj 94.86191; the test-CLI
context likewise 2/2 identical); (2) `_kv_diff_probe.py` proves the two
contexts' KV are IDENTICAL for every sim-relevant key; (3) therefore the
in-process 0.50–0.20 "survivors" rode the host process's BLAS
summation-order state — the documented stack bifurcation property
(AGENTS.md optimizer-insights; basin-gate caveat). **All conclusions in
this report use fresh-process CLI reproduction.** Lesson recorded for the
human track: verify rung winners in a clean process before believing a
wean depth.

## 5. WHAT THE TUNING DID (gain patterns, `round_best_r*.json`)

Weaning strengthens the sensory loops: heel-reset gain 0.50 (Ben) →
0.59–0.83 across rungs (1.15 at the marginal 0.20), contact-encoder
amplitude 2.0 → 0.7–3.0 (cref ~46–60 N), Ib current 1.0 → 1.3–1.8 nA at
low drive, commissural scale 1.0 → ~0.41 (the two RGs need LESS mutual
inhibition once heel feedback paces them), toe chain ×1 → ×1.1–2.2. The
kickoff + NaP half-centers + heel/toe/Ib feedback carry the rhythm at
60% drive; the tuned set is necessary (defaults die) and redundant
(no single pathway critical).

## 6. HONEST CAVEATS

- "kine" for this plant = the stack's own progression/antiphase metrics
  (sustained forward speed with >= 6-cycle + phase gates); NO kinematic
  reference exists for the W2L robot (kine_ref is human-stack).
- Duty 0.17/0.20 (combined heel|toe > 5 N) still sits below the campaign
  0.4–0.75 target; the test's own strict WALKS gate still reads
  SHUFFLES-IN-PLACE (needs heel duty >= 0.3).
- The harness is NOT weaned: wsup 0.7 vertical support + rig 0.25
  leash/tilt + yielding anchor (walker_yield 1.5 s) all still on — this
  ladder weaned the NEURAL drive only. (The 0.60 winner does carry less
  harness weight: 16% vs 28%.)
- Heels still double-tap (primary-cycle counting merges at 0.45 s);
  tilt 16.3 deg at the winner (13.4 in the test-CLI sibling context).
- r070's CLI audit is thin on the R side (18/11 cycles) — it survives,
  but 0.60 is the rung I would hand anyone.

## 7. FILES / CODE CHANGES (working tree only; see CHANGES_MANIFEST_20261016_w2lwean.md)

- `w2l_mujoco\test_w2l_ground.py` — g_heel/g_toe_in/g_toe_df/g_toe_df_mn/
  g_ib/comm KV knobs passed into `build_w2l_aff_net.build(gains=...)`.
  Defaults == Ben's rules values → build unchanged by construction;
  re-verified bit-exact post-edit: default walk gate (bursts 11/10,
  period 1.539, r −0.776, tilt 8.2, harness 49%, fwd +0.001) AND the
  bootstrap winner command (L33/R41, +1.676 m, 22/21, r −0.743) — both
  identical to the pre-edit runs.
- `reports_20261016\w2l_wean\` — everything else (driver scripts
  `_w2l_wean_eval.py` / `_w2l_wean_campaign.py` / `_w2l_wean_deep.py` /
  `_w2l_wean_mid.py` / `_w2l_wean_basin.py` / `_w2l_cli_verify_rungs.py` /
  `_w2l_finalize.py` / `_w2l_write_wean_params.py` / `_kv_diff_probe.py`,
  db, 402-row CSV, round jsons, basin jsons, `w2l_wean_best_params.json`,
  logs).
- Detach recipe for egg (measured tonight): git-bash nohup, cmd `start`,
  and PowerShell Start-Process ALL die with the ssh session (sshd job
  object); the working route is WMI `Win32_Process.Create` on a launcher
  .cmd (`w2l_wean_run.cmd`).

## 8. NOT RUN (honesty list)

- No AnimatLab realization (phase 9's job).
- No wean of the harness (wsup/rig ladder) — separate experiment, not in
  this ask.
- No optuna of protocol knobs (rig/sink/walker_yield/wsup pinned at the
  bootstrap winner's values).
- The toe→TOEDF path is the M5 DOCUMENTED REDUCTION (terminal inhibition
  lands on the DF MN, no dedicated DF half-center in this net) — tuned as
  a chain multiplier, not rebuilt.

— feedback optimizer (robot), 2026-10-16
