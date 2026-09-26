# figs/ — SCONE figures for the 2026-09-25 walker campaign

Regenerate any figure with: `"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" tools\make_scone_figs.py` (parser: `tools\scone_sto.py`; reads only `logs\motion_reeval_*.sto`, no sims).

- `scone_gait_overfit.png` — THE overfit figure: Tutorial 4a gait pelvis forward distance + height vs time, pretrained vs 4-s-window CMA-ES best vs 20-s full-horizon best; the D4 par falls at 4.68 s / 5.05 m (20-s fitness 83.3 = fell) while pretrained walks 21.48 m and D20 best 21.08 m, both the full 20 s / 36 steps.
- `scone_balance_gains.png` — Tutorial 3a balance pelvis height vs time: default falls at 1.04 s, 40-gen CMA-ES best at 2.81 s, 300-gen best stands the full 30 s at 0.93 m (fitness 101.8 → 94.5 → 2.19).
- `scone_tutorial_scores.png` — horizontal bar chart of the 15 first-ever tutorial evaluations at default params (Hyfydy, sconecmd -l 2): orange = maximize-style (jump 1/2a/2b/2c + ScriptMeasure 6a/6d, higher = better), blue = minimize-style (gait/balance, lower = better); scores not comparable across styles; values transcribed from goal1_scone_tutorials_and_optimization.md §1.

---

# figs/ — goal-4 VARIANT stage-5 winner figures (w2lvar + syn6), 2026-09-25 evening

Both stage-5 winners were REPLAYED through the exact curriculum objective code path
(env pin `AARL_NET` + figs-private `AARL_NPZ=tmp\figs_npz_{variant}.npz`, BASE_MUL merge from
`best_walk_params_v10.json`, `set_stage(5, p)`, `runner --eval --drive <winner>`; NO optuna db
opened, protected `spinal_run.npz`/`optuna_walk.db` untouched). BOTH REPLAYS PASS BIT-EXACT:
w2lvar RECOMPUTED −165.0405055213092 vs json −165.0405055213092 (delta 0.0000 %, kine −165.04,
kz 0.8297, tilt_max 27.72, duty 0.4597, n_r=3/n_l=3, T_r 1.173 s, lag_rl 0.545, bilateral=True);
syn6 RECOMPUTED −197.21456977796157 vs json −197.21456977796157 (delta 0.0000 %, kine −197.21,
kz 0.8396, tilt_max 28.07, n_r=4/n_l=0, T_r 1.877 s, lag_rl 0.395, bilateral=False — left leg
planted, the documented syn6 stage-5 outcome). NOTE: `--eval` rewrites the schedule to 16 s
(1 stand / 1 ramp / **walk 2–13 s** / 1.5 / 1.5, runner.py:901) — the figures use THAT window,
not the stock 5–15 s.

- `w2lvar_s5_walk_overlay.png` — w2lvar stage-5 winner ground-eval mean joint cycles (r+l, from
  each foot's own contact-loading onsets, mean ± 1 sd) vs the OpenSim subject01 reference cycle;
  score −165.0 annotated. Honest content: bilateral but shallow knee (sim ≈ −15° vs ref −70°),
  weak hip modulation, ankle +15° offset — the −165 is a pattern score, not a match.
- `syn6_s5_walk_overlay.png` — same for syn6 (−197.2): right leg cycles (n=4, T 1.88 s, duty 0.74),
  left leg planted (duty 1.00, 0 cycles), ankle PF-biased ≈ −45°.
- `w2lvar_s5_traces.png` — 16 s eval: hip/knee/ankle r+l, pelvis tilt + COM z, per-foot contact
  force, neural raster (DRIVE/POSTURE, RG_E/F r+l, PF_HIP-E/KNEE-E/KNEE-F/ANK-F_r merged
  joint-layer channels, BAL_*), walk window shaded.
- `syn6_s5_traces.png` — same; PF lanes are the **S1..S4 synergy channels** (runner.py:1148),
  burst rhythm visible per channel (S2/S3 near-degenerate per the chain report).
- `variant_stage_scores.png` — the 10 stage-winner scores from `curriculum_{w2lvar,syn6}_stage{1..5}.json`
  (80.8/60.2/25.1/−225.9/−165.0 w2lvar; 95.6/86.8/37.5/−237.3/−197.2 syn6); caption states
  higher-is-better for stages 1–3, kine ≤ 0 (closer to 0 = better) for 4–5.
- `w2lvar_s5_walk.gif` — 110-frame offscreen render (10 fps, walk window 2–13 s, side view,
  camera follows COM) replaying the saved `qfull` qpos into the raw cvt3 model; CHEAP
  (~0.01 s/frame with mujoco.Renderer — no re-simulation), so it was made, not skipped.

Replay + regenerate (cwd `Code\MuJoCo_SNS\spinal`, env `myo`):
```
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" reports_20260925\figs\tools\replay_s5.py w2lvar   (then syn6)
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" reports_20260925\figs\tools\make_figs.py all
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" reports_20260925\figs\tools\make_gif.py
```
Evidence: `tmp\replay_{w2lvar,syn6}_s5.log` (PASS verdicts) + `tmp\figs_npz_{variant}.npz`
(the replay runs themselves). Helpers: `tools\replay_s5.py`, `tools\make_figs.py`, `tools\make_gif.py`.

---

# figs/ — goal-2 M3–M6 w2l_mujoco port figures (single-LH-RG / split-RG / afferented air + ground), 2026-09-25 night

The five AnimatLab→MuJoCo port deliverables. Every sim REPLAYED the milestone gate loop verbatim
and reproduced the gate numbers bit-exactly before any figure was drawn: M3 (period 1.027 s, hip L
41.8°, 22 632 self-contacts), M4 (1.357 s both sides, RG-E r −0.750, 26 741 contacts), M5 control
(1.389 s, heel SN 1.50 mV, t* = 6.244 s) + causal pulse (shifts +54/+12/+10/+12 ms), M6 (fall=no,
pelvis z min 0.660 m, tilt 11.7°, harness +39%, heel duty 0.00/0.00, heel SN 0.00 mV, 14 bursts @
1.357 s). Runner = `tools\w2l_runs.py` (imports the milestone machinery read-only; `data_*.npz` +
`run_*.log` in tools\ are the evidence), figure builder = `tools\w2l_figs.py`.

- `w2l_air_stepping.gif` — M3: the 2023-ORIGINAL single-LH-RG net air-stepping the axis-fixed M1
  body (`w2l_air.xml`, pelvis welded +0.30 m). 200 offscreen MuJoCo renders (0.1 s of sim time per
  frame, 10 fps, 20 s), left = skeleton side view, right = live traces (hip flexion L/R window +
  L RG-E window + full-run hip-L strip with cursor) + knob footer. Backed by
  `goal2_m3_w2l_cpg_air.md`.
- `w2l_air_joints_vs_reference.png` — M3 static: hip/knee/ankle flexion L+R over the full 20 s vs
  the AnimatLab references (per-panel sim range vs ref ~38/61/16°, gray band = ref range around
  the trace mid); title carries the cadence comparison (measured ~0.97 Hz vs 2023-original
  0.77 Hz within 2× / modern 2.22 Hz). Same M3 report.
- `w2l_split_rg.gif` — M4: split-RG coupled air walk (same layout; traces = hip L/R + RG-E L/R,
  footer states the c1 3.0 / V3 0.08 coupling). Shows the locked pair at 1.357 s (0.74 Hz),
  L/R phase 0.507 cycle, band r −0.941. Backed by `goal2_m4_rg_split.md`.
- `w2l_afferented_air.gif` — M5: afferented split-RG walk with Ben's contact rules, canonical
  flags (amp 1.5 nA scripted heel/toe, pamp 4.0, pulse 0.5 s, ibnA 1). Traces = heel SN L/R +
  RG-E L/R; the EXTRA causal heel-L pulse (4 nA × 0.5 s at t* = 6.24 s) is visible as the tall
  spike with a pink dashed marker, and the RG-E onset shift (+54 ms first onset, then +10..12 ms
  persistent) shows right after it — the gate-b causal reset, animated. Backed by
  `goal2_m5_afferents.md`.
- `w2l_ground_supported_walk.png` — M6: the honest "does not walk" figure — harness-supported
  (S=1.0) in-place march, pelvis FREE: joint angles at the autonomous 1.357 s RG period, pelvis/COM
  height + tilt (harness carries +39% weight, pelvis z min 0.660 m, tilt max 11.7°), and plate
  forces showing ZERO heel loading (duty 0.00/0.00, heel SN max 0.00 mV → the M5 heel stance-reset
  never engages; free stand topples at 0.89 s because the COM sits 3.09 cm outside the support
  polygon). Backed by `goal2_m6_stand_ground.md`.

Render notes (what rendered vs fell back):
- **mujoco 2.3.7 `mujoco.Renderer` offscreen WORKS on EB475WS4** (first try, no EGL/OSMesa
  needed; probe kept as `tools\_smoke_render.py`). All three air GIFs are true GL renders
  (360×470 side view, azimuth 90, brightness/contrast lifted 1.22/1.05 in PIL for visibility);
  NO matplotlib-skeleton fallback was needed and NO pip installs were made (PIL 12.3.0 already in
  the `myo` env assembled the GIFs; imageio is not installed and was not needed).
- Frames captured every 100 physics steps (0.1 s) → 200 frames per 20 s run; GIF duration 100 ms/
  frame (10 fps). Camera fixed at lookat (−3.40, 0, 1.05), distance 2.35 (M6 is a static figure,
  no frames). Raw frames were deleted after GIF assembly — regenerate with
  `tools\w2l_runs.py m3|m4|m5` then `tools\w2l_figs.py all` (cwd any; env `myo` full path).
- One pitfall hit and fixed (recorded for the next session): `test_w2l_ground` re-wraps
  `sys.stdout` with a new TextIOWrapper over the shared buffer on import; the released wrapper is
  gc-closed mid-run and kills stdout ("I/O operation on closed file"). `tools\w2l_runs.py` keeps
  strong references to every stdout wrapper (`_STDOUT_KEEP`) and saves npz before printing.

