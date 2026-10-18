# CHANGES MANIFEST — bootstrap driver (robot) stage, 2026-10-16

Campaign stage 7 ("bootstrap driver (robot)"). Working-tree edits on
easteregg2 (`D:\GitHub\Bipedal_Robot`), NO commits, NO pushes. All other
modified files in `git status` belong to the earlier audit stage
(see `reports_20261016\w2l_audit\CHANGES_MANIFEST_20261016_w2l.md`).

## Code\MuJoCo_SNS\spinal\w2l_mujoco\test_w2l_ground.py

1. **alpha i-shadow BUG FIX** (run_phase, neural-update block):
   `for i, a_cur in i_alpha:` → `for i_a, a_cur in i_alpha:` with
   comment. The audit-stage loop shadowed the physics loop index `i`;
   with alpha ports built it corrupted every even-index log entry, fed
   the fall check stale tilt, and disabled the rig anchor capture
   (leash/tilt never engaged; body wandered 11 m off-platform).
2. **KV knob `legcoll=1.0` (default = shipped behavior)** + remask block
   after `MjModel.from_xml_path(GROUND)`: `legcoll=0` sets L-subtree
   geoms contype/conaffinity 1, R-subtree 2, Root+ground 3 — disables
   ONLY L<->R leg collisions (the planar aproj has none; the 3D MJCF
   tibias interpenetrate at spawn with 3368 N measured).
3. **KV knob `walker_yield=0.0` (default = fixed anchor)** + EMA block at
   the top of the leash section: anchor x follows the COM with time
   constant tau — the MOVING WalkingPath analog.

Gates re-run post-edit (default knobs): default walk run reproduces the
recorded baseline EXACTLY — finite, no fall, pelvis z min 0.951, tilt
8.2, harness 49 %, strikes 0/0, duty 0/0, fwd +0.001, L/R RG-E bursts
11/10 period 1.539 s, r −0.776, verdict FALLS/NO-GAIT (verified after
EACH of the three edits; logs in this folder's session transcript, run
live over ssh).

## Code\MuJoCo_SNS\spinal\w2l_mujoco\test_w2l_air_afferented.py

4. Same i-shadow fix (line ~166) — alpha loop renamed `i_a`. M5 canonical
   gate re-run post-fix: gate_a PASS, gate_b PASS, gate_c FAIL-as-
   specified (suppression nan) — identical to the audit's recorded
   state, so the no-alpha path is untouched.

## New files (all under reports_20261016\w2l_walk_sweep\ unless noted)

- `_w2l_drive_sweep.py` — sweep driver (subprocess + log-parse → CSV)
- `_w2l_probe_run.py` — instrumented run_phase probe (raw-log analysis,
  strike times, phases, npz dump; `--legcoll=0` monkeypatch predates the
  knob)
- `_w2l_static_pose.py` — spawn-pose geometry probe (found the tibia
  interpenetration)
- `_inspect_npz.py`, `_duty_calc.py`, `_write_best_params.py` — one-off
  analysis helpers (kept for provenance)
- `combos_{a,bc,d,e,f,g,h,w}.json` + `results.csv` + `logs\*.log` —
  128-run sweep record
- `probe_*.npz` — raw logs of the key runs (c05 pre-fix corruption
  evidence, d03, f05, f18, h07, w01 winner, w05)
- `w2l_walk_best_params.json` — WINNER param file (generated from
  probe_w01.npz, not hand-typed)
- `..\W2L_ROBOT_BOOTSWEEP_20261016.md` — the stage report
- this manifest
