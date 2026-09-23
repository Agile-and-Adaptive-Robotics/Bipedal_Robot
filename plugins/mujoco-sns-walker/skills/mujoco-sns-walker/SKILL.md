---
name: mujoco-sns-walker
description: Standing knowledge for MuJoCo + SNS-Toolbox spinal walker sessions (Code\MuJoCo_SNS\spinal) - env bootstrap, runner/curriculum contracts, tuning-study hygiene, Windows traps, vision/literature toolkit. Use whenever working on the gait2392 spinal-cord walker, its optuna studies, connectome spec, or circuit figures.
version: 0.1.0
---

# MuJoCo-SNS walker — session knowledge

## Environment bootstrap (every command)

- The SNS env python, by machine (check `hostname` FIRST — the AGENTS.md machine guard):
  EB475WS4 `C:\Users\Ben Bolen\.conda\envs\myo\python.exe` · easteregg2
  `D:\Anaconda\envs\myo\python.exe` · laptop `C:\Users\Ben\.anaconda3\envs\myoconv\python.exe`.
  The plugin's `scripts\walker.cmd` resolves all of this; `AARL_PYTHON` overrides.
- Set `CONDA_PREFIX` to the env dir (MuJoCo reads it). `walker.cmd` does.
- Python stdout on this Windows console: utf-8 wrapper needed or prints die on any
  non-ASCII (`wcommon.utf8_stdio()`).
- Bare `python` on PATH is the Microsoft Store stub — always full paths.
- Work dir: `Code\MuJoCo_SNS\spinal\` (the optuna db URL is RELATIVE: `sqlite:///optuna_walk.db`).

## Runner contract (runner.py, hand-parsed flags — no argparse)

- Flags: `--fitted`, `--best`/`--best5`…`--best10` (winner jsons), `--eval` (16 s schedule,
  quiet, RETURNS a metrics dict — prints nothing), `--drive N`, `--joint-pf [X]`,
  `--harness K|--no-harness`, `--no-afferents`, `--no-ground`, `--no-interleg`,
  `--leg-damping X`, `--rig-scale S`, `--phase-reset E F`, `--kneext-inh X`, `--renshaw X`,
  `--straight-start`, `--view`, `--scope`, `--realtime`, `--time S`. Unknown flags are
  SILENTLY IGNORED — typo'd flags don't error.
- Env vars: `AARL_MODEL`, `AARL_POSE` (normal|symmetric), `AARL_PELVIS_TY`, `AARL_KY`
  (lateral rig spring scale ×5e5), `AARL_DEAFF` (l|r|l,r), `AARL_NPZ` (output name),
  `RUNNER_DUMP_STATE` (config-diff json).
- npz layout: `t`, `q` = **13 joint angles in DEGREES** ordered by `key_joints`
  (pelvis_tilt, pelvis_list, pelvis_rotation, hip_flexion_r, knee_angle_r, ankle_angle_r,
  hip_adduction_r, subtalar_angle_r, then _l), `qfull` (rad), `act`, `com`, `neuro`
  (`neuro_names`), `contact` (col 0 = right), `footz`, `cfg`. Atomic write with
  os.replace retry — close NpzFile handles or WinError-5.
- Knee convention: NEGATIVE = flexion (converter preserved OpenSim's flexion-negative knee).
- Stdout grammar: startup lines (`loaded …`, `connectome_gains.json applied (N rules,
  full_rules=…)`, `start pose …`, `standing solve …`), 0.5 s progress lines, final
  `== summary ==` block (stayed up/FELL, per-leg `RG_r: cycles … E-duty …`, joint ranges).

## Curriculum / studies

- `python _curriculum.py <stage> [n]`; stages 1 air-deaff, 2 air-aff, 3 ground.
  Stage-3 study name lives in source (currently `curr_s3k_nocross`) — renaming for a new
  round is a source edit (Ben's wiring).
- Objective metric keys: `kine_score` (higher/less-negative better; kine_ref.compare —
  includes bilateral bonus +40 / frozen cost +30), `kz`, `tilt_max`, `duty`, `knee_min`.
- **JSON RULE**: every params.G knob a study uses must be saved in the winner json AND
  have a loader branch in runner.py (`/knob-check` audits the four places).
- Seeding: `study_manage.py seed <study> <trial>` rewrites curriculum_stage3.json;
  missing seed keys get SAMPLED by optuna (stage-2 lesson).
- Sentinels: no-rhythm/NaN/frozen configs must sit BELOW the worst genuine walker —
  recalibrate when the plant/objective/pose changes.
- Regression gate: winner reproductions are bit-exact only with identical binaries +
  seeds; late-run divergence between integrators is PHASE, not signal. The network is
  near a bifurcation — marginal winners may not transfer (basin gate: `basin_gate.py`,
  PASS = rhythm persists under ±1 % param perturbation).
- Known v10-winner baseline numbers (for quick sanity, NOT targets): e2_pf 0.222,
  f1_df 0.098, f1_kf 0.148, post_kneext 0.011, post_hipext 0.068, desc_f 1.400,
  e2_adapt 1.700, pelvis ty 0.909 m.

## Connectome ownership (hard rule)

- **Ben owns the wiring; the agent only tunes.** Wiring changes must trace to his spec
  (`connectome_gains.json` from connectome_editor.html) or a cited literature rule.
- runner consumes ONLY top-level `{"rules": {id: {enabled, gain_key, gain, hops}}}`;
  gain_key must exist in params.G or the rule is silently skipped; disabled rules are
  skipped-not-zeroed; any enabled hops≥1 sets global `full_rules`. `/connectome-check`
  before any spec-honoring run. Block-editor exports are design docs, NOT consumed.
- Supervisor gates M1/M2/M3 (`/supervisor-gate`) are mandatory before declaring results,
  editing figures/tex/circuit, or claiming literature fidelity.

## Windows traps (each has cost a session)

- Compound cmd chains (`cd && set && …`) silently fail — never chain; `walker.cmd` +
  `run_gate.py go` avoid shell entirely.
- Long waits: background tasks + poll logs (`== stage`, `DONE`, `FAILED`), never spin.
- File locks on `spinal_run.npz` (Spyder indexer/AV): unique `AARL_NPZ` names; atomic
  writes; close npz handles.
- After stopping runs: check for orphan python/bat processes (`status.py` lists them).
- `py_compile` the four core files before every launch (run_gate compile).

## Vision & literature toolkit (MCP tools)

- `mcp__zai-mcp-server__understand_technical_diagram` — circuit/architecture figures
  (paper connectomes, our draw_circuit outputs). `analyze_image` — general/markups.
  `ui_diff_check` — expected-vs-actual figure comparison. `analyze_data_visualization` —
  charts/plots. `extract_text_from_screenshot` — OCR of code/terminal screenshots.
  Images only (PNG/JPG) — render PDF pages first (MiKTeX `mgs.exe` on EB475WS4,
  ghostscript in the gs env on easteregg2).
- `/circuit-vision` is the guided workflow (paper figure ↔ our circuit diff, figure
  audit, Ben-markup transcription). Cross-check every vision claim against
  `_net_edges.py` compiled inventory.
- `web-search-prime` / `web-reader` — literature lookup and doc reading (other
  neuromechanical modeling software docs, paper pages). Full texts locally:
  `spinal\lit_*.txt`, `*_fulltext.txt`, `spinal\lit_pdfs\`, Zotero via
  `D:\Github\api_credentials_local.txt`.
- Figure standards are PROJECT-WIDE (AGENTS.md): 7.5×10 in usable, Arial ≥10 pt, no
  italics, Paul Tol palette from `Code\Matlab\Colors.m`, shape-coded synapses, alt-text
  file beside every figure.

## Plugin tools index

`/walker-status` `/walker-run` `/walker-probe` `/study-manage` `/connectome-check`
`/connectome-gui` `/circuit-vision` `/supervisor-gate` `/knob-check` — all go through
`plugins\mujoco-sns-walker\scripts\walker.cmd <script.py> [args]`.
