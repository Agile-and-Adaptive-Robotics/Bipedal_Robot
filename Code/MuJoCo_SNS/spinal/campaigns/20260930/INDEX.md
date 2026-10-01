# Campaign 20260930 — sustained modeling push (Ben's "everything" directive)

Index of everything this campaign produced or set running. Raw logs live on
easteregg2 unless noted; copies land here as jobs finish. Companion artifacts
from the 2026-09-28 prune matrix (the foundation for this campaign) are in
`spinal\` root: `prune_analysis.md`, `prune_results_*.jsonl`, `prune_matrix.py`,
`prune_combo.py`, `analyze_prune.py`.

## Running on easteregg2 (WMI-detached, survive disconnects)

| Job | Launched | Artifact | Status at last contact |
|---|---|---|---|
| Pruned-s3k air retune (30 trials) | ~00:00 | `optuna_pruned_pruned_air.db`, `curriculum_pruned_air.json` | running (after sqlite-lock fix) |
| Pruned-s3k walk retune ×3 parallel (batches of 40, stop rule 2×<+0.5, cap 5 batches) | ~00:00 | `optuna_pruned_pruned_walk_{1,2,3}.db` | running |
| Kinematic robustness sweep (4 variants × pelvis-ty ±1/±2 cm × WALK/STAND/PUSH) | ~00:00 | `robust_results_*.jsonl` | w2lvar COMPLETE (20 cells) before VPN drop; s3k/s3kpruned partial; syn6 relaunched after TY_BASE fix |
| Literature validation (4 walkers × 43 gait_refs) | ~00:00 | `gait_validation_20260930.csv` | s3k scored; s3kpruned in flight |
| Li chain: rebuild MJCF (freejoint + ΣB·r² joint damping + spawn keyframe) → axis fix → validate → 20 s gate → W2L air smoke | ~00:15 | `w2l_mujoco/li_chain.log`, regenerated `w2l_mjcf.xml` | running (group-redirect fix) |
| SNS Simscape (units test + KneeReflexDemo_R2025a) | ~00:20 | `SNS_Simscape/logs/easteregg2_runs_20260930/runs2.log` | running |

**VPN incident ~00:30:** the CECS OpenVPN tunnel dropped mid-campaign (adapter
gone, all lab hosts unreachable from the laptop). The OpenVPN Connect client
was driven to RETRY; it now waits at an **"Enter password" prompt for the CECS
split-tunnel profile** — Ben must type his PSU password in the client window
(dialog open + focused on the laptop screen). All easteregg2 jobs above are
detached and unaffected. First action after reconnect: run
`powershell -File D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\campaign_status.ps1`.

## Completed (local)

- **`lit_rules_conflict_audit.md`** — cross-comparison of rybak/shevtsova/
  shinohara rules vs Ben's own rules. 8 sign conflicts (the real ones:
  RG-E self-connection recurrent-exc vs lamination-inh; RG→PF gate sign;
  Ben's dual Ib/II/RC pathways both directions in one rule set), 6 gain
  disagreements >2× (supraspinal drive 1.0 vs 0.02–0.15; RG→PF 0.0075 vs 0.7).
- **`figs/prune_walk_deltas.png`, `figs/prune_air_scores.png`,
  `figs/prune_push_recovery.png`** — the 09-28 prune matrix as publication-grade
  bar charts (generator: `prune_figs.py`).
- **Connectome editor**: new templates `walker_s3k_pruned` (subsystem schematic,
  13/292) + `walker_s3k_pruned_flat` (233/390 vs s3k's 420 edges — the 30
  prunable-tagged edges removed). Dropdown labeled "s3k PRUNED". Generator now
  repo-relative (runs on any machine; old D:\Github hardcode fixed). Static
  checks pass; browser/js suite pending (no node on the laptop).
- **AnimatLab headless on easteregg2**: 3 asims ran clean into
  `campaigns/20260930/animatlab/` (BilateralRG Ground, W2L modern, Biped
  Standalone) with chart traces (Rhythm Generator.txt etc.) captured for
  figures.
- **W2L split-net V3 correction**: invented V3→contra-RG-E HC edges removed
  from `build_w2l_split_net.py` (09-26 audit ruling; Ben's 09-30 go) — now
  V3→contra RG ext IN only, matching the real BilateralRG.aproj.

## Code changes (uncommitted, laptop + easteregg2 both)

- `runner.py` — `--push/--push-time/--push-axis` (default-off; 09-28)
- `prune_matrix.py`, `prune_combo.py`, `analyze_prune.py`, `prune_robust.py`,
  `prune_figs.py`, `_curriculum_pruned.py`, `gait_validate_all.py`,
  `lit_rules_conflicts.py` (all new)
- `w2l_mujoco/make_w2l_mjcf.py` (freejoint + per-joint ΣB·r² damping + spawn
  keyframe + env-overridable aproj path), `w2l_mujoco/validate_body.py`
  (9-joint gate), `w2l_mujoco/test_li_stepping.py` (keyframe spawn + XML
  damping), `w2l_mujoco/build_w2l_split_net.py` (V3 fix)
- `make_editor_templates.py` (portable paths + pruned templates),
  `connectome_block_editor.html` (2 dropdown entries), `connectome_templates.json`
- Laptop env change: `optuna 5.0.0` pip-installed into myoconv (schema-parity
  with easteregg2's myo env)

## Next when compute returns

1. `campaign_status.ps1` → harvest all six job families.
2. `_curriculum_pruned.py finalize` (easteregg2) → winner + validation cells.
3. GIFs: `_render_gif.py` on the fresh `spinal_run_val_*.npz` per variant
   (4 walk gifs + combo gif into `figs/`).
4. Simscape `SNS_SpinalNetwork` verify needs an R2025b export (easteregg2 is
   R2025a) — export from the laptop R2025b via `export_slx_to_R2025a.m`
   inverse, or run the verify on the laptop.
5. Editor browser-verification pass (node machine).
