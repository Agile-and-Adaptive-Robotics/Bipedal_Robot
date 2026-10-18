# CHANGES manifest — 2026-10-16 connectome audit + L/R alpha/gamma drive

Stage: connectome auditor (human), overnight walk campaign
(dwfrun-5685cc77). NO commits, NO pushes. The five files below are
modified in the easteregg2 working tree (D:\GitHub\Bipedal_Robot,
branch KneeTestSetup_BenBo_stw) AND identically in the laptop authoring
tree (C:\Users\Ben\Documents\GitHub\Bipedal_Robot) — same content, both
uncommitted. Revert with
`git checkout -- Code/MuJoCo_SNS/spinal/{params.py,build_network.py,build_network_w2lvar.py,build_network_syn6.py,runner.py}`
(either machine).

## Modified files (git diff --stat: +133 insertions, 0 deletions)

1. `Code/MuJoCo_SNS/spinal/params.py` — five new params.G entries at the
   end of G: `gamma_l=0.0, gamma_r=0.0, alpha_l=0.0, alpha_r=0.0`
   (the four per-model hand-drive knobs; each > 0 builds ONLY that
   side's drive cells) and `alpha_to_mn=0.5` (conductance of the
   ALPHA_s→MN edges; edges exist only when alpha_s > 0).
2. `Code/MuJoCo_SNS/spinal/build_network.py` —
   a. `__post_init__` (after the VEST block): per-side
      `self.gamma_side` / `self.alpha_side` flags + creation of
      GAMMA_s / ALPHA_s neurons and GAMMA_c_s / ALPHA_c_s input ports,
      each only when that side's knob > 0.
   b. `_build_rg` (after the POSTURE→RG_E edge): GAMMA_side → RG_E_side
      (descend_to_rg_e) and → RG_F_side (descend_to_rg_f) when
      gamma_side.
   c. `_wire_balance` (end of the muscle loop): ALPHA_side → every
      ipsilateral MN_{act} at alpha_to_mn when alpha_side (inherited by
      the w2lvar/syn6 subclasses).
3. `Code/MuJoCo_SNS/spinal/build_network_w2lvar.py` — same per-side
   flags in its `__post_init__`, same cell creation after its
   descending-cells loop, GAMMA wiring after the DRIVE→rg_e/f edges in
   `_build_rg_w2l`. ALPHA fan-out inherited from the parent.
4. `Code/MuJoCo_SNS/spinal/build_network_syn6.py` — same pattern (flags
   + cells after the brainstem block; GAMMA wiring after the DRIVE
   edges in its `_build_rg`). ALPHA fan-out inherited.
5. `Code/MuJoCo_SNS/spinal/runner.py` — step loop, immediately after
   `u[iport["POSTURE"]] = posture`: port-name-guarded feed
   `u[GAMMA_c_s] = drive * G["gamma_s"]`, `u[ALPHA_c_s] = drive *
   G["alpha_s"]` for sides whose port exists. At defaults no ports
   exist → the u vector is untouched → bit-identical. NO other runner
   change (no new CLI flags; the knobs are params.G/json-driven per the
   campaign's four-knob spec).

## Semantics (for phase 4)

- gamma_s (RG/rhythm level): side s's RG pair receives
  drive·(1+gamma_s) through the winner-tuned descend_to_rg_e/f
  conductances (its own GAMMA cell in parallel with the shared DRIVE).
- alpha_s (tonic MN-pool level): every ipsilateral MN pool receives
  drive·alpha_s through a uniform alpha_to_mn synapse from that side's
  ALPHA cell (per-side analog of the shared POSTURE cell).
- JSON RULE: studies must save the four keys in winner jsons; loader =
  plain `P.G[k] = float(p.get(k, 0.0))` before build (set_stage-style
  branch). Wean ladder = shrink the gamma values toward 0 in the json.

## Gates run (egg, post-change) — all PASS

- `reports_20260925/audit/audit_defaults_gate.py` → (410, 376, 1186) PASS
- AARL_NET=w2lvar default build → (888, 382, 7866) (== pre-change)
- AARL_NET=syn6 default build → (794, 382, 3033) (== pre-change)
- `_s3k_anchor_eval_20261016.py A` → -161.56754173676563 BIT-EXACT
- `_s3k_anchor_eval_20261016.py B` → -25.0 (bit-exact vs its own
  pre-change value; flavor B = production full_rules=1 topology)
- `_lrdrive_probe_20261016.py` → winner+knobs (882, 386, 8492); GAMMA
  edges at descend_to_rg_e/f winner values; 46 ALPHA→MN per side.

## New untracked files (evidence + reusable gates; safe to delete)

spinal\_audit_human_20261016.py, spinal\_dup_probe_20261016.py,
spinal\_s3k_anchor_eval_20261016.py, spinal\_lrdrive_probe_20261016.py,
spinal\reports_20261016\{audit_human_edges_20261016.txt,
anchor_pre_change.log, anchor_post_change_A.log,
anchor_post_change_B.log, lrdrive_probe_20261016.txt,
connectome_audit_human_20261016.md, CHANGES_20261016_connectome_drive.md}.

## NOT changed

- No existing synapse, neuron, or gain anywhere (the audit's findings —
  stacked IBEXC→LBIN duplicates (N1), group-size-weighted LBIN drive,
  the rulings-doc undercount (N2), toe/heel routing gaps vs the drawing
  — are REPORTED in connectome_audit_human_20261016.md, not fixed;
  fixing N1 changes the tuned winner's dynamics and is Ben's call).
- The sibling stage's concurrent egg edits (spinal_run.npz,
  w2l_cpg/*, w2l_mujoco/*) were not touched.
