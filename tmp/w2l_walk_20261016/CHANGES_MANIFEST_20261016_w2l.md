# CHANGES MANIFEST — W2L robot audit + L/R drive (2026-10-16, egg)

Working tree only. NO commits, NO pushes (standing contract). Recovery =
`git checkout -- <file>` (all files below are tracked and were clean at
HEAD b497b67b before this session).

## Modified (7 tracked files, all under Code\MuJoCo_SNS\spinal\)

1. `w2l_cpg\build_w2l_net.py`
   - W2LNet dataclass: `alpha_l` / `alpha_r` fields (default 0.0).
   - Conditional ALPHA input-port block after the TONIC block (ports
     "ALPHA <MN label>", 6 per side, built only when that side's alpha
     != 0 — default builds nothing).
   - `_assert_complete`: input-count assert 10 -> 10 + len(alpha_ports).
   - `build(..., alpha_l=0.0, alpha_r=0.0)`; docstring notes gamma = the
     TONIC ports.
2. `w2l_cpg\smoke_w2l.py`
   - New knobs `--tonic-e-l/-r`, `--tonic-f-l/-r` (per-side gamma,
     default None = symmetric baseline), `--alpha-l/-r` (passed to
     build). Per-side tonic/alpha injection in the loop. Usage docstring
     updated.
3. `w2l_cpg\README.md`
   - Knobs paragraph: documents the new per-side drive knobs + the
     bit-identical default contract.
4. `w2l_mujoco\build_w2l_split_net.py`
   - W2LSplitNet: alpha_l/alpha_r fields + the same conditional ALPHA
     port block (inherited by W2LAffNet); input-count assert 6 -> 6 +
     len(alpha_ports); build() passthrough + docstring.
5. `w2l_mujoco\build_w2l_aff_net.py`
   - build(..., alpha_l, alpha_r) passthrough; aff census assert n_in
     +len(alpha_ports); docstring census corrected 188 -> 186 synapses.
6. `w2l_mujoco\test_w2l_air.py`
   - `write_air_xml()` REGRESSION FIX: the air variant now WELDS the Root
     (strips `<freejoint name="root"/>` and the `<keyframe>` block from
     the GENERATED w2l_air.xml only). Reason: the 2026-10-01 rebuild
     freejoint-ed the Root; lift-alone let the air body free-fall
     (77,646 ground contacts, gates M3/M4/M5 all broken at HEAD).
     w2l_mjcf_fixed.xml and the M6 ground path are untouched.
7. `w2l_mujoco\test_w2l_air_afferented.py` and
   `w2l_mujoco\test_w2l_ground.py`
   - KV knobs `te_l/te_r/tf_l/tf_r` (nan = symmetric baseline) +
     `alpha_l/alpha_r` (0 = ports not built); per-side tonic + alpha
     injection in run_sim/run_phase; afferented header prints the
     per-side values when any is set.

## New (untracked, campaign outputs)

- `reports_20261016\w2l_audit\` — this manifest's sibling logs/scripts:
  W2L_ROBOT_AUDIT_20261016.md (the report; copied to
  reports_20261016\ too), _w2l_rules_audit.py, _w2l_drive_gate.py,
  _w2l_alpha_debug.py, m5_gate_baseline.log, m5_gate_postedit.log,
  m5_gate_canonical.log, m6_gate_baseline.log, m6_walk_postedit.log,
  m6_walk_gammaL_demo.log, m6_walk_alphaL_demo.log, m4_air_postedit.log,
  m3_air_postedit.log, drive_gate_preedit.log, drive_gate_postedit.log,
  circuit-rules copy NOT left here (read from the repo paths).

## Not mine (attributed)

- `spinal\spinal_run.npz`, `spinal\build_network.py`,
  `spinal\build_network_syn6.py`, `spinal\build_network_w2lvar.py`,
  `spinal\params.py`, `spinal\runner.py`, untracked `spinal\_audit_human/
  _dup_probe/_lrdrive_probe/_s3k_anchor_eval_20261016.py` — the SIBLING
  campaign stage's L/R-drive work on the human/spinal side (no command
  in this audit touched them).
- `w2l_mujoco\w2l_air.xml` — a TRACKED but GENERATED artifact (written
  by every air-gate run); modified by this session's gate runs (now the
  welded variant). Harmless byproduct; regenerated on each run.
- `Solid_Models\...05_01_TI_R_008.SLDDRW` + untracked 8020/tmp items —
  pre-existing working-tree noise (per the session recon).

## Verified post-change (commands + results in the report §3)

- Default census: (95,208,10) / (87,164,6) / (97,186,12) unchanged.
- Bit-exact: probe hashes + smoke output + M6 walk log identical to
  pre-edit baselines.
- Functional: alpha/gamma knobs move exactly their own side's MN pools /
  the closed-loop walk; demonstrations logged.
