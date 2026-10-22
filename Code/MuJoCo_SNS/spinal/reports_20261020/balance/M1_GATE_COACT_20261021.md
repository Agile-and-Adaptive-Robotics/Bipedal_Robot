# M1 GATE — co-contraction stand descent stream (2026-10-21)

Gatekeeper execution record per SUPERVISOR_AGENT.md at milestone M1
(result claim: "no uniform co-contraction level free-stands; measured
ankle PF:DF torque imbalance is the blocking mechanism"). DEVIATION
(same as M1_GATE_20261021.md and M1_GATE_S3Q_20261021.md): no subagent
spawn facility exists in this session's environment, so the checklist
was executed in-session with every verification command re-run
firsthand.

## Checklist

1. FIGURE STANDARDS — PASS (n/a). No figures produced this stream.
2. DISSERTATION TEX — PASS. No .tex touched by this stream. The
   working tree's modified Main.tex/Main.pdf/Main.log under
   Documentation\...\Knee_Torque_Test are Ben's own pre-existing
   uncommitted edits (noted in the s3q status entry), verified
   untouched. git status: only spinal\DESIGN.md,
   spinal\reports_20261020\SUPERVISOR_STATUS.md, spinal\runner.py
   modified by this session (spinal_run.npz overwritten by the anchor
   eval was git-restored per the AGENTS convention).
3. CONNECTOME OWNERSHIP — PASS. build_network.py and params.py have
   ZERO diff. The drive rides the existing 2026-10-16 per-side ALPHA
   hand-drive cells (Ben's overnight-campaign spec). The only code
   change is runner.py `--alpha-stand` (+41 lines: parse, pre-build
   setup, port-current override) — default-off, and the original
   `drive * G[alpha_s]` path is preserved verbatim in the else branch.
   No new neurons or synapses at defaults.
4. LITERATURE FIDELITY — PASS (n/a). No circuit figures drawn.
5. BEN'S STANDING INSTRUCTIONS — PASS. Trunk tilt reported every run
   (tilt_max in every JSON record); per-leg joint envelopes recorded
   per run (hip_flexion_r/l, knee_angle_r/l, ankle_angle_r/l in the
   joints dicts); units stated (N·m, deg, %BW). Note: antagonist-group
   activation means pool both legs — per-leg joint data preserved in
   the JSONs.
6. HONEST REPORTING — PASS. The verdict is negative and stated as
   negative; the per-level descent table is decomposed (activation
   level, act_sum, stands/falls); the feet-%BW numbers on fallen runs
   are explicitly labeled inverted-heap artifacts, not standing loads;
   the descent's nonexistence is stated, not softened.

## Verification commands re-run firsthand this session

- Anchor: `D:/Anaconda/envs/myo/python.exe _s3k_anchor_eval_20261016.py A`
  → -161.56754173676563 BIT-EXACT after the runner edit
  (`_anchor_coactstand_20261021.log`).
- Census: `D:/Anaconda/envs/myo/python.exe _s3p_census_20261020.py`
  → defaults (410, 376, 1186) GATE PASS; coact=0.8 adds exactly 4
  synapses (1190) and restores (1186) — default build untouched.
- Fresh-process reproduction of headline rungs:
  `one stock 2.0 5.0 0.0 8` → fell, act_sum 92.0, sway 0.903,
  COM -0.641, feet 104.3 (identical to the in-process matrix rung);
  `one stock 2.0 2.0 0.0 8` → fell, act_sum 83.2, feet 104.8,
  COM -0.66 (identical to the 40% matrix rung). One-command
  verification of the central claim stands.
- Static torque probe re-run with corrected moment-matrix indexing:
  ankle_r net -182.7 N·m (PF -236.0 vs DF +53.3), knee balanced
  (-8.3/+11.9), hips -152.9/-17.8.

## VERDICT: SHIP

Advisories (non-blocking):
- The stock matrix JSON's buckle timestamps use the earlier t=2.0-s
  reference convention (the session upgraded forensics mid-flight);
  the fresh `one` run and the syn6/w2lvar JSONs use the t=1.0-s
  reference — the table cites the consistent fresh numbers.
- Per-run npz scratch (`_coactstand_*_20261021.npz`, ~14 files) is
  reproducible from the driver; kept as per-rung evidence for Ben to
  commit or delete.
- The W2L/BPA robot stacks have no ALPHA/coact mechanism (stated as
  N/A, not dodged — their unrigged capability records already exist).
