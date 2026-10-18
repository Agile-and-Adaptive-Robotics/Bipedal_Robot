# W2L ROBOT SIDE — RULES AUDIT + L/R ALPHA/GAMMA DRIVE (2026-10-16, egg)

Overnight-campaign stage: "connectome auditor (robot)". All work ran ON
easteregg2 (repo `D:\GitHub\Bipedal_Robot`, branch KneeTestSetup_BenBo_stw,
HEAD b497b67b) via ssh; python `D:\Anaconda\envs\myo\python.exe`.
Logs live in `spinal\reports_20261016\w2l_audit\`. NO commits, NO pushes —
working tree + this report + the manifest only.

Audited:
- `spinal\w2l_cpg\` — build_w2l_net.py (bilateralrg contact-driven CPG)
- `spinal\w2l_mujoco\` — M4 split-RG + M5 Ben's-rules afferents on
  w2l_mjcf_fixed.xml (12-muscle MJCF robot)

Rules reference (same as the sibling audits): `spinal\ben_rules_20260924.json`
== `Neuromechanical_Models\Mujoco_SNS_models\
Circuit_rules_CONNECTOME_md__connectome.json` — payload byte-identical
(sha256-prefix daf0017baeb6e56a; verified by `_w2l_rules_audit.py` §1).

---

## 1. RULES VERDICT (details + evidence per item below)

| rule | verdict |
|---|---|
| RG+PF feedback PER PF LAYER | **Compliant in kind, ankle layer absent** — the 2023 template has hip+knee PF layers only (ankle MNs ride the knee PFs, `pf_to_mn`); heel and Ib land per-layer on the layers that exist. Documented reduction, not silent. |
| HEEL = stance-phase reset AT THE PF LAYER | **VERIFIED** — 4 of the drawing's 5 edges per side, all exc g 0.5: heel→ipsi RG-ext-IN (InE), heel→contra RG-flx-IN (InF), heel→ipsi Hip/Knee PF-ext-IN. The 5th (ankle PF-E IN) is N/A (no ankle PF layer). Gain key `aff_heel`=0.5 = the rules' number. |
| TOE = dorsiflexion inhibition ONLY | **VERIFIED** — toe→toe IN (exc 1)→TOEDF IN (exc 5)→inh 2.749 onto the dorsiflexion drive. Terminal target reduced from the drawing's HC-PF-Dorsiflexion to the DF MN pool (`{S} Ank MN flx`) because the 2023 net has no DF PF half-center (documented in build_w2l_aff_net.py:31-37). No other toe-driven edge exists anywhere in the net. |
| Ib load (PF-E 0.5, RG-E) | **VERIFIED with one extra edge** — see FINDING 3: Ben's drawing gives `ib_load` THREE out-edges (RG-E HC, InE, ONE PF-E HC); the implementation adds a FOURTH per side (Knee PF-E HC) that is not in the JSON. |
| Autogenic afferents | **Structure present, silent** — the template's per-pool autogenic chains (SN-Ib→homonymous MN `ib_auto` 12 edges; SN-Ia→IaIN→antagonist `ia_recip` 24 edges) are built and asserted but receive ZERO input in every current runner. The ACTIVE Ib signal is a per-side extensor-GROUP force (matches Ben's own "Ib grp" group node); no per-muscle autogenic afferent is active in this robot. |
| No UNINTENDED afferent→RG edges | **VERIFIED** — full enumeration (`_w2l_rules_audit.py` §6): the only afferent-source edges onto RG-family cells are the rules-intended heel→InE/InF and Ib grp→RG-E(HC)/InE, plus (in w2l_cpg only) the real aproj's `aff_HC_excite` family (24 edges, g 0.08, hip/knee, from Ben's BilateralRG.aproj per the 2026-09-26 correction). The 2023-original base (w2laproj) has ZERO afferent→HC edges — matching the Deng-audit doctrine. |

Ben's rules JSON also carries flexor Ia/II→{RG-F HC, InF, PF-F HC} exc-0.5
edges (`SN-Ia_2843`, `SN-II_2844`) and per-muscle motifs
(`ia_homo` monosynaptic exc 2, `ib_auto` via IbIN **inhibitory** 0.5,
`ii_exc`/`ii_inh`) — none of these are implemented in w2l_mujoco (the M5
scope was heel/toe/Ib only); the template's own `ib_auto` is excitatory
(Ben's aproj wiring), opposite sign to the rules' ib_auto motif. These are
GAPS (not wrong edges) — listed for Ben's next-rules pass.

## 2. FINDINGS (all verified this session on egg)

1. **[HIGH — FIXED this session] M5/M3/M4 air gates were broken at HEAD.**
   The 2026-10-01 body rebuild (commit 935ac739) added
   `<freejoint name="root"/>` to w2l_mjcf_fixed.xml; `write_air_xml()`
   (test_w2l_air.py) only lifted the Root pos, so the "air" body
   free-fell onto the ground: measured 77,646 ground contacts (must be
   0), joints pinned at 0.0° range, RG latched in-phase (r +0.681), and
   the causal gates crashed on `onA[2]` (no onsets)
   [m5_gate_baseline.log]. FIX: the AIR variant now welds the Root
   (freejoint + keyframe stripped IN THE GENERATED w2l_air.xml ONLY;
   w2l_mjcf_fixed.xml and the M6 ground path untouched). Post-fix:
   M5 gate (a) all-11 PASS at the canonical knobs [m5_gate_canonical.log],
   gate (b) causal heel reset PASS; M4 split air all-8 PASS
   [m4_air_postedit.log]; M3 orig air all-6 PASS [m3_air_postedit.log].
2. **[MEDIUM] M5 file defaults do not pass gate (a) — and never did.**
   `test_w2l_air_afferented.py` defaults (amp=4.0, ibnA=6.0) saturate the
   Ib group (Ib grp max 9.62 mV) whose direct Ib→RG-E excitation latches
   both RGs; the 2026-09-25 PASS record used CLI overrides
   `--amp=1.5 --pulse=0.5 --ibnA=1.0` (log knobs line), and the 09-25
   smoke at the file defaults itself crashed
   (reports_20260925/logs/m5_smoke1.log). Canonical invocation to
   reproduce: `test_w2l_air_afferented.py --amp=1.5 --pulse=0.5
   --ibnA=1.0 --pamp=4.0 --causal=both`.
3. **[MEDIUM] Ib load family carries one edge per side beyond Ben's
   drawing.** The rules JSON's `ib_load: Ib grp` node has exactly THREE
   out-edges (→RG-E HC, →IN-InE_2803, →HC-PF-E_2805, all exc 0.5); the
   implemented IB_EDGES add a fourth per side
   (`Ib grp {S} → {S} Knee PF ext`, tag `ib_load`) — and the
   build_w2l_aff_net.py docstring (lines 39-44) MISQUOTES the drawing as
   containing both Hip and Knee PF-E edges. Defensible under the
   "feedback is per PF layer" ruling, but it is not in the drawing —
   Ben's call: keep (spirit) or trim to 3/side (letter). Evidence:
   `_w2l_rules_audit.py` §2-4 dump vs built net §5-8.
4. **[LOW] Toe chain misses the rules' second TOEDF-IN input** — the
   drawing also has HC-PF-E_2809 → IN-PF_dorsiflexion_inhibit (exc 2.749)
   (the ankle-PF-E drive of the DF inhibitor); N/A here for the same
   no-ankle-PF-layer reason as the terminal reduction (same reduction
   family, this facet was not separately documented in the docstring).
5. **[LOW] M5 gate (c) FAIL-as-specified reproduced** — the DF pool never
   reaches positive drive in air (−2.625 → −2.625 mV over the pulse
   window, suppression nan), identical signature to the 09-25 record.
   Known structural cause (no DF PF layer / DF MN idle in air).
6. **[INFO] Muscle-map TODO STILL OPEN** — w2l_cpg README "What is NOT
   done: No connection to the 92-muscle MuJoCo model" stands at HEAD.
   Filesystem scan: the only spinal-root consumers of w2l_cpg are
   build_network_w2lvar.py and make_editor_templates.py (neither wires
   the 12 MN pools into gait2392); w2l_mujoco drives its own 12 template
   muscles (MUSCLE_MAP incl. the R-hip name crossing) only.
7. **[INFO] Ask-premise drift**: the port is now **95 neurons / 208
   synapses / 10 inputs** (the 91/214 figure predates the 2026-09-26
   template correction). Smoke baseline REPRODUCED bit-exact on egg:
   `W2L_SMOKE PASS period=1.000 antiphase_r=-0.525 rg_e_bursts_L=22
   rg_e_bursts_R=22` (pre- and post-drive-edit).
8. **[LOW, fixed] build_w2l_aff_net.py docstring census said 188
   synapses**; the asserted/built count is 186 (164+8+6+8). Docstring
   corrected.
9. **[INFO] `spinal_run.npz` was modified on egg during this session by
   another process** (sibling campaign stage; none of this audit's
   commands write it). Recorded so the working-tree diff is attributable.

## 3. L/R ALPHA (MN) + GAMMA (RG) DRIVE — implemented, default-off bit-identical

Motivation: w2lvar stage-5 walks asymmetrically (left duty 0.81 vs right
0.56) — per-side drive is the natural knob. Implemented this session
(see CHANGES_MANIFEST for the exact hunks):

- **GAMMA (RG-level)** enters through the four pre-existing TONIC current
  ports (`TONIC {L,R} RG {ext,flx}`) — value-level only, NO topology
  change; per side via new knobs. Defaults reproduce the symmetric
  baseline bit-exactly.
- **ALPHA (MN-level)** enters through NEW CONDITIONAL current ports
  `ALPHA {L,R} {Hip,Knee,Ank} MN {ext,flx}` (6 per side), built ONLY when
  `build(alpha_l=…, alpha_r=…)` != 0 (lab law: extra ports change
  compiled shapes even at zero current). Runner knobs feed them.

Knobs (all default = baseline): smoke_w2l.py `--tonic-e-l/-r`,
`--tonic-f-l/-r`, `--alpha-l/-r`; test_w2l_air_afferented.py and
test_w2l_ground.py `--te-l/-r`, `--tf-l/-r`, `--alpha-l/-r`.
**JSON RULE flag for the next study author**: any curriculum that tunes
these must save them in the winner json and add `if key in best` loader
branches.

Gates run (all on egg, post-edit):
- Bit-exact hashes pre/post-edit, 3 fresh-net repeats each
  [drive_gate_preedit.log == drive_gate_postedit.log]:
  bilateralrg 01929488aad3ec02… (95/208/10); split 57637001168ca2d1…
  (87/164/6); aff 6834c8d250363444… (97/186/12).
- smoke_w2l.py output byte-identical pre/post (PASS line above).
- M6 ground walk (`--phase=walk`) post-edit log IDENTICAL to the pre-edit
  baseline: bursts 11/10, period 1.539 s, r −0.776, tilt 8.2°, harness
  49%, ranges 33.9/29.3/63.3/66.8/38.3/36.9 [m6_gate_baseline.log,
  m6_walk_postedit.log].
- ALPHA net-level sanity [via _w2l_alpha_debug.log output in the session
  transcript; script kept in this folder]: aff alpha_r=0.5 → R MN
  window-mean 1.5933→1.9000 mV (L 1.5115→1.5101, untouched); bilateralrg
  alpha_l=0.6 → L MN 1.5123→1.8809 mV (R 1.5885→1.5874).
- GAMMA closed-loop: `--tf-l=4.75` walk → tilt 6.9°, knee-L 65.2°, L
  RG-E period 1.513 s (vs 8.2 / 63.3 / 1.539 baseline) — knob reaches the
  plant [m6_walk_gammaL_demo.log]. ALPHA closed-loop: `--alpha-l=0.3`
  walk → regime change (harness 0%, shuffle verdict) — reaches the plant;
  0.3 nA OVERDRIVES this configuration (no RG bursts in the analysis
  window); usable ranges are a tuning-study question
  [m6_walk_alphaL_demo.log].
- Free-running asymmetry probe: te_r 3.0→2.0 stops the pair's alternation
  (1 burst/12 s) — the two RGs are c1-locked, so single-side gamma moves
  BOTH periods; usable asymmetry window in closed loop is narrower than
  the open-loop one.

## 4. NOT RUN (honesty list)

- No 92-muscle gait2392 wiring from w2l_cpg (the TODO item) — unchanged,
  out of scope.
- No retune/optimization of the new knobs (demonstrations only).
- The rules' Ia/II central edges and per-muscle motifs (ia_homo, ib_auto
  inhibitory, ii_exc/ii_inh) — not implemented in w2l_mujoco; not built
  here either (Ben's rules-pass decision).
- w2lvar's own builder (build_network_w2lvar.py) was NOT touched — its
  L/R-drive parity is a follow-up.

## 5. Reproduce (egg)

```
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg
D:\Anaconda\envs\myo\python.exe smoke_w2l.py
cd /d ..\w2l_mujoco
D:\Anaconda\envs\myo\python.exe build_w2l_aff_net.py
D:\Anaconda\envs\myo\python.exe test_w2l_air_afferented.py --amp=1.5 --pulse=0.5 --ibnA=1.0 --pamp=4.0 --causal=both
D:\Anaconda\envs\myo\python.exe test_w2l_ground.py --phase=both
cd /d ..\reports_20261016\w2l_audit
D:\Anaconda\envs\myo\python.exe _w2l_rules_audit.py
D:\Anaconda\envs\myo\python.exe _w2l_drive_gate.py
```

— connectome auditor (robot), 2026-10-16
