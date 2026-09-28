# Goal 4 — easteregg2 generations campaign: bring-back report

**Date:** 2026-09-27 · **Host:** easteregg2 (verified via `hostname` this session) · **Branch HEAD:** `b084a6ab` (2026-09-26 00:25:50, unchanged by the campaign — zero commits made)
**Campaign:** continue the w2lvar/syn6 curricula on easteregg2 under Ben's stop policy, after Ben ruled the kine_ref left-reference fix ("option a": fix first, scores move). Everything below was produced by the goal-4 campaign + its confirmer/gate passes; this bring-back pass verified paths and spot-read the logs named inline. Items marked **[verified this session]** were read/run by this pass; everything else is quoted from the campaign's claims record, its logs on disk, or the supervisor gate, and is attributed as such.

**NO git commits were made.** All campaign output is uncommitted working-tree state for Ben's GitHub Desktop.

---

## 1. Executive summary

- **kine_ref left-reference fix APPLIED** (uncommitted in the working tree; `git status` = ` M Code/MuJoCo_SNS/spinal/kine_ref.py` **[verified this session]**). Every kine score quoted before 09-26 shifts once it is committed. Old→new anchors: s3k −160.23425729850192 → −161.56754173676563; w2lvar s5 winner −165.040506 → −165.335963; syn6 s5 seed −197.2146 → −200.3796; syn6 s4 winner −237.2534 → −234.2860.
- **Convergence headline: w2lvar stage-5 has a NEW BEST ground walk, −144.464403 (trial 68)** — a +20.9 like-for-like gain over the rebased old winner (−165.336), from 100 new trials under Ben's two-consecutive-sub-0.5-batch stop rule. **syn6 stage-5 refused to beat its seed** (−200.3796 rebased) after 40 more trials — the block is architectural, not a sampling deficit. syn6 stage-4 improved to −180.844766 (trial 109, +140 trials, +53.4 like-for-like).
- **Connectome schematics restructured into nested subsystems** (walker_s3k / w2lvar / syn6 templates), with three real generator wiring errors corrected and every claimed error re-verified against fresh mines of Ben's `.aproj`/`.asim`. Two **new** findings in Ben's AnimatLab files (report-only): the RE clone-type bug is **.aproj-only**, and `patch_comm_types.pl`'s .aproj link-repoint **never landed** — the GUI project and the runnable `.asim` encode different commissural strengths.
- **Numbers audit: every quoted Xi/GoF/torque number reproduced from source** by the confirmer pass (own MATLAB R2025a runs + scipy mat dumps), independently re-derived by the supervisor gate. One material discovery: `minimizeFlxPin.m` executes **pitch-only for BOTH flexor brackets** while its header and the AGENTS 09-13 note say origin two-rotation — Ben's fix-code-vs-fix-docs ruling is still open. Two small tex corrections are owed (+5.02 %, Appendix C line refs).
- **add_mag3_r ruled per Ben's torque rule: NEITHER-PASSES-BOTH.** Tie-break applied on Ben's behalf (overrulable): the rounded repo P2 is a corruption artifact (duplicates P2_0); the **full-precision thumb P2 was APPLIED** (one-line change to `gait2392_robotbody.osim`), `gait2392_robot.osim` regenerated and probe-verified bit-exact.
- **BiPulley: DONE, nothing left running.** Both sanity gates PASS; both production runs exited 0 (16:41→17:02 on 09-26). Flexor pulley result `Bifemsh_20mm_Result_pulley_20260926_1652.mat`: min torque margin **+5.38 %** (clears the 5 % requirement). BiPulley 5-muscle batch = first-pass baselines (biarticular routes do not reach target at the straight G=1 config; clearance violations reported by design).
- **AnimatLab: 12 "missing" synapses restored** byte-verbatim into `Biped_2xCPG_wSubs.aproj` (they had been purged mechanically by Ben's Sep-16 GUI re-save, not deleted by anyone). **Li closed loop: still collapses** — root cause measured: the MJCF Root has **no joint** (welded base, heels can never touch → contact-driven CPG never triggers).
- **Supervisor gate: GO** (6/6 checks; the single harvest "FAIL" line is a stale-reference false positive, explained in §3).

---

## 2. kine_ref fix and old→new baselines

**The fix (working tree, uncommitted).** `kine_ref.load_reference()` used to cut each leg's cycle at the first GRF onset; for the left leg that onset (t = 0.005 s) precedes the IK file's support (starts 0.50 s), so `np.interp` edge-fill froze `ref['l']` for the first ~39.6 % of the left cycle (proof `reports_20260925\figs\ref_left_bug_proof.png`, 09-25). Ben ruled option a (fix first). Fix now in `kine_ref.py`: `_first_pair_in_ik()` picks the first onset pair **fully inside IK support**, with a coverage assert. Fix-liveness probes (campaign): `tmp/probe_ref_now.py` → T_l = 1.2133000000000003 s; `tmp/diff_ref_fix.py` → right cycle 0.00° max diff (left-only change), previously frozen leading left segment now 16.52° spread.

**Old→new baseline table** (like-for-like anchors under the fixed reference):

| Eval | Old (buggy left ref) | New (fixed left ref) | Δ | Provenance |
|---|---|---|---|---|
| s3k stock-chain ground eval (curriculum_stage3.json = curr_s3k_nocross t34 params, full_rules 0, drive 2.0773198435110425) | −160.23425729850192 | **−161.56754173676563** | −1.333284 | Campaign `tmp/rescore_s3k_fixed_ref.py`, run twice bit-exact; supervisor independently reproduced **both** values bit-exact (fixed tree and HEAD kine_ref via shim) |
| w2lvar s5 winner (trial 4) | −165.040506 | **−165.335963** | −0.295457 | **[verified this session]** `easteregg2/ee2_rebase_w2lvar_s5.log` lines 10–12; kz 0.8297 / tilt 27.72 (no penalties), duty 0.460, 3 bilateral cycles/leg, knee_min −18.14°, contact_frac 0.566/0.815, RG_E_r rises(t≥5)=2 |
| syn6 s5 winner (trial 0, the seed) | −197.2146 | **−200.3796** | −3.1650 | **[verified this session]** `easteregg2/ee2_rebase_syn6_s5.log` line 9; guard rails all pass (contact 0.955/1.0, kz 0.8396, tilt 28.07, rises=3) — log verdict "GENUINE GROUND WALK (rebased)" |
| syn6 s4 winner (trial 14) | −237.2534 | −234.2860 | +2.9674 | Campaign `tmp/rebase_syn6_s4_winner.py`; log `easteregg2/ee2_rebase_syn6_s4.log` (on disk **[verified this session]**); supervisor re-derived in CHECK 6 |

Notes:
- The s3k reproducer of record is `_pf_layer_variant_test.py` gate 3's flow — but the driver had to use `C.set_stage(4, ...)` because the 09-24-night reorder (standing=3, walking=4) moved the ground-walk key set to the `stage in (4, 5)` block; verified line-by-line against pre-reorder commit `2b1d679a` that old set_stage(3) applies the same keys. `git diff 2b1d679a..HEAD` over runner/params/build_network/optuna_walk_v10/muscle_map/kine_ref/best_walk_params_v10/curriculum_stage3 is empty, so kine_ref was the only delta vs the documented baseline. `spinal_run.npz` (mtime 2026-09-25 22:53:39, 8,238,797 B) and `optuna_walk.db` untouched; the rescore npz went to `scratch_s3k_rescore.npz`.
- **Not done by the campaign or this pass:** re-verifying the old −160.23425729850192 by reverse-patching the fix in place (it is uncommitted on disk; the old value stands as documented and was separately reproduced bit-exact from HEAD kine_ref by the supervisor's shim).
- **Stale gate:** `_pf_layer_variant_test.py` gate 3 is stale **twice over** (renumber: `set_stage(3)` now loads the balance stage → prints kine −25.0 FAIL; plus the fixed ref). Supervisor re-ran the full story: gates 1–2 PASS ((410,376,1186); variant delta (2,0,52) TOEDF r/l), gate 3 FAIL as-is, set_stage(4)+walk params = −161.56754173676563 on the fixed tree and −160.23425729850192 bit-exact on HEAD kine_ref. Repoint gate 3 at stage 4 and re-record post-commit.

---

## 3. Convergence per variant/stage (trials added, best before→after, stop reason, harvest verdicts)

Stop policy (Ben): **stop a stage after two consecutive batches improving the study best by < +0.5** (max 10 batches); STEP-3 air-stage rule: one +20 top-up if the in-flight improvement since the 09-25 reference > +1.0. Mixed-units caveat applies to s4/s5 of both variants: trials 0–22 (s4) / 0–21 (s5) were scored under the OLD left reference, later trials under the FIXED one, inside the same studies; stop-rule bookkeeping used the rebased winner replays (§2) as same-units baselines, and DB `best_value` is reported as-is.

### w2lvar (`optuna_w2lvar.db`; chain npz `spinal_run_w2lvar.npz`)

| Stage | Best before → after | Trials added (total) | Stopped because | Harvest verdict |
|---|---|---|---|---|
| s1 air-deaff | 127.752845 → 127.752845 (t36) | +20 (71) | STEP-3 single top-up (in-flight +46.9 > 1.0); batch changed nothing | REAL RHYTHM, rises=30; **objective saturates the rises=30 flutter cap** (score = 3·30 + 0.5·75.5 = 127.75 is the ceiling) — further s1 batches can only chase knee-flexion depth |
| s2 air-aff | 80.589509 (t22) → 94.924507 (t49) | +20 (68) | STEP-3 single top-up (in-flight +20.4 > 1.0); batch added +14.33 | REAL RHYTHM, rises=20; no winner sits at the ≥3 gate minimum |
| s3 balance | 25.613422 (t35) → unchanged | 0 (48) | in-flight improvement vs 09-25 reference = +0.536 < 1.0 → left as-is | genuine stander (exploit gate PASS) |
| s4 walk-no-contact | −225.890865 (t17) → unchanged | +40 (63) | two consecutive batches +0.0 (all 40 new trials ≤ −290.11); stopped well before the 10-batch cap | winner is a genuine eval (fixed-ref −230.860, tilt 1.07, duty 0.64; the study's 3 real-kine trials: −225.89 / −271.96 / −290.11). The harvest log's "exploit gate: FAIL" line is a **stale-reference false positive** — repro under the fixed ref gives −230.860209 vs the json −225.890865 recorded under the old ref (Δ 4.969 = exactly the kine_ref shift); **[verified this session]** `easteregg2/ee2_harvest_w2lvar_s4.log` lines 10–11 |
| s5 walk-contact | −165.040506 (t4) = −165.335963 rebased → **−144.464403 (t68)** | +100 (122) | two consecutive <0.5 batches (b4, b5); 5 of 10 batches used | batch path: b1 −163.906253 (t40, +1.13) → b2 −145.744493 (t60, +18.16) → b3 −144.464403 (t68, +1.28) → b4 0 → b5 0; **NEW BEST ground walk of the variant program** (bilateral, kz 0.830, 3 cycles/leg) |

`curriculum_w2lvar_stage5.json` now records trial 68 at −144.46440311078322 **[verified this session]**. Supervisor CHECK 3 re-played the rebaseline target (trial 4) from the pre-campaign snapshot `pre_w2lvar_stage5.json`: kine_score −165.33596254795964 vs the tuner log — **delta +0.000e+00, bit-exact**.

**Fragility caveat (carry to Ben):** the s5 winner sits near search-space edges — drive 3.82 (ub 4.0; trial 60 reached 3.99), ib_to_mn_inh 0.014, contact_onset 0.003 (near 0) — and trial 68 landed only two batches after a +18 jump. Treat as fragile until independently replayed/verified by Ben.

For orientation, the 09-25 report-time w2lvar winners (trials 17/13/0/17/4 = 80.835586695 / 60.230492583 / 25.077366853 / −225.890864513 / −165.040505521) were confirmed in the DB to <1e-9 by the program-numbers confirmer (sqlite mode=ro); s1–s3 were superseded by the concurrent in-flight extension **before** the campaign's top-up batches, which is why the table's "best before" for s1–s3 already reflects the in-flight values.

### syn6 (`optuna_syn6.db`; chain npz `spinal_run_syn6.npz`)

| Stage | Best before → after | Trials added (total) | Stopped because | Harvest verdict |
|---|---|---|---|---|
| s1 air-deaff | 95.639 (t17) → 130.238130 (t57) | +50 | STEP-3 single top-up (in-flight +33.668 > 1.0; top-up added +0.931) | winner replay bit-exact (+0.0000), RISES=30, REAL RHYTHM |
| s2 air-aff | 86.751 (t12) → 110.517696 (t59) | +50 | STEP-3 single top-up (in-flight +20.873; top-up added +2.894) | winner replay bit-exact (+0.0000), RISES=28, REAL RHYTHM |
| s3 balance | 37.463 (t11) → unchanged | +30 (48) | in-flight +30 left the best EXACTLY at the 09-25 winner (+0.000 < 1.0) | floor census: −150(fall)=5, stood=43, 0 NaN |
| s4 walk-no-contact | −237.2534 (t14) = −234.2860 rebased → **−180.844766 (t109)** | +140 (163) | STOP after 7 of 10 batches: b6 and b7 both +0.000 | batch path (post-fix trials under FIXED ref): b1 −201.682 (t32, +32.60 vs rebased baseline) → b2 −199.638 (+2.04) → b3 −182.432 (t79, +17.21) → b4 0 → b5 −180.845 (t109, +1.59) → b6 0 → b7 0; net like-for-like **+53.44**; harvest: 56 frozen sentinels of 163, 0 NaN |
| s5 walk-contact | −197.2146 (t0) = −200.3796 rebased → **unchanged — SEED STILL WINS** | +40 (62) | STOP after 2 batches, both improving the best by <0.5 | batch bests: t30 −216.262, t51 −242.010 — both BELOW the rebased seed; study best remains trial 0; harvest: 36 frozen/penalized of 62, 0 NaN. Confirms the documented syn6 weakness (winners at the rises=3 gate minimum) — **architectural, not a sampling deficit**; no batches were forced |

`curriculum_syn6_stage4.json` = −180.8447656274239 (t109); `curriculum_syn6_stage5.json` still records trial 0 at −197.21456977796157 (old units — it remains the study argmax) **[both verified this session]**. Supervisor CHECK 2: all 10 studies' DB argmax over COMPLETE trials equals the recorded curriculum winner (exact trial numbers, 0.000e+00 value deltas, params match to 1e-12).

Chain discipline: batches strictly sequential per stage, never two chain stages concurrent; AARL_NPZ pinned to the chain npz for chain runs and to `scratch_*.npz` for replays. Protected files untouched: `optuna_walk.db`, `spinal_run.npz`, `reports_20260923/`, `reports_20260924/` (supervisor CHECK 5; git status carries no entries for them; mtimes/sizes identical).

---

## 4. Connectome schematic changes + wiring errors found and disposition

### Schematic change (landed, verified)

The three walker templates in `connectome_templates.json` were restructured into **nested subsystems** (editor v3.3 `node.sub`): top-level censuses are now **walker_s3k 13 nodes / 322 edges** (10 SUB wrappers + 3 flat: DRIVE/POSTURE/PM), **w2lvar 18 / 2726**, **syn6 18 / 2123** (10 layer SUBs + 8 flat each); 30 wrapper SUBs per template (10 semantic layers + 20 per-pool, e.g. "hip_ext pool (L)"); three new `_flat` entries preserve the pre-schematic wiring verbatim (== HEAD). Conservation independently checked: flattened leaves 233/980/886 == HEAD with (label,type,grp) multisets equal and x/y preserved; edges 420/7526/3129 == HEAD with (sign,gain,tag) signature multisets equal. Subsystem load test `tmp/_tpl_subsys_load_test.js`: **42/42 PASS** (3,910 node loads / 19,344 edge loads / 92 subsystem enters). Shipped guard `tmp/_tpl_diff_schematic_20260926.py`: **21/21 PASS**. Diff vs git HEAD: 12 pre-existing templates canonically byte-identical; `li`/`w2laproj` differ only in `_note`/`_mine_provenance`; `bilateralrg` intentionally changed 105n/226e → 109n/222e (V3 correction below). Editor suites re-run this campaign: `_editor_edges_test.js` 17/17, `_editor_template_test.js` 4/4, `_editor_tree_test.js` 9/9, `_editor_static_check.py` (duplicate ids NONE; missing `{crumbBack}` pre-existing dynamic id), `node --check` on the last script block PASS (75,527 chars by two extraction methods — the claims' 75,554 was a 27-char bookkeeping delta on their side).

### Wiring errors found → disposition (all verified against fresh mines of Ben's `.aproj`/`.asim` and the `tools/build_*.pl` chain)

1. **V3→contra-RG-E HC edge was a generator invention.** Real model (build_comm.pl:12-13,46,48): V3 → contralateral **RG ext IN (InE)** only, 2 out-edges per side. Template corrected to `V3_L→'R RG ext IN'` / `V3_R→'L RG ext IN'` exc 0.1 (make_editor_templates.py:309-312). FIXED in templates.
2. **heel/toe→ANK MN ext was invented.** build_contact.pl defines exactly 20 links to {RG ext, Hip/Knee PF ext, Hip/Knee MN ext} — zero ankle. Template filter corrected; `contact_C20` = 20 edges with exact targets. FIXED.
3. **Afferent→HC set now matches build_aff.pl @AFF 1:1**: 24 edges, sources = Ia relay + NEW SN-II relays + Ib, own-joint PF + RG half-center, all exc 0.01. FIXED.
4. **R RG↔RG direct exc pair** now in template step 1 (make_editor_templates.py:272-273); real .aproj type census has "RG to RG Excite" = 4 links; the `build_w2l_net.py` mirror block adds 0 (verified conditional). FIXED.
5. **RE clone-type bug** ("R Hip MN flx RE → R Hip MN ext RE" typed "V3 Commissural Excite" equil −0.04, vs L mirror "RE to RE Inhibit" −0.07): CONFIRMED in the .aproj — but **.APROJ-ONLY**; the `.asim` connexion is correctly "RE to RE Inhibit" (Equil −70, SynAmp 0.5). Earlier ".aproj/.asim" phrasing overreached. Practical impact narrower than first claimed: the runnable standalone is clean. REPORT-ONLY (Ben's file); his one-type GUI fix applies to the .aproj alone.
6. **Contact-layer semantics vs Ben's 2026-09-24 ruling** (heel = stance reset at PF layer via InE/InF + PF-layer INs; toe = DF-inhibition only): the 09-16 W2L build instead has heel AND toe uniformly exciting the 5 extensor layers directly, no DF-inhibition chain, no ankle PF layer. Faithful transcription of the 09-16 build — REPORT-ONLY; the walker-side implementation of the ruling lives in `build_network.py` default-0 keys, outside this template.
7. **`w2l_mujoco/build_w2l_split_net.py` still hardcodes the invented edges** (lines 119, 122: V3_L→'R RG ext', V3_R→'L RG ext'; stale docstring lines 20-22). Its 09-25 march results were produced with them; removing them changes a live tuned net. REPORT-ONLY, pending Ben's go.
8. **%TEMP%-loss fix verified real**: make_editor_templates.py:1219-1241/1263-1285 load `w2l_aproj_mine.json` / `li_aproj_mine.json`, both verified faithful vs fresh parses (w2l 92n/163e = 150 functional links 1:1; li 58n/98e). **CAVEAT: both mine files are UNTRACKED in git** (`??` status) — "tracked" becomes true only when Ben commits; until then the regeneration contract rests on uncommitted files.
9. **Gain-convention fuzz** (SESSION_NOTES_20260916.md:57-61: .asim effective strength = type SynAmp, per-connexion G inert; README documents live-net gains as calibration): REPORT-ONLY.

**NEW discrepancy found by the confirmer (report-only, Ben's file):** `patch_comm_types.pl`'s **.aproj link-repoint never landed** in the committed `.aproj`. All four commissural out-links (c1→contra RG flx ×2, V3→contra RG ext IN ×2) still ride the strong "RG Excite"/"RG Inhibit" types (maxcond 2.749e-6) even though the dedicated weak types exist in the same file with the intended numbers ("c1 Commissural Inhibit" 2.749e-7 — 0 links; "V3 Commissural Excite" 1e-7 — 1 link, the mis-typed RE link). The `.asim` IS fully patched. SESSION_NOTES' own latch warning ("with RG-Excite (SynAmp 2.749) the V3 path latches both RGs") applies to the GUI project. **One GUI session can fix all 5 link types** (the RE clone + the 4 commissural out-links).

**No additional discrepancies:** the name-normalized template-vs-.aproj diff is 208/208 functional links with exactly ONE residual (the RE clone sign flip). Ia reciprocal, Ia↔Ia, Ib autogenic EXCITATION, II excitation, Renshaw (MN→RE/RE→MN/RE↔RE), PF drive/cross and RG lamination all match Ben's real model exactly.

**Live net re-verified:** `w2l_cpg/build_w2l_net.py` → 95 neurons / 208 synapses / 10 inputs / 12 outputs, per-tag counts as claimed; `smoke_w2l.py` → "W2L_SMOKE PASS period=1.000 antiphase_r=−0.525 rg_e_bursts_L=22 rg_e_bursts_R=22". Note for quoting: **antiphase_r is −0.525 on easteregg2** (deterministic across thread counts/cwds, numpy 1.21.6); the 09-24 report's −0.532 reproduces only on the machine that produced it — gate verdict unaffected (threshold r < −0.5).

---

## 5. Numbers-audit corrections, with provenance

Two confirmer passes ran (Xi/GoF/torque; program numbers), and the supervisor gate re-derived the decisive ones. Provenance per item; nothing was unreproducible except where stated.

### Xi / GoF / torque (evaluators, mats, logs)

- **Flexor pick-77 GoF reproduced exactly** (own MATLAB R2025a U1 runs of read-only reeval scripts; logs `Dig_out/xi_audit_tmp/confirmer_reeval_gof.log`, `supgate_reeval_20260927.log`): per-test RMSE 2.16465/1.36817/2.28691/1.61102/1.08363; FVU 0.171825/0.0468165/0.164141/0.0426721/0.0277468; MaxRes 5.4668/3.61397/7.8177/5.89053/1.7662; g77 = 3.93736985e-3 m, 39978.7285 N/m, 14733.7523 N/m (scipy dump of `minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat` row 77, xCols [4 5 6]).
- **Old row-107 appendix numbers (1.94/1.51/2.29/1.63/1.21; FVU 0.14) are NOT reproducible** under the current evaluator. Fresh row-107 truth: RMSE 2.01171/1.49305/2.314/1.73522/1.08103 (FVU 0.1484/0.0558/0.1681/0.0495/0.0276). The dissertation's replaced GoF paragraph is correct as edited.
- **Extensor flxr77 pick-32** (pool {1,2,5,6,7,8}): RMSE 1.24371/1.34986/0.977636/1.52455/0.48955/0.699429; FVU 0.130/0.888/0.163/1.327/0.053/0.171 — FVU>1 only ke6 inside the pool (excluded test 4 = 1.22092 also >1); baselines 3.09–3.62 (baseline FVU row 6 = 4.7304). Reproduced to every digit. g32 = −6.378022e-3 m, 39978.7285, 14733.7523, 0.158289 (Xi0 inside the Ben-2026-09-20 bounds [−2,+0.5] cm).
- **Bio-extensor (l0 = 51.8 cm) at the adopted values: 1.1702/0.6681/2.2088** (FVU 0.67 < 1) vs old adopted 2.1950/2.3506/4.7268 (= the old dissertation numbers, pipeline check) and baseline 3.0322/4.4856/5.5870; attribution mix (old Xi0/Xi3 + new pair) 2.1930/2.3462/4.7268 → the gain is the new Xi0/Xi3. The old "open gap" sentence was factually superseded; tex replacement correct. Other drafts still quoting 2.20/2.35/4.73 as current remain flagged for Ben.
- **Margins:** 09-20 flexor adoption +0.050000 (+5.00 %) and +5.043 % at −20.202° — exact; Bifemsh 09-21 mat min +5.007 % at −106.869° — exact; **Vas 09-25 mat minMargin 0.05024928 = +5.0249 % → the tex "+5.03 %" (30-results.tex:355 + Appendix C narrative) should read "+5.02 %" (or +5.025 %)**. CORRECTION OWED (verification pass made no tex edits).
- **MATERIAL: both flexor brackets execute pitch-only.** `minimizeFlxPin.m`: line 306 `Thbr = RpToTrans(RhbrZ, Pbr2')` LIVE (two-rotation commented 305); pitch-only `pbrBnew` 326 / `pbrAnew` 339 live (variants commented 325/327/338); `K2 = [X1, X1, X2]` at line 400 — **not y-symmetric, so the frame choice is numerically material** — while header line 7 and the AGENTS 09-13 "mixed convention" note say origin two-rotation. Created in this state by `8e3c1b9d` (2026-09-13 22:26:44 -0700; working tree == HEAD). Everything since 09-13 inherits it. **Ben's ruling open: fix code and invalidate, or fix docs only.**
- **Tex line-ref corrections owed** (substance unaffected): 94-AppendixC.tex bracket-frame paragraph cites 305/304/318/331/330; actual lines are **306 / 305 / 326 / 339 / 338, plus commented 327**.
- Confirmed as-live: extensor `K = [X1, X2, X2]` (Ben 2026-09-16) at minimizeExtX3.m:575; Pbri = [−48.11, −107.81, 13.8]/1000 at :261 (old [−27.5, …, −0.54] commented 262); bio rows minimizeExt.m:352 [X1,X2,X1] / minimizeFlx.m:314 [X1,X2,X2]; settled-Xi row of record = flexor pick **77** (buildKneeFlexorContext20mm.m:178, 107 commented) + ext pick **32** (:165), LOCKSRC/flxRow/xi1lock/xi2lock verified in-mat; both builders `requiredTorqueMargin = 0.05` (:144/:152); appendix static tables exact (kf/ke constants; M = 2 flexor / 7 extensor path points; Ak [−120,+10] / [−125,+35]); extensor seed table (fcec p3 64.06/−426.80, p4 52.83/−441.88, p5 37.77/−452.14) and release schedule (eliminatedAngleD p3 +6.06, p4 −55.66, p5 −89.80, p6 −47.78, p7 −20.20, p8 +7.37; "−p8") verified vs `Vas_Pam_20mm_Result_20260925.mat`.
- **No LaTeX compile exists on easteregg2** (PATH, D:/MiKTeX, D:/Program Files/MiKTeX, user-profile MiKTeX all checked) — structural smoke checks pass; **compile on EB475WS4 before Overleaf upload remains required.**

### Program numbers (spinal / AnimatLab / SNS_Simscape / libraries)

- Verified from committed mats/logs: KneeReflexDemo settle **43.4516°** (41.6104–44.8709, t>3 s); BPACPGLegDemo **9.7066–48.0143°**; BeerCup ON **+1.9965°** / OFF **+7.8464°**, A 0.3924→0.4170, At 0.4000→0.3163; Deng CPG period **1.938 s** (10 bursts / 20 s).
- **Units test 6.07e-4 mV double-verified**: the confirmer rebuilt SNS_Library + actuators from the committed build scripts in a temp copy and ran `sns_units_test_2n` on R2025a → PASS (max|dV_A| 6.07e-04 mV).
- **SNS_SpinalNetwork "verify 4.2e-6 mV" still has NO in-repo artifact** (README_SNS_Simscape.md only; no log line anywhere; spinal_net_export.json not committed; the committed slx is R2025b-format, unloadable on this R2025a machine) — NOT re-runnable here. Committing a verify log would close this.
- BilateralRG ground run (L 10 / R 11 onsets; period 0.4640/0.4641 s; L→R lag +0.474 cycle) and W2L modern-asim numbers (PF max −55.83 mV, RG −56.41 mV, 2.222 Hz everywhere; lags +57.0 ms r=0.916 / −63.0 ms r=0.924 / −89.0 ms r=0.900) verified on BOTH the report backup and the live chart byproducts. 2023 asim: 8 bursts each, median 1.302/1.299 s; the "40 N push" is ForceInput 40 N over t∈[0,1] s inside a 10 s run (AGENTS' "survives 40 N push 10 s" is loose wording).
- Li reference: period **1.305 s** (L_foot, 8 onsets) / 1.290 s (R) — verified from `Li Model\DataTool_7.txt`.
- **Goal-5 library scores** (goal5_allrefs_results.csv, 44 rows): subject01_of_record −189.45 (best), Case_40_motion −200.46, ong_speed_050→200 −203.32 … −224.42. **Correction: "monotone in speed" → "near-monotone"** — ong_speed_100 (−207.37) scores better than ong_speed_075 (−208.15). AGENTS carries the same wording.
- w2lvar/syn6 DB trial values quoted in goal4 confirmed to <1e-9 (sqlite mode=ro). The studies were also extended concurrently past midnight during the audit (completions through 2026-09-27 00:21–00:22) — final totals above are the closed-chain values.
- `muscle_force_compare.csv`: 8/12 in band verified (0.8836–0.9770), out-of-band exactly as documented (soleus 0.742, semimem 0.565, rect_fem −1.150, bifemsh −0.261); in-band mj/os ratios 0.9096–1.0180, lags 0–2 frames. **Header bug is 6 names vs 8 fields — THREE unnamed trailing columns** (two force columns + one numeric); cosmetic fix before anyone re-parses.
- Dissertation ProofFinal muscle-force corrections verified present and CSV-consistent (CPG_spinal_section_draft.tex:408-414 "eight of twelve … lags of 0–2 frames"; muscle_model_appendix_draft.tex:99-104 "±10 %" + semimembranosus clause 119-120) — **Overleaf sync still pending** (ProofFinal is local-only; upload mirror untouched per house rules).
- **NOT re-verified (carried as unverified):** the dissertation's s3c trial-54 gait numbers (nine cycles, 1.29 s, duty 0.87, best −149) — outside both confirmer passes' run sets; `_run_t54.py` exists if Ben wants it checked.

---

## 6. add_mag3_r verdict per Ben's torque rule

Ben's rule: the muscle must meet-or-exceed human torque on BOTH hip axes. Method: OpenSim 4.6 (`D:/Anaconda/envs/opensim/python.exe`), two variant .osim copies built by XML edit (originals untouched), grid = the project's hip RoM for this muscle (Add_Mag_Mesh_Opt.m:41-52: flexion −25..+85°, adduction −45..+20°, 322 poses), moment arms validated vs central-difference dl/dq (8/8 exact), primary torque convention tauA = |arm| × 488 N (the project's own human target: Adductor Magnus 3, MIF 488 N, Add_Mag_Mesh_Opt.m:93-99).

**Verdict: NEITHER-PASSES-BOTH — the candidates trade torque between the axes** (tauA peaks):

| Variant | Hip flexion | Hip adduction |
|---|---|---|
| repo P2 (−0.059, −0.108, −0.03) | 19.14 N·m = **0.63× human 30.61 — FAIL** (arm crosses ZERO at ~45° flexion, sign flips +3.5→−2.8 cm across the RoM) | 49.96 N·m = 1.49× human 33.60 (above human at 319/322 poses) |
| thumb P2 (−0.1357, −0.0929, +0.0591) | 34.28 N·m = **1.12× — PASS** (above human at 252/322 poses, min ratio 0.91 only at the 85° corner; sign-consistent like human) | 27.62 N·m = **0.82× — FAIL** (below human at 310/322 poses) |

(Equilibrium-convention tauB splits the same directionally; its peaks are ~10× inflated for BOTH variants by passive overstretch — lce 0.253 m = 1.93× OFL at default pose → 4509 N vs 488 N MIF — so they do not discriminate; flagged as a separate model-quality issue.)

**Tie-break (escalated per instructions; applied on Ben's behalf, overrulable):** (1) the rounded repo point is disqualified as a **corruption artifact** — it numerically duplicates the femur-side P2_0, giving a zero-length via segment → the flexion null + sign flip; its 1.49× adduction win is broken geometry nobody chose; (2) the **full-precision thumb point is KEPT and APPLIED** — one-line change to `gait2392_robotbody.osim` `add_mag3_r-P2 <location>` (git diff: 1 insertion, 1 deletion; file was git-clean before); (3) the 18 % adduction shortfall (0.82×) is routed to the **Opt_run/Mesh_Optimization MIF/route-sizing lever** — a sizing gap, not a path choice.

Applied + regenerated + verified: `repair_robotbody_muscles.py` re-run ("mirrored 46 right->left GeometryPaths", "non-mirrored pairs: 0", "VALIDATION OK"); `add_mag3_l-P2` = exact z-mirror; OpenSim probe of the regenerated `gait2392_robot.osim` == thumb variant **bit-exactly (max dev 0.00e+00)** at 4 grid-spanning poses. Verdict note: `Solid_Models\OpenSim\Gait2392_Robotbody\add_mag3_r_torque_verdict_20260926.md` (on disk **[verified this session]**); full artifacts in `add_mag3_compare_20260926\` (variant .osim copies, arm_tau CSVs, logs — all present).

**Open / downstream (Ben):** overrule window (if he weighs adduction above flexion, revert the one line — variants preserved); sizing item for Opt_run; the ~10× passive-overstretch finding (robot P1/TSL/route) needs its own decision; **stale one path-point**: `gait2392_robotbody_hip.osim` (side copy), `gait2327.osim` (generated), the `mjc\` MuJoCo conversion and any spinal references built from it — regeneration Ben-gated (heavy MyoConverter route; protected spinal npz/db untouched).

---

## 7. BiPulley status (incl. anything left running)

**Sanity gates — both PASS** (MATLAB R2025a `-batch`, exit 0): `Opt_sanity_pulley` — six assertions (G=1 regression identity ~1e-16; G=2 closure 1.137e-13 N; nBPA=2 tension-sum 4.5e-13 N; reaction identity; bowden u_t constant; 52/92 frames infeasible flagged+NaN), Xi loaded live = the flxr77 locked pair (0.00393737 / 39978.7 / 14733.8). `Opt_sanity_BiPulley` — geometry vs independent recomputation exact; G=1 mono reduction ~1e-15; G=2 distal transmission; `biPulleySpecsFromOpenSim` 5 muscles with 2-BPA bundle centroid on-line; `Opt_run_BiPulley` dry-run resolves 5 muscles. Logs: `Testing_Data\2022_02_Festo\Dig_out\pulley_campaign_20260926\{sanity_pulley,sanity_BiPulley}_20260926.log`.

**Production — both runs COMPLETE, exit 0, nothing left running** (launch record `launch_record_20260926.txt` **[verified this session]**; sequential, one pool at a time):

- **LAUNCH1 `Opt_run_pulley`** (flexor; cwd `Testing_Data\2022_02_Festo`; in-file production settings kept verbatim, smoke/configs env explicitly unset): 16:41:44 → exit 0 at 16:52:24 (10.7 min). Saved `Code\Matlab\Mesh_Optimization\Results\Bifemsh_20mm_Result_pulley_20260926_1652.mat` (mtime verified) — fBest **0.0846509**, xBest rest 0.56251 m / tendon 0.025 m, config nPulleyBPA=2, G=1, moving_via; ctx Xi = flxr77 locked pair (+Xi3 0.158289); exitflags surrogateopt 0 (eval-limit, as budgeted) / patternsearch 1 (converged). **Headline: minimum torque margin vs human target (insertion torque) = +5.38 %, mean remaining shortfall 0.00 % — clears the required 5 % margin.**
- **LAUNCH2 `Opt_run_BiPulley`** (cwd `Code\Matlab\Mesh_Optimization`; the only code change of the campaign: lines 77/80 RUN_BATCH/FULL_RUN false→true — the file's own documented production values; parpool(6) pre-opened): 16:52:24 → exit 0 at 17:02:40 (10.3 min). All 5 OpenSim muscles completed, checkpoint 5/5. `BiPulley_batch_summary.csv` **[verified this session]**: med_gas_r obj 101 / worstC 10 / margin −1 / ps −2; soleus_r 0.996794 / 0 / −0.997 / 1; tib_ant_r 0.928794 / 0 / −0.929 / 1; rect_fem_r 3.42299 / 0.242 / −1 / −2; bifemsh_r 101 / 10 / −1 / −2. **Honest reading: margins near −1 mean the straight G=1 config at the batch's hand-tune Xi neighborhood does NOT reach the fmax×5 cm nominal-arm target on these biarticular routes** — first-pass baselines for Ben's discrete config comparison, not tuned winners. Route-clearance violations reported by design (need ≥ 25.0 mm): soleus vs med_gas 19.1, tib_ant vs soleus 15.6, rect_fem vs med_gas 7.0, bifemsh vs rect_fem 4.3 mm.
- **Deviations / not done:** the pulley run's pool opened with **10 workers, not 6** (an implicit-parallel construct during the ctx build auto-created the default pool before the in-script `parpool(6)` guard; `parpool(6)` itself probed working; the BiPulley run opened exactly 6) — not killed mid-flight per house rule; master's-source manifest unavailable on easteregg2 (`D:\Bipedal humanoid` absent — designed OpenSim fallback used); the `OPT_PULLEY_CONFIGS` 6-config sweep left OFF per in-file defaults.
- **Anything left running: NO.** Both productions + both sanity gates exited 0 (chain DONE 2026-09-26T17:02:40); no MATLAB processes remain from the campaign; all byproducts (2 sanity logs, 2 production logs, launch record, chain script) are in `Dig_out\pulley_campaign_20260926\` **[all 6 files verified this session]**. The tuner chains also closed at their stop rules (final study totals in §3). Nothing is running on easteregg2 from this campaign.

---

## 8. AnimatLab + Li outcomes

### AnimatLab (Biped_2xCPG_wSubs)

- **Discovery:** the 12 "missing" synapse blocks were no longer IN the working file to relocate. Evidence chain: working .aproj == HEAD `8eb0efe8` (Ben's Sep-16 09:57 GUI save, 227 synapses, git-clean); the fix-commit pause state `53f3d3ec` had 261; key-set diff shows Ben's save removed exactly 34 keys and added 0; a CDATA-opaque scan shows the fix commit still had 40 loose synapses, 34 of which = exactly the 34 purged. Per the GUI loader semantics documented in CONTINUE_HERE, blocks outside a fragment's `<Links>` list never load → **Ben's re-serializing GUI save mechanically purged them; nobody deleted them deliberately.**
- **Fix applied** to `Biped_2xCPG_wSubs.aproj` (18 exact-range edits, CRLF-preserving, temp-copy-verified first): the 12 blocks restored byte-verbatim from `53f3d3ec` with original GUIDs (8 cross-joint Ia → Knee MN F into LH_Knee/RH_Knee `<Links>`; 4 commissural into a new Neural Subsystem `<Links>`), 4 commissural OffPage nodes re-added, every link ID registered in origin OutLinks/dest InLinks, page arrows rebuilt (incl. re-attaching RH_Knee's 4 surviving dead arrows), AddFlow headers fixed. Verification on the real file: whole-file XML PASS; 15/15 page CDATAs PASS; `verify_final.pl` — all pages declared==actual, 0 mismatches, 0 dangling/corrupted, **syn=239**; synapse set = pre-fix + exactly the 12 expected keys, 0 removed; git diff on the aproj: 904 insertions / 10 deletions.
- **GUI check PASS** (zero Error dialogs + status-bar "Load project complete"; GUI killed before any save). **Sim check, honest scope:** `AnimatSimulator.exe` runs flattened `.asim` only (the .aproj load validation is the GUI check); the model's only runnable export `Biped_2xCPG_wSubs_Standalone.asim` runs exit 0 with fresh chart byproducts — but that asim is a **STALE export (243 connexions ≠ the edited 239)**, so this validates the simulator/chart pipeline, NOT the 12 restored synapses. **Ben must re-export the Standalone asim from the GUI before any RG-oscillation re-test** (the earlier "latch" result was measured on a half-wired connexion set and is invalid per the file).
- **NOT done, reported:** 22 further synapses lost in the same GUI purge remain unrestored (6 K&A PF→BFlh/Semimem/RF MN, 1 RH_AnkleZ-Dorsi II→MN, 15 RH-container afferent/RG-E→PF broadcast) — same evidence trail exists, outside the ask's scope; RH afferent routing is the area Ben has been reworking by hand.

### Li closed loop (w2l_mujoco) — still collapses; root cause measured

Retried on the axis-fixed body (`w2l_mjcf_fixed.xml`; net census reproduced 56 neurons / 84 synapses / 20 inputs / 12 outputs exactly; Li gains/tonics/encoders verbatim, nothing touched): 3 s smoke → 0 stance episodes; full 20 s gate → **VERDICT "not yet"**: 0 stance episodes, duty 0.00, heel-contact ON 0.0 % L / 0.2 % R (log `reports_20260925\logs\test_li_stepping_fixed_run1_20s.log` **[on disk this session]**).

- **PRIMARY cause (single suspected):** `w2l_mjcf_fixed.xml` ships the Root with **NO joint** (njnt=8, all leg hinges; body_jntnum=0) — a jointless body is **welded to the world**, so the "ground run" is a fixed-base rig with feet hovering 2.6–3.1 cm above the floor (probe: pelvis at exactly z = 0.993 m for all 20 s). Li's CPG is contact-driven — the heel SN is the stance trigger — so no stance episode can exist. This contradicts make_w2l_mjcf.py:49-50's own note (the freejoint was never emitted) and validate_body.py:52 enshrines njnt==8 as the M1 gate.
- **The old "pelvis tunnels to −1.44 m" was a METRIC ARTIFACT:** with a welded root, qpos[2] is ankle_L in radians (qposadr 0=hip_L, 1=knee_L, 2=ankle_L) — every "height" printed since run 2 was ankle angle; the body never fell or translated.
- **Second layer (causal A/B):** adding ONLY `<freejoint name="root"/>` on a TEMP copy → the body crumples in <0.5 s (pelvis z→0.100 m, knee_L −44°, heels still never register) — the freejoint is **necessary but not sufficient**; the known M1 deviations gate it: no muscle damping (AnimatLab LinearHill B 400–800 N·s/m vs the 1.5 N·m·s/rad stand-in), rigid tendons, spawn pose OUTSIDE its own limits (qpos0=0 vs ankle range [−20,−5]° → 5° outside at t=0; kN ankle actuators blow soft limits through 18.7–99.1°; femur_L×femur_R grinding 3.8e6 N·steps).
- Neural side verified ALIVE (MV relays −99 → −23..−35 mV loaded; ctrl>0 on knee-ext/ankle-flx/hip pools) — it drives, but never receives the heel trigger and its body cannot use what it drives. **Resume order (Ben's call, not started):** (1) emit the freejoint in make_w2l_mjcf.py + fix validate_body.py's njnt==8 gate and the gate's qpos addressing together; (2) M1 open items #2/#4 become hard prerequisites (settled spawn ~0.99 m; muscle-damping stand-in sized from ΣB·r², not a scalar knob); (3) reconcile ankle joint-zero vs range at spawn; (4) only then re-attempt the 20 s gate vs Li's reference (period 1.305 s, duty ~0.5, antiphase, 0.64 m/s, height 0.95–1.02). Diagnosis written as **section 8 of `reports_20260925\goal2_m2_li_architecture.md`**.

---

## 9. Supervisor-gate verdict

**GO.** Six checks, all executed on easteregg2 by the gate pass; `failures: []`.

1. **Build gates PASS** — `gate_ab_w2lvar.py`: (410, 376, 1186) exact, variant (888, 382, 7430), 23/23 spot-checks, no CIN leak, SNS_NumpyFixedTau; `gate_syn6.py b`: (794, 382, 3033) exact.
2. **Winners real PASS** — all 10 studies: DB argmax over COMPLETE trials == recorded curriculum winner (exact trial number, 0.000e+00 value delta, params to 1e-12), tonight's batches included; air rhythm gates honored (rises 30/20/30/28, REAL RHYTHM verdicts); ground floors hold frozen −320 trials in the bulk, no winner on a floor. One caveat, **not an exploit**: `ee2_harvest_w2lvar_s4.log` prints "exploit gate: FAIL" because the gate compares a fixed-ref repro (−230.860) against an old-ref json (−225.891) — Δ 4.969 is exactly the kine_ref shift; winner t17 is a genuine walker.
3. **Bit-exact replay PASS** — w2lvar s5 trial 4 replayed from the pre-campaign snapshot: −165.33596254795964 vs the tuner log, **Δ +0.000e+00**.
4. **Wiring PASS** — all four editor suites re-run (17/17, 4/4, 9/9, static check); live net rebuilt (95/208/10/12) and smoke reproduces; every named wiring error re-verified from persistent evidence (incl. the clone-bug-is-aproj-only correction via the gate's own XML probe); the split-net stale edges confirmed report-only.
5. **Protected files PASS** — git status carries no entry for `optuna_walk.db`, `spinal_run.npz`, `reports_20260923/`, `reports_20260924/`; mtimes/sizes identical (4,403,200 B @ 22:53:38; 8,238,797 B @ 22:53:39); HEAD `b084a6ab` unchanged throughout — the campaign made NO commits. (Factual note: Ben's own commit `58646a0e` "SNS Simscape" exists on another ref, touches only `Code/Matlab/SNS_Simscape` dev files, is not campaign content.)
6. **Numbers audit re-derived PASS** — independent scipy mat dumps + fresh MATLAB re-evals reproduced every Xi/GoF/margin digit (§5).

---

## 10. For Ben — action list

1. **Commit everything** (GitHub Desktop): the kine_ref fix, `gait2392_robotbody.osim` one-liner, the edited `Biped_2xCPG_wSubs.aproj`, the connectome set (generator + templates + editor HTML), `w2l_cpg` files, chain dbs/npzs — **and the two UNTRACKED mine files** `spinal\w2l_aproj_mine.json` / `li_aproj_mine.json` (until committed, the %TEMP%-loss fix's regeneration contract rests on uncommitted files). All pre-09-26 kine scores shift once kine_ref lands — expect it.
2. **Rule on the pitch-only-origin discovery** (`minimizeFlxPin.m`): fix code and invalidate post-09-13 numbers, or fix docs only. Everything since `8e3c1b9d` inherits pitch-only-origin.
3. **Two tex corrections owed** (I/confirmers made none): "+5.03 %" → "+5.02 %" (30-results.tex:355 + Appendix C narrative); Appendix C bracket-frame line refs → 306/305/326/339/338 (+327). Then compile on EB475WS4 (no TeX on easteregg2) and **sync Overleaf** (muscle-force corrections still local-only).
4. **One AnimatLab GUI session:** fix the RE clone type + the 4 commissural link types in `Walker_2_Layer_CPG_BilateralRG.aproj` (its .asim is already correct); then **re-export `Biped_2xCPG_wSubs_Standalone.asim`** before any RG re-test. Scope the 22 still-missing synapses when ready.
5. **Go/no-go on `w2l_mujoco/build_w2l_split_net.py`** (still hardcodes the invented V3→contra-RG-E edges, lines 119/122 — its 09-25 march results used them).
6. **add_mag3_r:** overrule window open (thumb point applied on your behalf; variants preserved in `add_mag3_compare_20260926\`); route the 0.82× adduction shortfall into Opt_run sizing; decide the ~10× passive-overstretch (P1/TSL/route) finding; schedule regeneration of `gait2392_robotbody_hip.osim`, `gait2327.osim`, and the `mjc\` conversion (one stale path-point each).
7. **w2lvar s5 winner t68 (−144.464)**: near search-space edges (drive 3.82/4.0, ib_to_mn_inh 0.014, contact_onset 0.003) — replay-verify before trusting; the s5 loop stopped at 5 of 10 batches and is resumable if you want more.
8. **syn6:** s4 improved +53.4 like-for-like (−180.845); s5 seed still wins after 40 trials — architectural block confirmed; decide continue vs shelve.
9. **Housekeeping:** repoint `_pf_layer_variant_test.py` gate 3 at stage 4 and re-record post-commit (stale twice over); commit an SNS_SpinalNetwork verify log (the 4.2e-6 mV number has no artifact); fix `muscle_force_compare.csv` header (6 names vs 8 fields); quote W2L smoke antiphase with the env (−0.525 easteregg2 vs −0.532 elsewhere); goal-5 "monotone" → "near-monotone".
10. **Li**, if pursued: resume order in §8 (freejoint + validator/qpos fix → M1 damping/spawn prerequisites → ankle zero-vs-range → 20 s gate).

---

### What this bring-back pass itself ran (verification honesty)

`hostname`; `git log -1` (HEAD b084a6ab); `git status` on kine_ref.py (` M`); directory listings of `reports_20260925\easteregg2\`, `Dig_out\pulley_campaign_20260926\`, `add_mag3_compare_20260926\`, `reports_20260925\tmp\`, `Mesh_Optimization\Results\`; full reads of `ee2_rebase_w2lvar_s5.log`, `ee2_rebase_syn6_s5.log`, `ee2_harvest_w2lvar_s4.log`, `BiPulley_batch_summary.csv`, `launch_record_20260926.txt`, one batch stdout log (`ee2_w2lvar_s5_b3.log`); winner-value greps of the three curriculum JSONs; existence checks for the add_mag3 verdict, LI report and LI gate log. Everything else in this report is quoted from the campaign's logs/claims record and the confirmer/gate passes, attributed inline. No files were modified by this pass except this report and the AGENTS.md appendix; no commits made.
