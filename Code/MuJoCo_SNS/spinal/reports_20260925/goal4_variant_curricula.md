# Goal 4 — Variant Curricula (`w2lvar` + `syn6`): architecture, wiring, 5-stage chain results

**Campaign:** 2026-09-25, EB475WS4, env `myo` (`C:\Users\Ben Bolen\.conda\envs\myo\python.exe`), cwd `Code\MuJoCo_SNS\spinal`.
**Status: COMPLETE — both variant chains finished 5/5 stages; supervisor verdict GO ("the winners are real and reproducible"); this is the final write-up.**
Sources: `goal4_build_w2lvar.md`, `goal4_build_syn6.md`, `goal4_wiring.md`, `goal4_chain_w2lvar.md`, `goal4_chain_syn6.md` (same folder), plus the supervisor audit verdict. Code anchors and every quoted number below were re-checked by the writer against the chain logs, the two variant dbs, and the source reports (see §0 for exactly which checks were re-run here vs attributed to the supervisor).

---

## 0. Verification provenance (what was run, by whom)

**Run by this writer during the write-up (2026-09-25 ~21:05–21:15):**

1. **Read-only db census** — `reports_20260925\audit\audit_db_census.py` (sqlite `mode=ro` on all three dbs). Result: **`NO PROBLEMS: db census matches all claims (counts, winners, jsons).`** All 10 studies present with exactly the claimed COMPLETE counts (w2lvar 20/18/18/23/22; syn6 20/18/18/23/22; every trial `COMPLETE`, 0 failed/running), each recorded winner IS the argmax of its study, and each winner json matches the db bit-exactly (score/trial/params). `optuna_walk.db`: 20 studies, zero `curr_w2lvar_*`/`curr_syn6_*` names (no leak).
2. **Defaults gate re-run** — `reports_20260925\audit\audit_defaults_gate.py` with `AARL_NET` unset: `counts neurons/inputs/synapses = (410, 376, 1186)` → **GATE: PASS**; `params.G['syn6'] default = 0.0`, `params.G['syn6_brainstem'] default = 0.0` (JSON-RULE (a) confirmed in code at `params.py:323,329`).
3. **Protected-file stat** — `spinal_run.npz` 8,238,797 B @ 2026-09-25 01:10:02 AM; `optuna_walk.db` 4,403,200 B @ same (the pre-campaign baseline vintage, unchanged).
4. **Spot replay (syn6 s1)** — `reports_20260925\audit\audit_replay_s1.py syn6` re-run by the writer (log `logs\writer_replay_syn6_s1.log`, npz redirected to an audit-private name; the supervisor's prior evidence npz was additionally backed up first as `audit\audit_replay_syn6_s1.sup.npz` because the script's scratch name collides). Result: **REPLAY PASS, bit-exact** (details in §3.2). Note: the script's AFTER-verdict housekeeping (`shutil.move` of the scratch npz) crashed with WinError 32 (file still held open by the script's own `np.load`), which is why the process exit code was 1 — the PASS verdict and score had already printed; the supervisor's evidence npz was never touched and the scratch npz is preserved as `audit\audit_replay_syn6_s1.writer.npz`.
5. **Code anchors verified first-hand:** variant selector `build_network.py:1020-1031` (`AARL_NET` env routes, lazy import); `params.py:323,329` (`syn6=0.0`, `syn6_brainstem=0.0`); syn6 PF watch branch `runner.py:1148` (`_pf_watch = ("PF_S1_r", "PF_S2_r", "PF_S3_r", "PF_S4_r")`). All 22 chain artifacts exist (`_curriculum_{w2lvar,syn6}.py`, `build_network_{w2lvar,syn6}.py`, `optuna_{w2lvar,syn6}.db`, `spinal_run_{w2lvar,syn6}.npz`, 10 winner jsons).
6. **Independent landscape cross-check:** the census's per-study value lists reproduce every per-stage landscape claim in the chain reports (e.g. w2lvar s1 = 11 real + 9 rhythm-gate floors, syn6 s1 = 12 real + 7 gate + 1 NaN, syn6 s3 = 5 falls, w2lvar s4 = 16 frozen + 5 clip + 2 real, w2lvar s5 = 12 real + 10 frozen, syn6 s5 = 8 real + 14 frozen). Gate arithmetic re-verified: −10 + 0.05·81.1 = −5.9456 and −10 + 0.05·93.0 = −5.3484.

**Attributed to the supervisor's audit (GO verdict, not re-run here unless listed above):** both s1 winner replays bit-exact, per-number cross-check of `goal4_chain_*.md` against logs/db (exit codes, seeded/columns-verified lines, harvest + exploit-check metrics, file sizes/mtimes), and the protected-artifact before/after sweep (only post-20:44 spinal-root write = `spinal_run_w2lvar.npz` 20:47 from the campaign's own stage-5 harvest repro). The verdict text as transmitted ends "Minor notes below are low-severity; none blocks the write-up" — the minor notes themselves were not included in the copy given to this writer; the GO is unconditional.

---

## 1. Architecture of the two variants

Both variants were built as NEW builder modules selected by a shared, default-inert selector at the top of `build_network.py` `build()` (lines 1018–1031: `AARL_NET=w2lvar` / `AARL_NET=syn6` or `G["syn6"]>0`; lazy import — with the env unset and the key 0 the standard path never touches variant code). Gating evidence: defaults build stays exactly (410, 376, 1186) after every edit (gate (a), re-run 5× across the campaign and again by the supervisor and by this writer).

| | **w2lvar** | **syn6** |
|---|---|---|
| Builder | `spinal\build_network_w2lvar.py` (`W2LVarnet`) | `spinal\build_network_syn6.py` (`Syn6Network` subclasses stock `SpinalNetwork`) |
| Selection | env `AARL_NET=w2lvar` | env `AARL_NET=syn6` **or** `G["syn6"]>0` (`params.py:323`) |
| Net size (variant build) | 888 neurons / 382 inputs / 7430 synapses | 794 neurons / 382 inputs / 3033 synapses |
| Layout idea | Walker_2_Layer_CPG (Ben's AnimatLab biped) transplanted onto the 92-actuator gait2392 body: 1 RG/side + TWO joint-layer PF pairs (hip; ONE knee+ankle **synergy**) | Data-driven 6-synergy walker: 1 RG/side + SIX synergy PF layers `PF_S1..S6_{r,l}`, one per channel of `synergy_basis.npz` (W [43×6], rank-6 NMF of back-solved activations) |
| RG | verbatim parent `_build_rg` mechanics: persistent-Na half-centers (`NonSpikingNeuronWithPersistentSodiumChannel`, fixed tau_h = 0.35 s via `SNS_NumpyFixedTau`) + InE/InF laminated mutual inhibition (g 4.0) + DRIVE→E/F + POSTURE→E | same stock `_build_rg` pattern, verbatim |
| PF→MN mapping | merged cells drive per-muscle MNs: KNEE-E → knee_ext + ankle_pf, KNEE-F → knee_flex + ankle_df (W2L biarticular-gas pattern); weights from the parent's fitted `joint_pf_weights.json` with anatomical fallback | measured W through the audited Szczecinski-2017 Eq-18 `analytical_conductance` imported from `fsa_backsolve.py:132-137` (`g = k·R·Gm/(ΔE − k·R)`, R=5 mV, Gm=1 µS, ΔE=8 mV); edges only where W>0 (conditional topology); 152 edges (l)/159 (r); the 6 runner-pruned actuators keep MN+reflex arc but get NO synergy drive |
| Per-muscle motifs | 92×: MN + Ia/II/Ib encoder ports, IaIN, IIX, IIIN, IBIN, RC (Renshaw), IBEXC stance reversal, KINH available at gain 0 | Shinohara autogenic motif per muscle (Ia→MN/IaIN/reciprocal; II→IIX/IIIN; Ib→IBIN) — proprioception autogenic per muscle, never group-broadcast (Ben's reading rule); NO Renshaw/KINH/AFF loops (renshaw forced back to 0 after `OW.set_params`) |
| Rule-file dress (both) | heel stance-reset at PF layer, toe→dorsiflexion inhibition, Ib load, Shevtsova V2a/V0V/V0D/V3-E commissurals (replaces the parent C1/V3 block), Shinohara afferents, master-rules crossed heel edge (w2lvar) | heel→InE/InF/PF_IN_E; toe→TOEDF 5.0→inh the dorsiflexion channel (S5); per-stance-muscle Ib→IBEXC; Shevtsova commissurals per side |
| Measured degeneracy handling | — | symmetrized families `[F,E,E,F,F,E]` both sides; tau stagger within family; MIXED channel S5 (stance_frac 0.41) excluded from lamination and given BOTH RG drives (final form: S5 distinct, |corr| ≤ 0.51; within-family channels remain near-degenerate — documented residual) |
| Known deviations | D1–D6 documented in `goal4_build_w2lvar.md` (headlines: RG→PF drive 2.4 not the drawing's SynAmp-scale 0.1; per-muscle afferent central gains = aggregate ÷ family count — fixed a measured co-latch; Shevtsova antiphase/sync rebalance — fixed measured in-phase stepping; toe-inhibition lands on the merged KNEE-F cell, inherent to the synergy layout) | D1–D3 documented in `goal4_build_syn6.md` (S3K-tuned knobs for RG→PF/lamination not the 0.1/2.749 drawing scale; Shevtsova IniE/IniF folded into InE/InF; brainstem gamma/alpha fold gated behind `syn6_brainstem`, default 0, after it E-latched the net) |

**Honest framing (per the 2026-09-25 early-AM correction in AGENTS.md):** these variants take Ben's rule files (`ben_rules_20260924.json` 90n/99e, `ben_shinohara_20260924.json`, `ben_shevtsova_20260924.json`) as the *dress* over two architecture families, with documented deviations — they are NOT a rebuild of Ben's per-micro-layer connectome drawing, and must not be described as "following his connectome rules."

**Builder-gate outcomes** (from the build reports): w2lvar gates (a)/(b)/(c) all PASS — gate (c) rhythm smoke passed after two documented tuning rounds (co-latch → per-muscle afferent division; L/R in-phase → Shevtsova rebalance; final: RG E/F antiphase r = −0.910/−0.916, period 1.80 s, E-duty 0.69/0.71, `logs\gate_c_smoke3.log`). syn6: gates (a)/(b) PASS and network-only rhythm (c1) PASS (period 2.133 s, E/F corr −0.989), but the full-runner air rhythm (c2) was **PARTIAL** (strict interleave `FEFEFFE`, slow ~2–3 s events that die after ~2 cycles; ≥5-peak criterion FAILed). Two E-latch mechanisms were found and gated off (brainstem fold → `syn6_brainstem`; un-gated force-proportional Ib→LBIN → gated on `G["ib_rge"]>0`); the residual slowness is documented as a tuning surface, not a builder defect.

## 2. Curriculum design (both chains)

Wired in NEW self-contained drivers `spinal\_curriculum_w2lvar.py` / `spinal\_curriculum_syn6.py` (per `goal4_wiring.md`): each pins its own env (`AARL_NET`, chain-private `AARL_NPZ` → `spinal_run_{variant}.npz`) as the first act of `main()`, uses its OWN sqlite db (`optuna_w2lvar.db` / `optuna_syn6.db` — the stock `_curriculum.py` db path is hardcoded), and never touches `optuna_walk.db`, `spinal_run.npz`, or any default-path behavior.

**Ben's 5 stages** (numbering per his "standing = stage 3, walking = stage 4" ruling):

| stage | study (w2lvar / syn6) | run | objective |
|---|---|---|---|
| 1 | `curr_w2lvar_s1_air_deaff` / `curr_syn6_s1_air_deaff` | `--no-ground --no-afferents --no-interleg --time 14` | air: `3·rises + 0.5·(−knee_min)` over t∈[5,17], rhythm-gated |
| 2 | `..._s2_air_aff` | `--no-ground --time 14` (afferents + interleg ON) | same air objective |
| 3 | `..._s3_balance` | `--stand-eval 8 --rig-scale S` | standing: `100 − 400·sway − tilt_max − 40·|0.5−contact_sym|`, NaN −200 / fall −150 |
| 4 | `..._s4_walk_nocontact` | `--no-ground --eval` (contact OFF) | walking kine `max(kine_score, −315)` vs the reference cycle in air (documented interpretation: pattern-match before facing contact; kz floor is stage-5 only) |
| 5 | `..._s5_walk_contact` | `--eval` (normal ground) | kine + kz<0.62 (−20) + tilt>40° (−10) penalties |

- **Trials:** defaults 18/18/18/22/22 (`optimize(n_trials=N)` INCLUDES the enqueued seed). CLI: `_curriculum_{variant}.py <stage> [n]`.
- **Sentinels** (stock discipline preserved): NaN −200 (air) / −400 (walk), fall −150, rhythm gate `rises<3 → −10 + 0.05·(−knee_min)`, frozen-walk −320, kine clip −315, flutter >30 rises −200, unphysical RoM −200.
- **Seeding:** FULL dicts (every searched key explicit — `enqueue_trial` samples anything missing); chaining via `curriculum_{variant}_stage{N}.json`, falling back to stage N−1. **Concurrency:** the two chains run concurrently safely (separate dbs + npz), proven by the concurrent stage-1 smokes; smoke trials remain as the first trials of each study (`load_if_exists=True` resume).
- **Search spaces** — only knobs each variant actually consumes, verified by `tmp\gkey_scan.py` regex over both builders + runner (not assumed):
  - **w2lvar:** S1 stock KEYS1 (`drive, rg_nap_h, desc_e, desc_f, rg_to_pf`); S2 + `ia_to_mn/ia_to_antagonist/ii_to_mn/ib_to_mn_inh` seeded AT the params defaults 0.6/0.4/0.4/0.35 (documented deviation from the [0,0.5] new-gain rule — these are pre-existing keys the variant's motifs read; seeding 0 would change stage-1 behavior); S4 + `rg_mutual_inh [2,6], rg_weak_exc [0,0.5], pf_to_mn [0.5,3]`; S5 + `contact_onset [0,1], contra_swing [0,1.5], ib_group_exc [0,0.8], ib_exc_to_mn [0,0.9]`.
  - **syn6:** S1 KEYS1; S2 + `heel_rge [0,0.5], toe_rge [0,0.5]`; S4 + `rg_mutual_inh, rg_weak_exc, pf_recip_inh [1,5]`; S5 + `contact_onset, contra_swing, ib_rge [0,0.5]` (build-time gate for the E-latch fix).
  - Excluded per variant: keys the builder's `VARIANT_G` overlay clobbers (w2lvar), keys the variant never reads (`f1_*` for syn6, `c1/v3/*_central`, `vest_*` — both build no VEST cells).
  - `G["syn6"]=1.0` pinned in every syn6 merged dict + seed (JSON RULE); `syn6_brainstem=0.0` pinned (the measured E-latch); w2lvar `renshaw` 0.5 matches `VARIANT_G`.
  - No s3k BASE merge (stock stage-4/5 move): s3k tuned the STOCK topology; the variant chains are self-contained (v10-multiplier shell + their own stage winners).

## 3. Per-variant, per-stage outcomes

All scores bit-verified two ways: (i) db-vs-json census (writer, §0) and (ii) exact objective-path replay (supervisor for both s1 winners; chain harvest scripts for every stage; writer additionally re-ran the syn6 s1 replay, below). **The 09-20 static-pose exploit pattern is absent everywhere.**

### 3.1 w2lvar (`optuna_w2lvar.db`) — CHAIN COMPLETE 5/5

| stage | study | trials (COMPLETE) | winner | score | exploit check |
|---|---|---|---|---|---|
| 1 air deaff | `curr_w2lvar_s1_air_deaff` | 20/20 (2 smoke + 18) | trial 17 | **80.836** | **PASS** — replay bit-exact (80.835587 vs 80.835587, delta 0.000000): rises 14, knee_min −77.67°, RG-E span 9.75 mV, 1.43 Hz, E-duty 0.45. 0 NaN, 9 gate floors, 11 real rhythms |
| 2 air aff | `curr_w2lvar_s2_air_aff` | 18/18 | trial 13 | **60.230** | **PASS** — bit-exact (60.230493 == 60.230493): rises 10, knee −60.46°, span 7.88 mV. Seed (stage-1 winner + afferents ON) LOST the rhythm (−9.804); TPE re-found one; reflex gains pushed DOWN (ii 0.048, ib 0.151, ia 0.356) |
| 3 balance | `curr_w2lvar_s3_balance` | 18/18 | trial 0 (random sample) | **25.077** | **PASS (genuine stander)** — bit-exact: bal_fell false, sway 0.0996 m, tilt_max 25.32°, contact_sym 0.2562. 1 fall, 17 standers. TPE could not beat the first random sample in 17 tries → standing map under-explored |
| 4 walk nocontact | `curr_w2lvar_s4_walk_nocontact` | 23/23 (1 smoke + 22) | trial 17 | **−225.891** | **PASS (real cycle pattern)** — bit-exact: kine −225.891, tilt 1.07°, duty 0.6404. Landscape: 16× frozen −320, 5× clip −315, only 2 real kine scores — air pattern-matching is HARD here, but the −320 floors prove motionlessness orders far below |
| 5 walk contact | `curr_w2lvar_s5_walk_contact` | 22/22 | trial 4 | **−165.041** | **PASS (genuine upright ground walk)** — bit-exact: kine −165.041 with ZERO penalties, kz 0.8297 (floor 0.62), tilt 27.72° (<40), duty 0.4597, contact_frac real. 12/22 real kine scores (vs 2/23 in air — ground contact HELPS this variant), 10 frozen. TPE adopted nonzero contact-pathway gains: contact_onset 0.134, contra_swing 0.298, ib_group_exc 0.087, ib_exc_to_mn 0.570; drive up 1.68→3.67 |

### 3.2 syn6 (`optuna_syn6.db`) — CHAIN COMPLETE 5/5

| stage | study | trials (COMPLETE) | winner | score | exploit/validity check |
|---|---|---|---|---|---|
| 1 air deaff | `curr_syn6_s1_air_deaff` | 20/20 (2 smoke + 18) | trial 17 | **95.639** | **PASS (real rhythm)** — bit-exact: rises 16 (gate ≥3), knee −95.3..+10.3°, span 9.74 mV, 1.57 Hz, E-duty 0.45. 12 real + 7 gate + **1 NaN (−200, trial 18)** — far below the winner, no inversion |
| 2 air aff | `curr_syn6_s2_air_aff` | 18/18 | trial 12 | **86.751** | **PASS (real rhythm)** — bit-exact: rises 21, span 9.78 mV, 2.06 Hz; knee swing shallowed to −47.5..−12.3° (trades flexion depth for burst count). heel_rge 0.011 / toe_rge 0.043 ≈ 0 — contact ports add nothing in air (same "tolerated, not additive" outcome as the stock s2b) |
| 3 balance | `curr_syn6_s3_balance` | 18/18 | trial 11 | **37.463** | **PASS (genuine stand)** — bit-exact: 100−400·sway−tilt−40·|0.5−sym| = 37.4630 == json; bal_fell false, sway max 0.0878 / rms 0.0610, tilt 18.76°, contact_sym 0.284, kz 0.903. rig_scale 0.946 (barely weaned); the mechanosensor gains dead in air WON here (heel 0.417 / toe 0.477) |
| 4 walk nocontact | `curr_syn6_s4_walk_nocontact` | 23/23 (resumed over 1 smoke) | trial 14 | **−237.253** | **PASS (genuine air-gait)** — bit-exact: bilateral=True, 3 cycles/leg, period_cv 0.0004; RG rises 3, knee −49.6..+10.1°. Pattern is real but SLOW/in-phase: T 2.45 s vs ref 1.23, lag_rl 0.0 (legs IN-PHASE vs ref 0.51), duty 0.41 vs 0.61, ankle PF-biased (−53.5° vs +4.2°). 10 real + 13 frozen. `rg_weak_exc` 0.177 > 0 — the conditional RG↔RG weak-excitation topology ENGAGED |
| 5 walk contact | `curr_syn6_s5_walk_contact` | 22/22 | **trial 0 = the SEED** | **−197.215** | **PASS (genuine ground walk)** — bit-exact: contact_frac_r 0.955 / _l 1.0 (real loading), right leg n_cycles 4 at T 1.88 s, lag_rl 0.395, kz 0.840, tilt 28.07°, zero penalties. 8 real + 14 frozen. **"Seed wins" outcome:** best searched trial only −261.95; the three stage-5 contact pathways never earned nonzero values (contact_onset 0.0, contra_swing 0.0, ib_rge 0.0). Left leg planted ~100% (bilateral=False) |

**Writer's own spot replay (check 4, §0):** `audit_replay_s1.py syn6` re-run through the exact stage-1 objective path (same BASE_MUL merge, same runner args `--no-ground --no-afferents --no-interleg --time 14 --drive 3.393724708721609`, audit-private npz) — outcome: **REPLAY PASS, recomputed 95.6388239983187 vs json 95.6388239983187, delta +0.000000 (0.0000%), rises 16, knee_min −95.28°, span 9.7424 mV** (log `logs\writer_replay_syn6_s1.log`: `finite=True rom_ok=True rises=16 span=9.7424 knee_min=-95.28 knee_max=10.29`), matching the supervisor's independent replay and the chain's exploit check bit-for-bit. The process exit code was 1 solely due to the AFTER-verdict `shutil.move` housekeeping crash (WinError 32, npz held open by `np.load`) — the PASS verdict had printed; supervisor evidence untouched, writer npz preserved as `audit\audit_replay_syn6_s1.writer.npz`. (The w2lvar s1 replay stands on the supervisor's run + the chain harvest — not re-run by this writer.)

## 4. Supervisor verdict and findings

**Verdict: GO** — "Audited the goal-4 variant curricula (w2lvar + syn6) before write-up; verdict GO — the winners are real and reproducible." Findings as transmitted:

1. **Trial counts / winner integrity:** both dbs read READ-ONLY (`mode=ro`) via `audit_db_census.py` — all 10 studies present with exactly the claimed COMPLETE counts (w2lvar 20/18/18/23/22; syn6 20/18/18/23/22; ALL trials COMPLETE, 0 failed/running); winner jsons match the db bit-exactly (score/trial/params); each recorded winner IS the argmax of its study; `optuna_walk.db` untouched (20 studies, zero `curr_w2lvar_*`/`curr_syn6_*` names). *(Independently re-confirmed by this writer, §0 check 1.)*
2. **No sentinel/exploit winners:** air winners have rises 14/10 (w2lvar s1/s2) and 16/21 (syn6 s1/s2) vs the ≥3 gate, scoring 60–96 vs a gate ceiling of −3; all four kine winners (−225.89 / −165.04 / −237.25 / −197.21) are genuine cycle-detected scores with the frozen floor at −320 ×16/10/13/14 proving ordering; standing winners are genuine stands (bal_fell False, sway 0.088–0.100 m); the one NaN (−200, syn6 s1 trial 18) sits far below its winner. The 09-20 static-pose exploit pattern is absent.
3. **Defaults gate re-run by the supervisor** (and again by this writer): `audit_defaults_gate.py` → (410, 376, 1186) PASS with `AARL_NET` unset; `params.G` `syn6`/`syn6_brainstem` default 0.0.
4. **Spot replays (supervisor):** one cheap winner per variant via the exact objective code path (`audit_replay_s1.py`, AARL_NPZ redirected to audit-private npz): w2lvar s1 80.83558669459252 vs json 80.83558669459252 (delta +0.000000, 0.0000%; rises 14, knee −77.67°) and syn6 s1 95.6388239983187 vs 95.6388239983187 (delta +0.000000; rises 16, knee −95.28°) — both bit-exact. *(syn6 s1 additionally reproduced by this writer, §3.2.)*
5. **Full cross-check of every quoted number** in `goal4_chain_w2lvar.md` and `goal4_chain_syn6.md` against logs/db: value lists match the db exactly, all 10 exit codes = 0, seeded lines / 'columns verified' lines / harvest + exploit-check metrics (rises, spans, knee ranges, duty, kz, tilt, contact_frac, T_r, lag_rl, sway, contact_sym), file sizes/mtimes, gate arithmetic (−5.9456 = −10+0.05·81.1; −5.3484 = −10+0.05·93.0), and the gate-c smoke metrics all verified.

**Protected artifacts** (supervisor sweep + writer stat): `spinal_run.npz` 8,238,797 B @ 2026-09-25 01:10:02 AM before/after; `optuna_walk.db` same vintage; `reports_20260923\` / `reports_20260924\` untouched during the goal-4 window; the only post-20:44 spinal-root write is `spinal_run_w2lvar.npz` 20:47 from the campaign's own stage-5 harvest repro. Minor notes referenced by the verdict are low-severity and none blocks the write-up (notes themselves not transmitted to this writer — see §0).

## 5. Honest assessment — what the new architectures do and do NOT yet do vs s3k

**What they DO (all writer-verified against the census/logs):**
- Both are complete, reproducible chains: 10/10 studies, 202 trials total, every trial COMPLETE, every winner bit-exactly replayable, zero exploits, zero leaks into protected stores.
- Both self-sustain deafferented air rhythms with deep knee swing (w2lvar 1.43 Hz / knee −78°; syn6 1.57 Hz / knee −95°) and survive the afferented-air stage.
- Both produce GENUINE ground walks with real foot loading, upright (kz 0.83–0.84, far above the 0.62 floor), moderate tilt (~28°), and zero penalty deductions at stage 5.
- w2lvar's contact pathways were EARNED by TPE at stage 5 (4 nonzero gains adopted, drive up 1.68→3.67) — unlike most new-pathway knobs in this lab's history, which tend to end at 0.
- syn6's stage-4 winner engaged the conditional `rg_weak_exc` topology (>0), and its stage-3 win confirms contact ports matter once feet load (heel 0.42 / toe 0.48 after ≈0 in air).

**What they do NOT yet do (context caveats, not hidden):**
- **Not apples-to-apples with s3k.** The variant kine scores (−165.0 w2lvar / −197.2 syn6) land numerically in the band the stock s3k walker scores against real gait references (−189.5…−241, goal-5 reports), but the chain reports flag these as **context only, different reference/scoring routes** — no goal-5-style scoring against subject01/Falisse/Ong has been run for the variants. No head-to-head claim is made.
- **Cadence/duty/phase gaps are the same ARCHITECTURAL family the stock chain documented.** w2lvar ground: duty 0.46 vs ref 0.61, tilt 27.7° high. syn6 ground: heavy double support (duty 0.745), left leg planted ~100% (no clean left cycle — the obvious next architecture problem, same family as the stock interleg/phase gaps), knee flexion shallow (−49.7 vs −69.7°), and its stage-4 air pattern was slow (T 2.45 s) and IN-PHASE (lag_rl 0.0). Neither has attacked the sensory-phase-reset piece the stock chain identified as the architecture-level bottleneck.
- **syn6's full-runner rhythm was PARTIAL at the builder gate** (slow events that die after ~2 cycles); the chain tuned around it, but its rhythm fragility remains a known tuning surface (DRIVE, `rg_weak_exc`, `rg_nap_h`, POSTURE→RG_E).
- **Standing is weak in both:** w2lvar 25.1 (winner = a random sample; 17 TPE tries could not beat it — under-explored map) and syn6 37.5 (at rig_scale 0.946 — barely weaned; the search keeps rig high because lower rig directly lowers score). Neither approached the weaning boundary from below (stock chain support boundary S≈0.8–1.0).
- **Search efficiency at stages 4–5 is low** (sentinel-heavy landscapes: 16/23 and 13/23 frozen in air; 10/22 and 14/22 frozen on ground) — expected for first contact-free/ground kine evaluations of brand-new networks, but it means the winners may not be near the achievable optimum.
- **syn6 stage-5 was won by the SEED** — 21 TPE trials could not improve the transplanted stage-4 config on the ground; more stage-5 trials (or re-seeding stage 4 from stage 2, documented as a legitimate correction) are the cheap next levers.
- **Within-family synergy channels in syn6 remain near-degenerate** in constant-drive air (S2/S3/S6 corr +0.98…+0.997); per-channel delay interneurons would be the fix if 6 visually distinct traces are wanted — not built.
- **Neither winner has been basin-gated** (`basin_gate.py` is tuned to the stock chain's dump state; the Simulink-realization bifurcation caution in AGENTS.md applies to any future port).
- **w2lvar's merged knee+ankle synergy loses knee-vs-ankle phase differentiation by design** (D4); the 3-stack layout is where that edge is exact, if Ben wants it back.
- Cosmetic known issue: `runner.py`'s final matplotlib panel hardcodes stock phase-cell names → "(plot skipped)" at the end of a variant run; the npz itself carries all watch columns.

## 6. Exact resume / extension commands

**All 10 studies are COMPLETE — no stage is unfinished.** The chains are resumable/extensible: each stage seeds from `curriculum_{variant}_stage{N}.json` (falling back to stage N−1) and `load_if_exists=True` means re-running a stage ADDS n trials to the existing study (the smoke trials remain as the first trials; seed = trial 0 everywhere). The two chains may run concurrently (different dbs, different npz — proven). Each script self-pins `AARL_NET` + `AARL_NPZ` as the first act of `main()`, so each stage is ONE detached command:

```
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
C:\Users\Ben Bolen\.conda\envs\myo\python.exe _curriculum_w2lvar.py <stage> [n]
C:\Users\Ben Bolen\.conda\envs\myo\python.exe _curriculum_syn6.py   <stage> [n]
```

stages 1..5, defaults 18/18/18/22/22. Highest-value extensions implied by the outcomes: `... _curriculum_w2lvar.py 3 18` (re-tune the under-explored standing map), `... _curriculum_syn6.py 5 22` (seed won; give TPE more shots / consider re-seeding stage 4 from stage 2), and stage-5 continuation for both (sentinel-heavy landscapes leave the optimum unexplored).

## 7. Artifact index (all under `Code\MuJoCo_SNS\spinal\` unless noted)

- Builders: `build_network_w2lvar.py`, `build_network_syn6.py`; shared selector `build_network.py:1018-1031`; `params.py:323,329`.
- Curriculum drivers: `_curriculum_w2lvar.py`, `_curriculum_syn6.py`; dbs `optuna_w2lvar.db` / `optuna_syn6.db`; chain npz `spinal_run_w2lvar.npz` / `spinal_run_syn6.npz`; winner jsons `curriculum_{w2lvar,syn6}_stage{1..5}.json`.
- Run + exploit-check logs: `reports_20260925\logs\curr_{w2lvar,syn6}_s{1..5}.{log,done}`, `logs\harvest_w2lvar_s{1..5}.log`, `logs\curr_syn6_s{1..5}_exploit_check.log`, smokes `logs\smoke_{w2lvar,syn6}_s{1,4}.log`, builder gates `logs\gate_c_smoke{,2,3}.log` + `logs\syn6_net_smoke.npy` / `logs\syn6_air_smoke.{npz,log}`.
- Audit (this campaign): `reports_20260925\audit\{audit_db_census.py, audit_defaults_gate.py, audit_replay_s1.py, audit_logs_scan.py, defaults_gate_output.log, replay_{w2lvar,syn6}_s1.log, audit_replay_{w2lvar,syn6}_s1.npz, audit_replay_syn6_s1.sup.npz}`; writer's replay log `reports_20260925\logs\writer_replay_syn6_s1.log`.
- Source reports: `reports_20260925\goal4_build_w2lvar.md`, `goal4_build_syn6.md`, `goal4_wiring.md`, `goal4_chain_w2lvar.md`, `goal4_chain_syn6.md`.
