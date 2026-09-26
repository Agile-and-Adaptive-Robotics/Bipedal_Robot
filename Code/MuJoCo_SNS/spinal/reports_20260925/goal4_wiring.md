# goal4 — VARIANT CURRICULA WIRED + SMOKED (`w2lvar` / `syn6`)

**Status: COMPLETE (all requested smoke gates PASS).** 2026-09-25 13:56–14:30
EB475WS4, env `myo` (`C:\Users\Ben Bolen\.conda\envs\myo\python.exe`), cwd
`Code\MuJoCo_SNS\spinal`. Everything below was executed this session; commands
and outputs are quoted.

## Files delivered (both NEW; no protected file touched)

| file | role |
|---|---|
| `spinal\_curriculum_w2lvar.py` | w2lvar 5-stage curriculum (fresh studies, own db/npz, self-contained env pin) |
| `spinal\_curriculum_syn6.py` | syn6 5-stage curriculum (same shape) |
| `reports_20260925\tmp\verify_smoke_goal4.py` | post-smoke verification (studies, jsons, npz columns, protected files) |
| `reports_20260925\tmp\gkey_scan.py` | G-key/TAU-key scan used to pick the per-variant search spaces |
| `reports_20260925\logs\smoke_{w2lvar,syn6}_s1.log` | stage-1 smoke logs (2 trials each) |
| `reports_20260925\logs\smoke_{w2lvar,syn6}_s4.log` | stage-4 route-proof logs (1 trial each) |
| `optuna_w2lvar.db`, `optuna_syn6.db` | the chains' OWN sqlite dbs (fresh this session) |
| `spinal_run_w2lvar.npz`, `spinal_run_syn6.npz` | chain-private npz (via `AARL_NPZ`) |
| `curriculum_{w2lvar,syn6}_stage{1,4}.json` | stage-winner chaining jsons (stage 1 + the stage-4 smoke) |

`build_network.py`, `runner.py`, `params.py`, `spinal_run.npz`,
`optuna_walk.db` were NOT modified by this task (the variant selector they
use was merged by the earlier goal-4 build tasks; gate (a) re-proven below).

## Ben's 5 stages (both chains)

| stage | study name (w2lvar / syn6) | run | objective |
|---|---|---|---|
| 1 | `curr_w2lvar_s1_air_deaff` / `curr_syn6_s1_air_deaff` | `--no-ground --no-afferents --no-interleg --time 14` | air: 3·rises + 0.5·(−knee_min), rhythm-gated |
| 2 | `curr_w2lvar_s2_air_aff` / `curr_syn6_s2_air_aff` | `--no-ground --time 14` (afferents + interleg ON) | same air objective |
| 3 | `curr_w2lvar_s3_balance` / `curr_syn6_s3_balance` | `--stand-eval 8 --rig-scale S` | standing: sway/tilt/symmetry, fall sentinel |
| 4 | `curr_w2lvar_s4_walk_nocontact` / `curr_syn6_s4_walk_nocontact` | `--no-ground --eval` (CONTACT DISABLED) | **WALKING kine objective** vs the reference cycle (documented interpretation: pattern-match the gait in air before facing contact); kz floor penalty is stage-5 only (the air rig holds the pelvis up) |
| 5 | `curr_w2lvar_s5_walk_contact` / `curr_syn6_s5_walk_contact` | `--eval` (normal ground) | walking kine + kz + tilt penalties |

Trial defaults 18/18/18/22/22; CLI `... _curriculum_w2lvar.py <stage> [n]`.
Sentinels preserved exactly as the stock discipline: NaN −200 (air) / −400
(walk), fall −150, rhythm gate rises<3 → −10 + 0.05·(−knee_min), frozen-walk
−320, kine clip −315, flutter >30 rises −200, unphysical RoM −200. Seeds are
FULL dicts (every searched key explicit — `enqueue_trial` samples anything
missing); chaining via `curriculum_{variant}_stage{N}.json`.

## Search spaces — only knobs each variant actually consumes

Verified by `tmp\gkey_scan.py` (regex over `G["…"]`/`TAU["…"]` in both
builders + runner), not assumed:

- **w2lvar** (`build_network_w2lvar.py` reads: descend_to_rg_e/f, rg_to_pf,
  rg_mutual_inh, rg_weak_exc, posture_to_*, pf_to_mn, ia_to_mn,
  ia_to_antagonist, ii_to_mn, ib_to_mn_inh, ib_group_exc, ib_exc_to_mn,
  f1_*, heel_rge, toe_rge, ib_rge, joint_pf, renshaw; TAU rg_nap_h):
  - S1: `drive, rg_nap_h, desc_e, desc_f, rg_to_pf` (stock KEYS1 ranges).
  - S2: + `ia_to_mn` [0,0.9], `ia_to_antagonist` [0,0.6], `ii_to_mn`
    [0,0.6], `ib_to_mn_inh` [0,0.5] — these are PRE-EXISTING params keys
    with defaults 0.6/0.4/0.4/0.35 that the variant's per-muscle motifs
    read, so they seed AT the defaults (seeding 0 would change stage-1
    behavior — not new default-0 gains; documented deviation from the
    [0,0.5] new-gain rule).
  - S4: + `rg_mutual_inh` [2,6], `rg_weak_exc` [0,0.5] (0 = conditional
    topology absent), `pf_to_mn` [0.5,3].
  - S5: + `contact_onset` [0,1], `contra_swing` [0,1.5], `ib_group_exc`
    [0,0.8], `ib_exc_to_mn` [0,0.9].
  - NOT searched (the builder's `VARIANT_G` overlay at
    build_network_w2lvar.py:128 CLOBBERS them every build): heel_rge,
    toe_rge, ib_rge, joint_pf, ia_in, renshaw, full_rules,
    f1_kneext_inh, f1_anklepf_inh. Also excluded: c1_gain/v3_gain/
    *_central (stock-builder keys the variant never reads),
    vest_ext/vest_flex_inh (variant builds no VEST cells,
    build_network_w2lvar.py:188 `self.vest = False`).
- **syn6** (`build_network_syn6.py` reads: descend_to_rg_e/f, rg_to_pf,
  rg_mutual_inh, rg_weak_exc, pf_recip_inh, posture_to_*, renshaw,
  ib_rge, syn6, syn6_brainstem):
  - S1: stock KEYS1. S2: + `heel_rge` [0,0.5], `toe_rge` [0,0.5]
    (runner-side port scaling; syn6 dress edges exist, stance_fb=True).
  - S4: + `rg_mutual_inh` [2,6], `rg_weak_exc` [0,0.5],
    `pf_recip_inh` [1,5].
  - S5: + `contact_onset` [0,1], `contra_swing` [0,1.5], `ib_rge` [0,0.5]
    (gates the variant's per-muscle Ib→LBIN pathway at build time —
    the gate-c2 E-latch fix).
  - Excluded: f1_* (syn6 builds no KINH), c1/v3/central (unread),
    vest_ext/vest_flex_inh (`self.vest=False`, build_network_syn6.py:294).
  - `G["syn6"]=1.0` pinned in EVERY merged param dict + seed (G-key
    selector, JSON RULE) and `G["syn6_brainstem"]=0.0` pinned (the
    measured E-latch, goal4_build_syn6.md §4.1). `renshaw` forced back
    to 0.0 after `OW.set_params` (which hard-sets 0.5): the syn6 dress
    list has NO RC population — leaving 0.5 would silently add RCs to
    every trial and change the gate-verified topology (794/382/3033).

Both scripts: `renshaw` note for w2lvar = `OW.set_params`'s hard 0.5
matches `VARIANT_G`. No s3k BASE merge (stock stage-4/5 move): s3k tuned
the STOCK topology; the variant chains are self-contained (v10-multiplier
shell that `OW.set_params` requires + their own stage winners).

## SMOKE EVIDENCE (all run this session)

### Stage-1 smoke (2 trials each, real end-to-end runs)

Commands:
```
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum_w2lvar.py 1 2
    > reports_20260925\logs\smoke_w2lvar_s1.log 2>&1        (exit 0)
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum_syn6.py 1 2
    > reports_20260925\logs\smoke_syn6_s1.log 2>&1          (exit 0)
```
(both launched concurrently — the concurrency itself is part of the test)

1. **Switch selects the variant net** — proven by the npz's OWN
   `neuro_names` (the watch tuple is chosen by the net branch at
   runner.py:1145-1153 and saved at runner.py:1637):
   - w2lvar npz: `['DRIVE','POSTURE','RG_E_r','RG_F_r','RG_E_l','RG_F_l',
     'PF_HIP-E_r','PF_KNEE-E_r','PF_KNEE-F_r','PF_ANK-F_r','BAL_*']`
     (w2lvar-merged-cell names only exist in the w2lvar net; the stock net
     would show `PF_E1_r…`, and runner would KeyError on `PF_HIP-E_r`
     otherwise) — `cfg: air`.
   - syn6 npz: `[..., 'PF_S1_r','PF_S2_r','PF_S3_r','PF_S4_r', 'BAL_*']`
     — `cfg: air`.
2. **Objective reads THIS net's columns (verified, not assumed)** — the
   objective derives indices from the npz arrays at every trial and
   prints them once per process (log lines):
   ```
   [w2lvar] npz columns verified: neuro_names=[…PF_HIP-E_r…] -> RG_E_r idx 2; key_joints -> knee_angle_r idx 4
   [syn6]   npz columns verified: neuro_names=[…PF_S1_r…]    -> RG_E_r idx 2; key_joints -> knee_angle_r idx 4
   ```
3. **Study + json write** — `optuna.load_study` on each fresh db:
   `curr_w2lvar_s1_air_deaff: 2 trials COMPLETE`,
   `curr_syn6_s1_air_deaff: 2 trials COMPLETE`; stage jsons
   `curriculum_w2lvar_stage1.json` (score 56.002, trial 0, db/npz/net
   fields recorded) and `curriculum_syn6_stage1.json` (61.908).
4. **No sentinel inversion** — seed trial 0 (the gate-c-proven rhythm
   point) scored **56.002 (w2lvar)** / **61.908 (syn6)**; trial 1
   (TPE startup sample, drive 1.78) lost the rhythm and hit the gate:
   −5.9456 = −10 + 0.05·81.1 (w2lvar knee min) and −5.3484 = −10 +
   0.05·93.0 (syn6 knee min), arithmetic re-verified by command. The
   static/rhythmless config orders far BELOW the genuine rhythm;
   −10-floor < seed; NaN −200 and flutter −200 never approached.

### Stage-4 route proof (1 trial each — the NEW objective route)

Commands:
```
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum_w2lvar.py 4 1
    > reports_20260925\logs\smoke_w2lvar_s4.log 2>&1        (exit 0)
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum_syn6.py 4 1
    > reports_20260925\logs\smoke_syn6_s4.log 2>&1          (exit 0)
```
`--no-ground --eval` + walking-kine ran end to end for both variants; the
full-dict seeds printed in the logs (w2lvar 12 keys; syn6 10 keys +
`"syn6": 1.0`); both stage jsons written: `curriculum_w2lvar_stage4.json`
(−315.000 = the clip of a sub-−315 kine_score — a real score, route works,
pattern simply bad at the seed) and `curriculum_syn6_stage4.json`
(−254.936, a genuine kine_score above the clip). Both sit ABOVE the −320
frozen sentinel and BELOW any future tuned walker, as designed. Stages 2/3/5
share their route with stock lines (2 = stage-1 args minus the afferent/
interleg switches; 3/5 = stock metric keys `bal_*` / `kine_score` verified
against runner.py:1699-1724) and were NOT separately executed — labeled
not run.

### Defaults gate re-run (after all edits)

Command: `"C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
reports_20260925\tmp\gate_ab_w2lvar.py` (asserts `AARL_NET` unset):
```
counts neurons/inputs/synapses = (410, 376, 1186)
GATE A: PASS
```
(GATE B also re-passed: variant 888/382/7430, 23/23 spot-checks, aliases OK.)

### Protected artifacts

- `spinal_run.npz`: size 8,238,797, mtime 2026-09-25 **01:10:02** before AND
  after all runs (baseline captured 14:10, re-checked in
  `tmp\verify_smoke_goal4.py` output: `size=8238797 mtime=1790323802.16`).
- `optuna_walk.db`: mtime 1790323802.33 (01:10 AM, pre-session), still
  **20 studies**, zero `curr_w2lvar_*`/`curr_syn6_*` names leaked.

## Exact relaunch commands (self-contained; no env setup needed)

The scripts pin `AARL_NET` + `AARL_NPZ` themselves as the first act of
`main()` (w2lvar: build_network.py:1026 env route; syn6: env + `G["syn6"]`
routes), so each stage is ONE detached command:
```
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
C:\Users\Ben Bolen\.conda\envs\myo\python.exe _curriculum_w2lvar.py <stage> [n]
C:\Users\Ben Bolen\.conda\envs\myo\python.exe _curriculum_syn6.py   <stage> [n]
```
stages 1..5, defaults 18/18/18/22/22. The two chains may run concurrently
(different dbs, different npz; proven above). Stage order is 1→2→3→4→5;
each stage seeds from `curriculum_{variant}_stage{stage}.json` falling back
to `stage{stage-1}`. Resumable: `load_if_exists=True`; the smoke trials
remain as the first 2 (stage 1) / 1 (stage 4) trials of each study — the
seed trial is trial 0 everywhere.

## Notes / honest caveats

- The w2lvar S2 reflex-motif ranges bracket the params defaults (see
  above) — a deliberate, documented deviation from the [0,0.5] new-gain
  rule because those keys are not new and not default-0.
- Stage-4 seed scores are poor by construction (first contact-free kine
  evaluations for both variants); stage-4/5 tuning is the campaign's next
  step, not this task's.
- syn6's full-runner air rhythm was PARTIAL at the builder gates
  (goal4_build_syn6.md §1); its stage-1 smoke passed here because the
  seed config sustained ≥3 RG-E bursts in the air window (score 61.9),
  but its rhythm fragility is a known tuning surface for stages 2+.
- Stages 2/3/5 were not executed (only stages 1 and 4 were smoked); their
  runner invocations and metric keys are line-verified against stock
  `_curriculum.py` / runner.py, labeled here as not run.
