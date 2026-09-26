# Easteregg2 continuation handoff — w2lvar/syn6 curriculum generations

**Written:** 2026-09-26 (EB475WS4 session, after the 09-25 campaign + Ben's figure review).
**Purpose:** Ben's takeaway from the 09-25 campaign is that the variant curricula were
*shallow* (18–23 trials/stage vs s3k's 80-trial study). The heavy-box (easteregg2,
10 cores / 128 GB) should run more generations on the same studies. This file is the
complete recipe. Easteregg2 sessions: read this plus `EXEC_SUMMARY_20260925.md` first.

## TL;DR

```bat
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
:: sanity gates FIRST (see below) — then, per stage (n = ADDITIONAL trials):
C:\Users\Ben Bolen\.conda\envs\myo\python.exe ...   :: NO — easteregg2 python is:
D:\Anaconda\envs\myo\python.exe _curriculum_w2lvar.py <stage> <n>
D:\Anaconda\envs\myo\python.exe _curriculum_syn6.py <stage> <n>
```

- The scripts **self-select the variant**: `_curriculum_w2lvar.py` / `_curriculum_syn6.py`
  pin `AARL_NET` and their own `AARL_NPZ` (`spinal_run_w2lvar.npz` / `spinal_run_syn6.npz`)
  at `main()` — no env setup needed. Run from the spinal cwd (the optuna db paths are relative).
- Studies are **resumable** (`load_if_exists=True`); the dbs are **committed to git**
  (`optuna_w2lvar.db`, `optuna_syn6.db`), so after a pull the studies continue with their
  full TPE trial history. Do NOT create new study names.
- Suggested budgets per variant: s1 +30, s2 +30, s3 +30, s4 +60, s5 +60.
  Kine evals ran ≈1–2 min each on EB475WS4; the two chains ran concurrently there
  (separate dbs + npz, no interference) and can do the same on easteregg2.
  Never run two stages of the SAME chain concurrently (shared db + npz).
- Each stage seeds from `curriculum_<variant>_stage<N-1>.json` (committed), so you can
  go straight to deepening s4/s5 without touching s1–s3.

## Prereqs

1. Repo at commit **5f1890d7** ("Figures", 2026-09-25 22:43) or later — that commit contains
   the whole campaign: both builders (`build_network_w2lvar.py`, `build_network_syn6.py`),
   the selector in `build_network.py`, `params.py` (`syn6`, `syn6_brainstem` keys),
   `_curriculum_w2lvar.py` / `_curriculum_syn6.py`, both dbs, all winner jsons, reports.
2. Env: `D:\Anaconda\envs\myo\python.exe` (optuna 5.0.0, schema matches — db is schema 12).
3. Post-commit EB475WS4 additions (optional for tuning; in the next push):
   `reports_20260925\tmp\probe_ref_flat.py`, `reports_20260925\tmp\prove_ref_left_bug.py`,
   `reports_20260925\figs\ref_left_bug_proof.png`, and the connectome-editor template
   additions for w2lvar/syn6.

## Sanity gates — run these BEFORE launching on easteregg2

Different machine ⇒ re-verify counts before trusting anything:

```bat
D:\Anaconda\envs\myo\python.exe reports_20260925\tmp\gate_ab_w2lvar.py
:: Gate A (AARL_NET unset)  -> (410, 376, 1186) PASS   [the defaults gate — mandatory]
:: Gate B (w2lvar)          -> (888, 382, 7430)
D:\Anaconda\envs\myo\python.exe reports_20260925\tmp\gate_syn6.py b
:: -> (794, 382, 3033)
```

Then a 1-trial smoke per chain: `... _curriculum_w2lvar.py 1 1` and `... _curriculum_syn6.py 1 1`
(exit 0, trial recorded in the db, npz written, objective log line "npz columns verified").

## OPEN DECISION before burning s4/s5 generations — the kine_ref left-reference bug

Ben caught this reviewing the 09-25 overlay (2026-09-26): `kine_ref.load_reference()` cuts
the LEFT reference cycle at the first left GRF onset (t = 0.005 s), but the IK file only
covers t ≥ 0.50 s — `np.interp` edge-fill therefore **freezes `ref['l']` for the first
39.6% of the left cycle** (hip −14.4°, knee −2.1°, ankle +7.4°). Affects every
`ref['l']`-based score since the per-side left reference landed (~`a4738f92`, Sep 21 —
the curr_s3* → s3k → variants lineage). The v1 (Sep 12) reference was right-leg-only and clean.
Proof: `figs\ref_left_bug_proof.png`; probe: `tmp\probe_ref_flat.py`.

- Option (a) **fix first**: cut the left cycle at the first onset pair fully inside IK
  support (1.255→2.468 s) + assert full coverage in `load_reference`. Changes the objective —
  winners from 09-25 stop being strictly comparable, but further generations stop
  optimizing against a biased target.
- Option (b) **keep as-is**: strict comparability with the recorded winners; known bias
  (left early-stance under-activity partially rewarded — one contributing factor in syn6's
  planted-left outcome).

**Ben owns this call.** EB475WS4 did not modify `kine_ref.py`. The fix is ~5 lines and
should land as a reviewed commit before option (a) generations start.

## Harvest discipline (unchanged from 09-25)

- Air stages: winners must satisfy the rhythm gate (rises ≥ 3; a static pose must never
  outscore a walker). Every stage on 09-25 has an `*_exploit_check.log` in
  `reports_20260925\logs\` — reproduce that check on easteregg2 harvests.
- Stage winners overwrite `curriculum_<variant>_stage<N>.json` in the spinal cwd.
  Preserve the 09-25 winners first (copy to `reports_20260925\easteregg2\pre_*.json`
  or let git show the diff) so improvements are provable.
- Known caveats to carry: syn6 rhythm fragility (its s4/s5 winners sit AT the rises=3
  gate minimum) and left-planned leg; w2lvar left leg genuinely slower
  (duty 0.81 vs 0.56, T 1.40 vs 1.17 s).

## Do NOT touch

`optuna_walk.db` (the legacy 20 studies), `spinal_run.npz`, `reports_20260923\`,
`reports_20260924\`, and the default (AARL_NET-unset) build behavior — the defaults gate
above is the tripwire.

## Bring-back

Copy updated dbs + stage jsons + logs into `spinal\reports_20260925\easteregg2\` and
append a dated section to `goal4_variant_curricula.md` (or a new
`goal4_easteregg2_generations.md`): trials added, new winners, score deltas vs the 09-25
record, exploit checks. Ben commits via GitHub Desktop as usual.

## Connectome editor note

`w2lvar` and `syn6` now have editor templates (added 2026-09-26; regenerate via
`make_editor_templates.py` if the builders' topology changes). The editor is a **snapshot**
of the architecture, not a live view — live traces are `neuro_scope.py` / `runner --scope`.
