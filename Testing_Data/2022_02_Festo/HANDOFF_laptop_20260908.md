# HANDOFF — laptop zcode session → easteregg2 zcode session (2026-09-08)

**STATUS: ADJUDICATED AND CLOSED (easteregg2 session, 2026-09-08 late, per Ben).** Ben ran
two chats in parallel; this session (laptop, DESKTOP-5Q16KE9) executed zero optimization
runs — it only mined existing results and staged scripts, per his instruction ("mine the
tests for information but don't actually execute new ones"). Everything below that was
STAGED has since been renamed into the Collect/Dig scheme or deleted; see "Fate of the
staged files". The staging caveats no longer apply to what remains, except that neither
kept script has been run end-to-end yet — smoke them first.

## What Ben asked this session (his words, condensed)

1. Mine the minimize results for patterns: is one flexor test obviously garbage and should
   be left out of ALLBPA? Does the data favor numHoldout = 1, numBPA-1, or an even split?
   His hunch: tests 2/3/4 (46cm, 47cm, 40cm-tendon) are good validation tests; tests 1/5
   (48cm, 41cm) are "quirky" and should stay in training to capture more behavior.
2. Try different picks from `minimizeFlxPin10_results_20260907_2brkt_2trans_smoke.mat`.
3. Rerun the 2-bracket 1-transform CV that was killed early on 09-07.

## What the EXISTING data already answers (mining only, no new runs)

- **No test is garbage.** Legacy 0817 leave-2-out CV (`Dig_out\mining_patterns_20260907.txt`):
  every test is predicted well when held out, in every fold — norm RMSE vs baseline:
  48cm 0.25–0.32, 46cm 0.10–0.20, 47cm 0.43–0.56, 40cm-tendon 0.12–0.25, 41cm 0.10–0.19.
  Held-out means: 48cm 0.29 / 46cm 0.15 / 47cm 0.49 / 40t 0.21 / 41cm 0.17.
  47cm is the hardest but rock-stable; nothing blows up anywhere → no test poisons training.
  (Its ~5° encoder lag was later IDENTIFIED as the 47cm test and is corrected in
  `minimizeFlxPin2brk.m`'s kf-build — see AGENTS.md.)
- **Training-set diversity matters mildly**: the worst row is when 46cm+41cm (the two
  easiest) are BOTH held out — everything degrades ~20%. Fits Ben's intuition that the
  quirky/easy extremes carry information.
- **numHoldout: legacy data only covers leave-2 (train=3).** The kept
  `Dig_allbpaNumHoldScan.m` completes the spectrum on the 2brk evaluator at smoke scale:
  leave-1 (train=4), leave-2, Ben's split (train=2), and train=1.
- **Pick structure of the full result (fold-1 front: 265 filtered candidates, all pass
  baseline; the mat stores 3445 filtered across 10 folds):**
  - Xi0 sits at the lower bound 0.0000 for every top pick — and the extensor refit wanted
    −0.0024 (published extensor Xi0 = −0.0122). Ben may want to reconsider the flexor lb=0.
  - **Xi2 is tightly pinned (~1.06–1.11e4 across the top of the front) while Xi1 is a flat
    valley (2e5–2.1e6 with validation distance 0.054–0.059).** "Pick" is effectively a
    choice of Xi1, and the val-distance sort barely discriminates.
  - Consistency story vs the extensor published solution (Xi1=1.508e4, Xi2=1.218e4):
    **Xi2 agrees across configurations (~1.1e4 vs 1.2e4); Xi1 does not (≥2e5 vs 1.5e4).**
    So the advisor's one-(Xi1,Xi2) requirement hinges on Xi1 — which is exactly the
    poorly-constrained parameter. Low-Xi1 picks are the ones worth cross-predicting.
  - crossPredict 09-07 (pick 1, Xi1=2.107e6): pinned-extensor pool mean RMSE
    pub 1.08 | pure-substitution 2.50 | refit(Xi0/Xi3 re-solved, Xi1/Xi2 locked) 1.35.
    Biomimetic flexor 7.45 → 2.91. Biomimetic extensor: baseline 3.03 | pub 1.22 |
    pure 1.21 | refit 0.67. I.e. pick-1's high Xi1 does NOT transfer to the pinned
    extensor, but transfers fine to the biomimetic extensor.

## RESOLVED: the "_smoke" mat question (was an open trap)

The git-tracked `minimizeFlxPin10_results_20260907_2brkt_2trans_smoke.mat` was the FULL
09-07 production CV under a smoke-mode name (no TRANSMODE field — saved by the
pre-416f08f driver; its pick 1 [0, 2.107e6, 1.065e4] matches the 09-07 crossPredict's
source exactly). During the 2026-09-08 Collect/Dig rename it was renamed to the canonical
**`minimizeFlxPin10_results_20260907_2brkt_2trans.mat`**, and easteregg2 verified those
contents (10 folds, ALLBPA [1..5], pick 1 = Xi0 1.65e-6 / Xi1 2.11e6 / Xi2 1.06e4 — the
restored-two-rotation, corrected-Pbr2 chain run). Provenance is therefore CLOSED: there is
no separate smoke mat to worry about, and no smoke-named file should ever be resurrected.

## Fate of the staged files (Ben's adjudication, 2026-09-08)

| staged file | fate |
|---|---|
| `minimizeFlxPin2brk_1trans.m` | **DELETED** — Ben: easteregg2's `transMode` flag on `minimizeFlxPin2brk.m` is the vehicle of record (and 1trans==2trans is proven for these y-symmetric K arrays; the arm-2 campaign already produced the 1trans CV: 4 mats `_20260908_2brkt_1trans_*`) |
| `minimizeFlxPin10mm_2brk_1trans.m` | **DELETED** — same ruling |
| `picksScan_2brkt.m` | **KEPT → `Dig_FlxPin_2brkt_picksScan.m`** (1trans branch routed through the transMode arg; outputs now go to `Dig_out\`) |
| `nightBatch_20260908.m` | **DELETED** — superseded: the 1trans rerun it queued is moot (proof + arm-2 mats), the full 2trans CV already exists (canonical mat), and its picks scan lives on in the kept script |
| `allbpaNumHoldScan_20260908.m` | **KEPT → `Dig_allbpaNumHoldScan.m`** — old (0) provenance check and old (B) angle-lag scan REMOVED (closed/resolved, see above); (A) stored-result per-test table + (C) E-fold numHoldout spectrum remain; Robot_Data addpath added (the 2brk kf-build needs it) |
| `mine_smoke_20260908.m` | **DELETED** — its target name no longer exists; its findings are preserved in this file |
| `HANDOFF_laptop_20260908.md` | this file, updated in place |

## Still runnable (both UNTESTED end-to-end — smoke first)

1. `Dig_allbpaNumHoldScan` — (A) ~1 min, then (C) 14 ga folds at 25×30, ~15–30 min on 6
   workers. Answers Ben's numHoldout question; its E-fold results should set ALLBPA/NUMHOLD
   for any future full runs.
2. `Dig_FlxPin_2brkt_picksScan('minimizeFlxPin10_results_20260907_2brkt_2trans.mat', 12)` —
   the "different picks" answer. Watch whether low-Xi1 picks transfer to the pinned
   extensor better than pick 1 did (~2–3 min fast pass + ~1–2 min per deep pick).

## Still-valid gotchas

- `Dig_crossPredict` saves `Dig_crossPredict_<stamp>.mat` → two calls on the same day
  overwrite; rename/copy aside the first result before rerunning.
- Pre-TRANSMODE results files (the canonical 09-07 mat) have no TRANSMODE field; scripts
  default them to '2trans' — VERIFIED correct for that mat (see RESOLVED above).

## Do NOT

- Don't port flexor `Pbr2` [-52.61, 0, 75.06] to the extensor evaluators (pinned-flexor
  only, per AGENTS.md).
- Don't resurrect any smoke-named mat or the separate-file 1trans evaluator — both are
  retired by Ben's ruling; the canonical mat + the transMode flag replace them.
