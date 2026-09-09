# HANDOFF — laptop zcode session → easteregg2 zcode session (2026-09-08)

**Read this first if you are the other zcode chat.** Ben ran two chats in parallel and they
got mixed up. This session ran on the laptop (DESKTOP-5Q16KE9) and **executed zero
optimization runs** — it only mined existing results and staged scripts, per Ben's
instruction ("mine the tests for information but don't actually execute new ones").
Everything below is STAGED and UNTESTED. **You own this folder's compute — rename, rewrite,
or delete any of it freely.** Ben commits this via GitHub Desktop.

## What Ben asked this session (his words, condensed)

1. Mine the minimize results for patterns: is one flexor test obviously garbage and should
   be left out of ALLBPA? Does the data favor numHoldout = 1, numBPA-1, or an even split?
   His hunch: tests 2/3/4 (46cm, 47cm, 40cm-tendon) are good validation tests; tests 1/5
   (48cm, 41cm) are "quirky" and should stay in training to capture more behavior.
2. Try different picks from `minimizeFlxPin10_results_20260907_2brkt_2trans_smoke.mat`.
3. Rerun the 2-bracket 1-transform CV that was killed early on 09-07.

## What the EXISTING data already answers (mining only, no new runs)

- **No test is garbage.** Legacy 0817 leave-2-out CV (`mining_patterns_20260907.txt`):
  every test is predicted well when held out, in every fold — norm RMSE vs baseline:
  48cm 0.25–0.32, 46cm 0.10–0.20, 47cm 0.43–0.56, 40cm-tendon 0.12–0.25, 41cm 0.10–0.19.
  Held-out means: 48cm 0.29 / 46cm 0.15 / 47cm 0.49 / 40t 0.21 / 41cm 0.17.
  47cm is the hardest but rock-stable; nothing blows up anywhere → no test poisons training.
- **Training-set diversity matters mildly**: the worst row is when 46cm+41cm (the two
  easiest) are BOTH held out — everything degrades ~20%. Fits Ben's intuition that the
  quirky/easy extremes carry information.
- **numHoldout: legacy data only covers leave-2 (train=3).** The staged
  `allbpaNumHoldScan_20260908.m` completes the spectrum on the 2brk evaluator at smoke
  scale: leave-1 (train=4), leave-2, Ben's split (train=2), and train=1, plus the encoder
  angle-lag test below.
- **Pick structure of the full result (265 filtered candidates, ALL pass baseline):**
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

## ⚠ File gotcha: the "_smoke" mat is actually the FULL production run

`minimizeFlxPin10_results_20260907_2brkt_2trans_smoke.mat` contains ALLBPA=[1..5],
NUMHOLD=2, POP=150, MAXGEN=600, SOLVER=gamultiobj, 265/265 candidates passing the baseline
filter, and NO `TRANSMODE` field → saved by the pre-416f08f driver (old naming; the 09-07
crossPredict sourced it on D: as `minimizeFlxPin10_2brk_results_20260907.mat`, pick 1
[0, 2.107e6, 1.065e4] matches this file's pick 1 exactly).

**Provenance unverified:** whether its stored scores came from the restored TWO-rotation
evaluator or the superseded single-transform one. Step 0 of
`allbpaNumHoldScan_20260908.m` resolves this in ~1 min (re-evaluates one stored candidate
with both evaluators and compares). Also note `crossPredictFlx.m` auto-pick skips any file
containing "_smoke" — pass this file explicitly or rename it.

## Staged files (all UNTESTED — smoke them before trusting)

| file | what |
|---|---|
| `minimizeFlxPin2brk_1trans.m` | evaluator: exact recreation of the superseded single-transform (yaw-only frames, (d) pbrBnew/pbrAnew) version from **git commit 3374837's successor 3374847**; function renamed with `_1trans` suffix |
| `minimizeFlxPin10mm_2brk_1trans.m` | driver copy routed to that evaluator, `TRANSMODE='1trans'` → saves `..._2brkt_1trans.mat` (no clobber) |
| `picksScan_2brkt.m` | **fast pass scores the WHOLE filtered front** (one 5-test evaluator call per candidate → CSV + plateau stats), then deep cross-prediction (bio flexor, extensor pure + coarse refit, bio extensor) on picks 1..N. Evaluator variant from the file's TRANSMODE or 3rd arg |
| `nightBatch_20260908.m` | queued sequence: picks scan → full 1trans CV (the killed rerun, surrogateopt) → crossPredict+picks → full 2trans CV → crossPredict+picks. Self-locating paths (safe on D: or C:). ~4–5 h on 6 workers, faster on your 10 |
| `allbpaNumHoldScan_20260908.m` | step 0 provenance check + per-test held-out table from the stored full result + angle-lag encoder scan (the ±5° discriminator, flexor 5 tests across 5 solutions incl. baseline models) + new smoke folds E1/E3/E4/E4rev/E5 |
| `mine_smoke_20260908.m` | read-only dump of the full mat + crossPredict summary (already run once on the laptop; findings are in this file) |

## Suggested order when Ben says go

1. `allbpaNumHoldScan_20260908` (~25 min on 6 workers, less on 10): resolves provenance,
   answers Ben's test/numHoldout questions, identifies the encoder-shift suspect if there
   is one. Its E-fold results should set ALLBPA/NUMHOLD for the full runs.
2. `picksScan_2brkt` on the 265-candidate result (nPicks 12): the "different picks" answer.
   Watch whether low-Xi1 picks transfer to the pinned extensor better than pick 1 did.
3. `nightBatch_20260908` (or your own runner): the full 1trans rerun + a full 2trans
   production run for comparison. Solver note: abSolver A/B on 09-07 gave gamultiobj
   333 s/dist 0.309 vs surrogateopt 231 s/dist 0.350 → the rule picks surrogateopt; full
   gamultiobj at POP 150 × MAXGEN 600 × 10 folds is a DAY-length run on the laptop.

## Known bugs/gotchas found (not fixed — existing files left untouched for you)

- `runBatch.m` step 2 globs `minimizeFlxPin10_2brk_results_*.mat` but the 416f08f driver
  saves `minimizeFlxPin10_results_*_2brkt_*.mat` → step 2 can't find step 1's output.
- `runBatch.m` and `smoke_2brk.m` hardcode `D:/GitHub/...` paths (fine on easteregg2,
  wrong on the laptop — no D: here).
- `crossPredictFlx` saves `crossPredict_<stamp>.mat` → two calls on the same day
  overwrite; nightBatch renames after each call.
- Pre-TRANSMODE results files (like the misnamed one above) have no TRANSMODE field;
  picksScan defaults them to '2trans' — verify with the provenance check first.

## Do NOT

- Don't port flexor `Pbr2` [-52.61, 0, 75.06] to the extensor evaluators (pinned-flexor
  only, per AGENTS.md).
- Don't treat the 1trans evaluator as the method of record — two-rotation is restored and
  supersedes it; the 1trans rerun is comparison data for the run Ben saw killed.
- Don't trust the staged files sight-unseen: none were executed. `minimizeFlxPin2brk_1trans.m`
  body is byte-for-byte from git 3374847 except the function name and a 2-line header.
