# BW Tier 2 verdict — conductance-class solve (2026-10-24, easteregg2)

Study: `curr_w2l_condsolve`, 80 trials, `spinal\reports_20261024_condsolve\`
(db `optuna_w2l_condsolve.db`, 15 class dims + pinned drive te=3.0/tf=4.0,
metric + survival gate VERBATIM via `_w2l_wean_eval`; per-term diagnostics
in `condsolve_trials.csv` + optuna user_attrs).

## Headline results

1. **The seed won; the search added nothing.** Best trial = seed 1 =
   **Ben's literal corrected-aproj defaults** (ch_hc 0.08, ch_ank 1.0,
   ch_contra 0.08, ctn_amp 4.0, cref 50, ibnA 1.0, comm 1.0; every new
   class scale s_* = 1.0 = its documented W2L_GAINS anchor). obj20 = 90.26
   (matches the 2026-10-20 record "Ben's UNTUNED defaults survive the walk
   gate, obj 90.0"), fresh-process 30 s verification SURVIVES (33 cycles L,
   phase healthy). No searched configuration beat the documented values at
   pinned drive — 80 TPE trials over the 15 conductance classes.
2. **Strict basin gate: 10/12 — FAIL vs the 11/12 bar.** Both failures
   (trials 5 and 9) are phase excursions (phase_med 0.27 / 0.13 with both
   legs still cycling 25-33 cycles) — the same documented fragile piece as
   the contact2 campaign (best there: 9/12). 10/12 improves on 9/12 but
   does not clear the bar.

## Verdict (honest, per the batch's own rule)

- **NOT accepted as a winner.** The plateau below 11/12 is reported as
  **DESIGN.md open-issue-#4 evidence**: the phase-lock fragility of the
  corrected W2L network is structural at the conductance-class level — a
  15-class magnitude search at pinned drive neither improves the objective
  over Ben's defaults nor hardens the phase robustness past 10/12. A
  restructure decision (Ben's) is the documented next step.
- The tuned contact2 winner's obj 91.4 (drive te 4.35/tf 6.79) remains the
  family's best score — with drive FREE. At pinned drive the defaults are
  the optimum this search found.
- Anti-lottery: winner = documented defaults (not searched weights), so
  remove-and-retest for lottery weights is not applicable; the ±1% basin
  IS the local gain-sensitivity probe (10/12 pass); per-term diagnostics
  recorded per trial; determinism shown by fresh-process replay matching.
- No wean-ladder claims. No harness-free walking claim.

## Bit-identity and contract checks

- New s_* knobs default 1.0 = anchors; **default build proven
  bit-identical** (G-signature + full-net signature comparison, 96/204
  census, this folder's run record).
- JSON RULE: winner record `condsolve_best.json` carries every knob incl.
  all s_*; loader branch = `test_w2l_ground.py` KV → gains plumbing
  (CLI `--k=v` sets KV; documented in the driver docstring).
- s3k default census + anchor gate re-run after the Tier-1 generator edit:
  **PASS (410/376/1186 + −161.56754173676563 bit-exact)**; the Tier-2
  edit touches only the w2l build path (`test_w2l_ground.py` bilateralrg2
  branch), which the s3k gate does not exercise.

## Files

`spinal\reports_20261024_condsolve\`: `_w2l_condsolve_campaign.py`,
`_verify_finalist.py`, `_basin_condsolve.py`, `deviation_table.md`,
`condsolve_best.json`, `condsolve_trials.csv`, `verify_trials.csv`,
`basin_condsolve.json`, `condsolve_run.log`, `basin_run.log`.
Ben's originals untouched: the .aproj untouched (flagged F1/V3 items remain
for his GUI pass); edits limited to `make_editor_templates.py` (bilateralrg2
carry-through, Tier 1) and `test_w2l_ground.py` (s_* knobs, default-off) —
both documented above for his commit.
