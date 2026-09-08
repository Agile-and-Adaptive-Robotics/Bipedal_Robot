# ChatGPT handoff — Knee minimizer / Xi-factor work (Bipedal_Robot, Sept 2026)

This file briefs any AI assistant (ChatGPT or otherwise) working on this project while the
primary assistant (ZCode) is unavailable. Read it fully before doing anything. The repo-wide
`AGENTS.md` next to this file has additional standing context — both are plain markdown.

## The person you are helping

Ben Bolen, Mechanical Engineering PhD researcher (PSU, Agile and Adaptive Robotics Lab).
He uses GitHub Desktop (not git CLI), and thinks in terms of the physical test hardware.
Be concrete, show numbers, don't pad.

## The science, briefly

He identifies correction factors (Xi0–Xi3) that map an ideal rigid-body + pneumatic-muscle
(BPA) model onto measured knee-test data, then feeds those into placement-optimization
scripts (`Code\Matlab\Mesh_Optimization\Opt_run*.m`) to design legs whose two 20 mm BPAs
meet or exceed human monoarticular muscle torque. Data: pinned-knee flexor and extensor
tests (5 and 9 length-variants), plus two biomimetic-knee validations.

**Advisor's requirement:** one (Xi1, Xi2) pair consistent across all four configurations
(pinned flexor, pinned extensor, biomimetic flexor, biomimetic extensor). Ben believes a
higher Xi1 is needed for the mesh optimization to reach human-level torque.

**Critical nuance:** Xi1/Xi2 are *effective system-stiffness* parameters — they lump the
bracket, fixtures, AND the cable winch of the test mechanism. "Bending" is a simple
Hooke-law spring in N/m (no EI/L length dependence). Ben shorthand-calls them "bracket
axial/bending stiffness" but they are not literal beam stiffnesses. The bracket reference
points (Pbr, Pbri, Pbr2) and frame conventions are interpretation choices; moving a point
changes the identified stiffness dramatically (one 4 cm move swung flexor Xi1 by 16x).

## Current state (2026-09-07, after today's runs)

Two competing "worlds" fit the pinned-flexor data equally well (mean RMSE ~1.6 vs baseline ~7.6):

| | World A (1trans) | World B (2trans, current code) |
|---|---|---|
| Flexor model | yaw-only bracket frames, pbrB/A = [norm(xy),0,z]+defl | two-rotation (Z then Y) frames, pbrB/A = [norm,0,0]+defl |
| Flexor Xi0 / Xi1 / Xi2 | +7.5 mm / 5.4e4 / 2.7e4 | ~0 / 2.1e6 / 1.1e4 |
| Extensor bracket point Pbr | [-2.65,-54.71,75.06]/1000 | [-3.84,-46.44,62.5]/1000 (rib midpoint, ACTIVE in code) |
| Extensor pick | Xi3=0.294, Xi0=-12.2 mm | Xi3=0.159, Xi0=-6.4 mm |
| Cross-prediction (extensor pool RMSE: published / pure-substitution / refit) | 1.09 / 2.10 / 1.23 | 1.08 / 2.50 / 1.35 |

Landscape finding: raising Xi1 above the flexor value monotonically worsens the extensor
fit (tested to 100x). Xi3 is the helpful dial: at flexor stiffnesses, extensor pool RMSE
improves from 1.41 (Xi3=0) to ~1.01 (Xi3=0.35). Fold-mining of the old CV mats shows the
extensor's best-generalizing folds carry Xi3 = 0.29–0.42 (never trust folds holding out
test 47cm), and Xi0 stable at -10 to -16 mm.

## Active code/model choices (do not silently change)

- `minimizeFlxPin2brk.m` — two-bracket pinned-flexor evaluator. Both brackets use
  two-rotation frames. K (tibia bracket) = [X1, X2, X1]; K2 (origin bracket) = [X2, X1, X2].
  `Pbr2 = [-52.61, 0, 75.06]/1000` (pinned-FLEXOR only — never port it elsewhere).
  `USE_BRACKET2=false` gives the single-bracket ablation.
- `minimizeExtX3.m` — pinned-extensor evaluator, one bracket (two-rotation),
  `Pbr = [-3.84, -46.44, 62.5]/1000` (rib midpoint). Line ~236 passes `Xi0` (not `[]`)
  into Contraction — that fix matters, do not revert it.
- Results naming convention (Ben-mandated): `minimizeFlxPin10_results_<yyyymmdd>_2brkt_<1trans|2trans>.mat`
  and `minimizeExt10mmX3_results_<yyyymmdd>_2trans.mat`. Many old scripts load results by
  exact filename — never rename historical files.

## Rules of engagement

1. Do NOT commit or push anything. Ben reviews via GitHub Desktop.
2. Do NOT launch long optimizer runs (gamultiobj/surrogateopt CVs take 45 min–2 h each)
   without saying so explicitly in your reply.
3. MATLAB is R2025a on this machine. Scripts must run with cwd =
   `Testing_Data\2022_02_Festo`, with `Code\Matlab\Functions`,
   `Code\Matlab\Functions\ModernRobotics`, and `Code\Matlab\Robot_Data` on the path.
4. Never open: point-cloud `.txt` files, `*.mat` in `Previous Optimization Code`, `.asv`
   files (stale autosaves), or any log file beyond tail/grep.
5. If you change a model constant (bracket point, K order, bounds), say so loudly and
   record the OLD value in your report.
6. Plot-quality bar: journal-publication ready if you make plots; otherwise make none.

## Open questions you may be asked to work on

- Which transform convention (1trans vs 2trans) is "right"? Current evidence: both fit the
  pinned flexor equally well (non-identifiable from flexor data alone). Discriminators:
  Lm_p vs Lm_h length-match plots (1trans gave Xi0=+7.5 mm which matched well), physical
  plausibility of Xi1=2.1e6 N/m for an Onyx FDM bracket, and cross-configuration consistency.
- Known data concern: one pinned test used an encoder whose true measured angles may be
  ~+5 degrees offset from reported. Unconfirmed — if fits look angle-shifted on one test
  only, this is the first suspect. Ask Ben which test before assuming.

## When you are done — leave this for ZCode

ZCode (the primary assistant) reads this repo on Ben's main machine and will integrate your
work. Before finishing, either tell Ben to paste this to ZCode, or (preferred) write a file
`CHATGPT_REPORT.md` at the repo root containing:

1. **Files created or modified** — exact paths, and for each: what and why.
2. **Model changes** — any bracket point, transform, K order, bound, or flag that changed,
   with OLD and NEW values.
3. **Runs performed** — script names, configs, result file names (following the naming
   convention), and wall-clock time.
4. **Numbers and conclusions** — the Xi values found, fit metrics vs baseline, and what
   you recommend.
5. **Unfinished business** — anything attempted but not completed, with the error if any.

Do not commit; Ben handles git. ZCode will verify, re-test, and fold the work in.
