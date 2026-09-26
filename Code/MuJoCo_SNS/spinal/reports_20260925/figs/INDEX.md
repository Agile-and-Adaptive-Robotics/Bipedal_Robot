# figs/ — SCONE figures for the 2026-09-25 walker campaign

Regenerate any figure with: `"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" tools\make_scone_figs.py` (parser: `tools\scone_sto.py`; reads only `logs\motion_reeval_*.sto`, no sims).

- `scone_gait_overfit.png` — THE overfit figure: Tutorial 4a gait pelvis forward distance + height vs time, pretrained vs 4-s-window CMA-ES best vs 20-s full-horizon best; the D4 par falls at 4.68 s / 5.05 m (20-s fitness 83.3 = fell) while pretrained walks 21.48 m and D20 best 21.08 m, both the full 20 s / 36 steps.
- `scone_balance_gains.png` — Tutorial 3a balance pelvis height vs time: default falls at 1.04 s, 40-gen CMA-ES best at 2.81 s, 300-gen best stands the full 30 s at 0.93 m (fitness 101.8 → 94.5 → 2.19).
- `scone_tutorial_scores.png` — horizontal bar chart of the 15 first-ever tutorial evaluations at default params (Hyfydy, sconecmd -l 2): orange = maximize-style (jump 1/2a/2b/2c + ScriptMeasure 6a/6d, higher = better), blue = minimize-style (gait/balance, lower = better); scores not comparable across styles; values transcribed from goal1_scone_tutorials_and_optimization.md §1.
