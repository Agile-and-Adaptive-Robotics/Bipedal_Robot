# Goal-5: s3k winner scored against two references (2026-09-23)

Run: `spinal_run_s3k.npz` (eval walk window from t=2.0 s).
kine_score is kine_ref.compare's objective: LOWER = better match to that reference.

| metric | subject01 (of record) | Falisse Case_40 (predicted) |
|---|---|---|
| n_cycles_r | 11 | 11 |
| n_cycles_l | 0 | 0 |
| duty_r | 0.25 | 0.25 |
| duty_r_ref | 0.61 | 0.57 |
| duty_l | nan | nan |
| duty_l_ref | nan | nan |
| knee_min_r | -48.95 | -48.95 |
| knee_min_r_ref | -69.69 | -59.21 |
| T_r | 0.90 | 0.90 |
| T_r_ref | 1.23 | 1.11 |
| range_hip_r | 8.17 | 8.17 |
| range_hip_r_ref | 43.29 | 51.60 |
| rmse_knee_r | 21.41 | 18.35 |
| mean_hip_r | 19.98 | 19.98 |
| mean_hip_r_ref | 4.34 | 14.84 |

- subject01 (of record): kine_score **-189.45255762990192**

- Falisse Case_40 (predicted): kine_score **-200.45678058907487**

## Reading
- The gap between the two kine_scores measures how reference-specific the s3k tuning is: a walker trained on one reference cycle (subject01) should score visibly better against it than against an independent predicted gait.
- Falisse reference is one periodic cycle (T 1.113 s, duty 0.57, knee_min -59.2 deg) - faster and stiffer than subject01 (1.233 s / 0.61 / -69.7 deg).
- No retraining was run (goal-4 next-study config awaits Ben); this is the measurement half of multi-reference training.
