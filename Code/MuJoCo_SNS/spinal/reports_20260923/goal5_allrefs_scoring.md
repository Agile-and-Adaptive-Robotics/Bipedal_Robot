# s3k production walker vs the full thumb-drive library (2026-09-24)

Run: `spinal_run_s3k.npz`, walk window from t=2.0 s. 44 scored references (lower kine_score = better).

| ref | kine_score | sim duty/T/knee | ref duty/T/knee |
|---|---|---|---|
| subject01_of_record | -189.5 | 0.25 / 0.90s / -49deg | |
| Case_40_motion | -200.5 | 0.25 / 0.90s / -49deg | |
| ong_selfsel_Init050 | -205.5 | 0.25 / 0.90s / -49deg | |
| ong_selfsel_Init125 | -216.1 | 0.25 / 0.90s / -49deg | |
| ong_selfsel_Init200 | -215.5 | 0.25 / 0.90s / -49deg | |
| ong_speed_050 | -203.3 | 0.25 / 0.90s / -49deg | |
| ong_speed_075 | -208.2 | 0.25 / 0.90s / -49deg | |
| ong_speed_100 | -207.4 | 0.25 / 0.90s / -49deg | |
| ong_speed_125 | -211.9 | 0.25 / 0.90s / -49deg | |
| ong_speed_150 | -215.7 | 0.25 / 0.90s / -49deg | |
| ong_speed_175 | -218.4 | 0.25 / 0.90s / -49deg | |
| ong_speed_200 | -224.4 | 0.25 / 0.90s / -49deg | |
| subject01_Run_20002 | -197.9 | 0.25 / 0.90s / -49deg | |
| subject01_Run_30002 | -231.7 | 0.25 / 0.90s / -49deg | |
| subject01_Run_40002 | -260.9 | 0.25 / 0.90s / -49deg | |
| subject01_Run_50002 | -288.2 | 0.25 / 0.90s / -49deg | |
| subject02_Run_20002 | -205.7 | 0.25 / 0.90s / -49deg | |
| subject02_Run_30002 | -236.6 | 0.25 / 0.90s / -49deg | |
| subject02_Run_40005 | -279.1 | 0.25 / 0.90s / -49deg | |
| subject02_Run_50002 | -307.1 | 0.25 / 0.90s / -49deg | |
| subject04_Run_20002 | -192.2 | 0.25 / 0.90s / -49deg | |
| subject04_Run_30002 | -237.2 | 0.25 / 0.90s / -49deg | |
| subject04_Run_40002 | -277.9 | 0.25 / 0.90s / -49deg | |
| subject04_Run_50004 | -311.8 | 0.25 / 0.90s / -49deg | |
| subject08_Run_20002 | -188.2 | 0.25 / 0.90s / -49deg | |
| subject08_Run_30001 | -222.5 | 0.25 / 0.90s / -49deg | |
| subject08_Run_40002 | -265.1 | 0.25 / 0.90s / -49deg | |
| subject08_Run_50003 | -293.6 | 0.25 / 0.90s / -49deg | |
| subject10_Run_20002 | -203.1 | 0.25 / 0.90s / -49deg | |
| subject10_Run_30002 | -243.7 | 0.25 / 0.90s / -49deg | |
| subject10_Run_40002 | -303.8 | 0.25 / 0.90s / -49deg | |
| subject10_Run_50002 | -332.2 | 0.25 / 0.90s / -49deg | |
| subject11_Run_20002 | -190.0 | 0.25 / 0.90s / -49deg | |
| subject11_Run_30002 | -224.1 | 0.25 / 0.90s / -49deg | |
| subject11_Run_40002 | -246.8 | 0.25 / 0.90s / -49deg | |
| subject11_Run_50002 | -276.4 | 0.25 / 0.90s / -49deg | |
| subject17_Run_20002 | -196.4 | 0.25 / 0.90s / -49deg | |
| subject17_Run_30001 | -251.9 | 0.25 / 0.90s / -49deg | |
| subject17_Run_40001 | -316.4 | 0.25 / 0.90s / -49deg | |
| subject17_Run_50001 | -326.9 | 0.25 / 0.90s / -49deg | |
| subject20_Run_20001 | -208.0 | 0.25 / 0.90s / -49deg | |
| subject20_Run_30001 | -256.4 | 0.25 / 0.90s / -49deg | |
| subject20_Run_40001 | -297.1 | 0.25 / 0.90s / -49deg | |
| subject20_Run_50002 | -320.3 | 0.25 / 0.90s / -49deg | |

All 32 Arnold trials are RUNNING (duty 0.31-0.46, knee_min -80..-140 deg across 4 speeds x 8 subjects); Falisse Case_40 is predicted WALKING. The s3k walker is a WALKER (duty ~0.69, knee ~-53): it should match walking refs far better than running refs - the spread quantifies gait-specificity of the current tuning.

Training-side use: swap kine_ref.REF_CACHE per trial (gait_lib_pilot2.py pattern) to build a multi-reference / multi-gait curriculum.
