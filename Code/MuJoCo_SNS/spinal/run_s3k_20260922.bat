@echo off
rem 2026-09-22 s3k: retune the SEVERED-L/R configuration (no_cross
rem fixed on + full_rules on) - the only config that walked bilaterally.
rem Seeded from the s3i winner (t20) via stage3.json.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3k start %date% %time% === >> curriculum_s3k_20260922.log
%PY% _update_stage3.py curr_s3i_ground 20 >> curriculum_s3k_20260922.log 2>&1
%PY% _curriculum.py 3 40 >> curriculum_s3k_20260922.log 2>&1 || goto :fail
echo === s3k DONE %date% %time% === >> curriculum_s3k_20260922.log
exit /b 0
:fail
echo === s3k FAILED %date% %time% === >> curriculum_s3k_20260922.log
exit /b 1
