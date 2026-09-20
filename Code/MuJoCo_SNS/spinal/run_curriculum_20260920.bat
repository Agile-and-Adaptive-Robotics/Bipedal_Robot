@echo off
rem 2026-09-20 curriculum relaunch: rhythm-gated stage 2 (fresh study
rem curr_s2b_air_aff, seeded from the valid stage-1 winner) then stage 3
rem (curr_s3b_ground, now searching contact_onset). Resumable: each
rem stage load_if_exists + enqueue-seed only on empty study.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === chain start %date% %time% === >> curriculum_rerun_20260920.log
%PY% _curriculum.py 2 25 >> curriculum_rerun_20260920.log 2>&1 || goto :fail
%PY% _curriculum.py 3 30 >> curriculum_rerun_20260920.log 2>&1 || goto :fail
echo === chain DONE %date% %time% === >> curriculum_rerun_20260920.log
exit /b 0
:fail
echo === chain FAILED %date% %time% === >> curriculum_rerun_20260920.log
exit /b 1
