@echo off
rem 2026-10-03 tuning-seed: stock s3k-lineage stage 4 under a fresh
rem dated study name (post mutual-inhibition-fix). Same idiom as
rem run_curriculum_20260920.bat. Resumable: load_if_exists + seed only
rem on an empty study.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
set LOG="D:\Github\Bipedal_Robot\tmp\tuning_seed_20261003\s4_mutfix_20261003.log"
echo === s4-mutfix start %date% %time% === >> %LOG%
%PY% "D:\Github\Bipedal_Robot\tmp\tuning_seed_20261003\run_s4_mutfix_20261003.py" >> %LOG% 2>&1 || goto :fail
echo === s4-mutfix DONE %date% %time% === >> %LOG%
exit /b 0
:fail
echo === s4-mutfix FAILED %date% %time% === >> %LOG%
exit /b 1
