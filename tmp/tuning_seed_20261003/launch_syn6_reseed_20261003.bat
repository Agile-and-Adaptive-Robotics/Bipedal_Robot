@echo off
rem 2026-10-03 tuning-seed: syn6 stage 5 fresh-seed rescan under a
rem fresh dated study name. Same idiom as run_curriculum_20260920.bat.
rem Resumable: load_if_exists + seed only on an empty study.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
set LOG="D:\Github\Bipedal_Robot\tmp\tuning_seed_20261003\syn6_reseed_20261003.log"
echo === syn6-reseed start %date% %time% === >> %LOG%
%PY% "D:\Github\Bipedal_Robot\tmp\tuning_seed_20261003\run_syn6_reseed_20261003.py" >> %LOG% 2>&1 || goto :fail
echo === syn6-reseed DONE %date% %time% === >> %LOG%
exit /b 0
:fail
echo === syn6-reseed FAILED %date% %time% === >> %LOG%
exit /b 1
