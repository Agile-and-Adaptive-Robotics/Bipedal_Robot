@echo off
rem 2026-09-21 s3d: both-leg objective + the t54-diagnosis levers.
rem Searches pf_gain [0.3,3] (log) + contra_swing [0,1.5] on top of the
rem s3c space; seeded from the CORRECTED s3c best (trial 54).
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3d start %date% %time% === >> curriculum_s3d_20260921.log
%PY% _curriculum.py 3 40 >> curriculum_s3d_20260921.log 2>&1 || goto :fail
echo === s3d DONE %date% %time% === >> curriculum_s3d_20260921.log
exit /b 0
:fail
echo === s3d FAILED %date% %time% === >> curriculum_s3d_20260921.log
exit /b 1
