@echo off
rem 2026-09-21 s3f: per-side contact-reset PHASE MACHINE (heel-strike
rem reset, antiphase coupling, swing-window MN gating) + the s3e space.
rem Seeded from the s3e best (stage3.json) + pm_gain 0.5, pm_T 1.23.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3f start %date% %time% === >> curriculum_s3f_20260921.log
%PY% _curriculum.py 3 40 >> curriculum_s3f_20260921.log 2>&1 || goto :fail
echo === s3f DONE %date% %time% === >> curriculum_s3f_20260921.log
exit /b 0
:fail
echo === s3f FAILED %date% %time% === >> curriculum_s3f_20260921.log
exit /b 1
