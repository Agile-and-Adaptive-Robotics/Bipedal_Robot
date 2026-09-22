@echo off
rem 2026-09-21 s3g: phase machine v2 - WEIGHT-SHIFT (pm_ws: load-gated
rem swing window + abductor prep asymmetry). Seeded from s3f t32.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3g start %date% %time% === >> curriculum_s3g_20260921.log
%PY% _curriculum.py 3 40 >> curriculum_s3g_20260921.log 2>&1 || goto :fail
echo === s3g DONE %date% %time% === >> curriculum_s3g_20260921.log
exit /b 0
:fail
echo === s3g FAILED %date% %time% === >> curriculum_s3g_20260921.log
exit /b 1
