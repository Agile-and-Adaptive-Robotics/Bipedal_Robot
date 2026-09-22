@echo off
rem 2026-09-21 s3e: adds contra_kinh (crossed KINH - opposite heel
rem strike suppresses this side's extensor MNs) on top of s3d space.
rem Seeded from the s3d best (t33).
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3e start %date% %time% === >> curriculum_s3e_20260921.log
%PY% _curriculum.py 3 40 >> curriculum_s3e_20260921.log 2>&1 || goto :fail
echo === s3e DONE %date% %time% === >> curriculum_s3e_20260921.log
exit /b 0
:fail
echo === s3e FAILED %date% %time% === >> curriculum_s3e_20260921.log
exit /b 1
