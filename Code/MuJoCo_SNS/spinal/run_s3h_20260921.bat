@echo off
rem 2026-09-21 s3h: pm_add (ADDITIVE flexor burst in the swing window)
rem on top of the weight-shift/soft-rig space. Seeded from s3g t0.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3h start %date% %time% === >> curriculum_s3h_20260921.log
%PY% _curriculum.py 3 40 >> curriculum_s3h_20260921.log 2>&1 || goto :fail
echo === s3h DONE %date% %time% === >> curriculum_s3h_20260921.log
exit /b 0
:fail
echo === s3h FAILED %date% %time% === >> curriculum_s3h_20260921.log
exit /b 1
