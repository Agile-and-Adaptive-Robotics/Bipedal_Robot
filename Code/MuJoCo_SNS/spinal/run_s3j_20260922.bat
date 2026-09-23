@echo off
rem 2026-09-22 s3j: FULL Deng-style connectome (full_rules=1 fixed)
rem retune of the operating point. Seeded from s3i t20.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3j start %date% %time% === >> curriculum_s3j_20260922.log
%PY% _curriculum.py 3 40 >> curriculum_s3j_20260922.log 2>&1 || goto :fail
echo === s3j DONE %date% %time% === >> curriculum_s3j_20260922.log
exit /b 0
:fail
echo === s3j FAILED %date% %time% === >> curriculum_s3j_20260922.log
exit /b 1
