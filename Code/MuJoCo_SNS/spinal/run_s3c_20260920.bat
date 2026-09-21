@echo off
rem 2026-09-20 s3c retune: both-leg objective (kine_ref v2), searches
rem pelvis_ty + f1_anklepf_inh, seeded from the s3b ground winner.
rem Study curr_s3c_ground is fresh; resumable via _curriculum.py 3 N.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3c chain start %date% %time% === >> curriculum_s3c_20260920.log
%PY% _curriculum.py 3 40 >> curriculum_s3c_20260920.log 2>&1 || goto :fail
echo === s3c DONE %date% %time% === >> curriculum_s3c_20260920.log
exit /b 0
:fail
echo === s3c FAILED %date% %time% === >> curriculum_s3c_20260920.log
exit /b 1
