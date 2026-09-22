@echo off
rem 2026-09-21 s3i: pm_aff afferent disfacilitation (silence the swing
rem leg's load inputs during its own swing window). Seeded s3h t39.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set CONDA_PREFIX=C:\Users\Ben Bolen\.conda\envs\myo
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
echo === s3i start %date% %time% === >> curriculum_s3i_20260921.log
%PY% _seed_s3i.py >> curriculum_s3i_20260921.log 2>&1 || goto :fail
%PY% _curriculum.py 3 40 >> curriculum_s3i_20260921.log 2>&1 || goto :fail
echo === s3i DONE %date% %time% === >> curriculum_s3i_20260921.log
exit /b 0
:fail
echo === s3i FAILED %date% %time% === >> curriculum_s3i_20260921.log
exit /b 1
