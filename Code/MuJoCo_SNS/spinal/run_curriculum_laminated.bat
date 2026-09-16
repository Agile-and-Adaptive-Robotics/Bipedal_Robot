@echo off
rem Laminated-architecture curriculum re-run (2026-09-16): stage 1 -> 2 -> 3,
rem 100 trials each, sequential seeding via curriculum_stageN.json.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
set LOG=curriculum_lam_20260916.log

echo === STAGE 1 START %date% %time% === >> %LOG%
%PY% _curriculum.py 1 100 >> %LOG% 2>&1 || goto :fail
echo === STAGE 2 START %date% %time% === >> %LOG%
%PY% _curriculum.py 2 100 >> %LOG% 2>&1 || goto :fail
echo === STAGE 3 START %date% %time% === >> %LOG%
%PY% _curriculum.py 3 100 >> %LOG% 2>&1 || goto :fail
echo === ALL STAGES DONE %date% %time% === >> %LOG%
exit /b 0

:fail
echo === CHAIN FAILED %date% %time% === >> %LOG%
exit /b 1
