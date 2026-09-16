@echo off
rem Resume the laminated curriculum on FIXED wiring (2026-09-16):
rem stage 1 tops up to 100 total (86 done when the process died mid-trial,
rem no traceback - one-off native crash; its trials are valid on the fixed
rem code because both fixes live behind conditionals absent at stage-1
rem gains, proven by _fix_check.py), then stages 2 and 3 fresh at 100.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set PY="C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
set LOG=curriculum_lam_20260916.log

for /f %%i in ('%PY% _stage_remaining.py 1') do set REM1=%%i
echo === STAGE 1 RESUME %date% %time% (remaining %REM1%) === >> %LOG%
%PY% _curriculum.py 1 %REM1% >> %LOG% 2>&1 || goto :fail
echo === STAGE 2 START %date% %time% === >> %LOG%
%PY% _curriculum.py 2 100 >> %LOG% 2>&1 || goto :fail
echo === STAGE 3 START %date% %time% === >> %LOG%
%PY% _curriculum.py 3 100 >> %LOG% 2>&1 || goto :fail
echo === ALL STAGES DONE %date% %time% === >> %LOG%
exit /b 0

:fail
echo === CHAIN FAILED %date% %time% === >> %LOG%
exit /b 1
