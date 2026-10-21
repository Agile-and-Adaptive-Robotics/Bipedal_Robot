@echo off
rem 2026-10-20 stage-3 BALANCE study - MAIN human model (stock stack).
rem The missing curriculum rung: curr_s3_balance was never launched
rem (2026-10-20 audit; DESIGN.md + SUPERVISOR_STATUS note). 40 trials.
rem AARL_NET cleared; AARL_NPZ redirected off the protected default.
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
set AARL_NET=
set AARL_NPZ=spinal_run_curr_s3bal_20261020.npz
echo START %date% %time% >> curr_s3bal_20261020.log
D:\Anaconda\envs\myo\python.exe -u _curriculum.py 3 40 >> curr_s3bal_20261020.log 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> curr_s3bal_20261020.log
