@echo off
rem 2026-10-20 stage-3 PERTURBATION+METABOLIC study - MAIN human model.
rem CHAINED: waits for the quiet curr_s3_balance study (same DB, no
rem concurrent sqlite writers) and for its curriculum_stage3.json seed,
rem then runs curr_s3p_perturb (20 trials x 7 evals).
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
echo WAIT-START %date% %time% >> s3p_stock_20261020.log
:wait
findstr /C:"EXITCODE" curr_s3bal_20261020.log >nul 2>&1
if errorlevel 1 (
  timeout /t 60 /nobreak >nul
  goto wait
)
echo CHAIN-FIRED %date% %time% >> s3p_stock_20261020.log
D:\Anaconda\envs\myo\python.exe -u _curriculum_s3p.py stock 20 >> s3p_stock_20261020.log 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> s3p_stock_20261020.log
