@echo off
rem 2026-10-21 stage-3Q FIXED-OBJECTIVE study - MAIN human model (Ben's
rem order: rig-share term + metabolic term added to _curriculum_s3p.py;
rem see that docstring). REAL budget: 60 trials x 7 evals. New study
rem name curr_s3q_perturb - the recorded curr_s3p_* study is untouched.
rem Seeds: trial 0 = quiet incumbent defaults, trial 1 = the 2026-10-20
rem s3p winner params (curriculum_s3p_stock.json).
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
echo LAUNCH-START %date% %time% >> s3q_stock_20261021.log
D:\Anaconda\envs\myo\python.exe -u _curriculum_s3p.py stock 60 >> s3q_stock_20261021.log 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> s3q_stock_20261021.log
