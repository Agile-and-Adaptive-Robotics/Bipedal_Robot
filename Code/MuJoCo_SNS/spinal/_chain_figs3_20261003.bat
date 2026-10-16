@echo off
rem FIGS round 3: w2lvar fresh stage-1 air winner (the fresh chain
rem reached stage 3 = balance; the walk stages are future work).
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe
echo [%DATE% %TIME%] FIGS3 start >> CURRICULUM_RUN_20261003.log
"%PY%" _replay_best_20261003.py w2lvar optuna_w2lvar.db curr_w2lvar_s1_fix_20261003 spinal_run_w2lvar_s1_best_20261003.npz 1 >> figs_replay3.log 2>&1
"%PY%" _figs_20261003.py spinal_run_w2lvar_s1_best_20261003.npz w2lvar_fix_s1air >> figs_replay3.log 2>&1
echo [%DATE% %TIME%] FIGS3 end rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
