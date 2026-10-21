@echo off
rem 2026-10-20 stage-3 PERTURBATION+METABOLIC study - w2lvar variant
rem (wave 2, Ben's expansion). Fresh study curr_w2lvar_s3p_perturb in
rem optuna_w2lvar.db; script pins AARL_NET/AARL_NPZ itself.
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
echo START %date% %time% >> s3p_w2lvar_20261020.log
D:\Anaconda\envs\myo\python.exe -u _curriculum_s3p.py w2lvar 18 >> s3p_w2lvar_20261020.log 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> s3p_w2lvar_20261020.log
