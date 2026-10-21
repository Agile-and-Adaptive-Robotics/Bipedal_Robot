@echo off
rem 2026-10-21 resume session: fresh-process re-verification of the
rem three quiet stage-3 balance winners (first-hand numbers for the
rem supervisor gate) + one fresh Ib-ablation reproducibility check.
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=D:\Anaconda\envs\myo\python.exe
set LOG=verify_20261021.log
echo START %date% %time% >> %LOG%
%PY% -u _s3_verify_20261020.py stock >> %LOG% 2>&1
echo RC_stock %ERRORLEVEL% >> %LOG%
%PY% -u _s3_verify_20261020.py w2l_bal >> %LOG% 2>&1
echo RC_w2l_bal %ERRORLEVEL% >> %LOG%
%PY% -u _s3_verify_20261020.py syn6_bal >> %LOG% 2>&1
echo RC_syn6_bal %ERRORLEVEL% >> %LOG%
%PY% -u _s3_ib_ablate_20261020.py stock >> %LOG% 2>&1
echo RC_ibablate_stock %ERRORLEVEL% >> %LOG%
echo EXITCODE 0 %date% %time% >> %LOG%
