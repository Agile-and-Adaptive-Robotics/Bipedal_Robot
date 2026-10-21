@echo off
rem 2026-10-21 resume session: fresh-process verification of the three
rem stage-3p perturbation winners (battery re-run, score vs DB value).
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=D:\Anaconda\envs\myo\python.exe
set LOG=s3p_verify_20261021.log
echo START %date% %time% >> %LOG%
%PY% -u _s3p_verify_20261021.py stock >> %LOG% 2>&1
echo RC_stock %ERRORLEVEL% >> %LOG%
%PY% -u _s3p_verify_20261021.py syn6 >> %LOG% 2>&1
echo RC_syn6 %ERRORLEVEL% >> %LOG%
%PY% -u _s3p_verify_20261021.py w2lvar >> %LOG% 2>&1
echo RC_w2lvar %ERRORLEVEL% >> %LOG%
echo EXITCODE 0 %date% %time% >> %LOG%
