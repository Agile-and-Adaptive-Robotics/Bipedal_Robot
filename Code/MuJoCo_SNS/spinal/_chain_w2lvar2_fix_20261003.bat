@echo off
rem SLOT A v2: w2lvar fresh post-fix chain - NO fail_stale (a concurrent
rem fail-stale externally failed a live trial on 2026-10-03 23:59; stale
rem trials are now failed manually BEFORE launching this chain).
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe
echo [%DATE% %TIME%] SLOT A2 w2lvar chain start >> CURRICULUM_RUN_20261003.log
"%PY%" _fresh_20261003.py w2lvar 1 21 >> w2lvar_s1_fix_20261003.log 2>&1
echo [%DATE% %TIME%] w2lvar s1 done rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
"%PY%" _fresh_20261003.py w2lvar 2 30 >> w2lvar_s2_fix_20261003.log 2>&1
echo [%DATE% %TIME%] w2lvar s2 done rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
"%PY%" _fresh_20261003.py w2lvar 3 30 >> w2lvar_s3_fix_20261003.log 2>&1
echo [%DATE% %TIME%] SLOT A2 w2lvar chain end rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
