@echo off
rem SLOT C v2: spiking s1 (remainder to ~25) then a short s2. No
rem fail_stale inside chains.
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe
echo [%DATE% %TIME%] SLOT C2 spiking chain start >> CURRICULUM_RUN_20261003.log
"%PY%" _curriculum_spiking.py 1 20 >> spk_s1_20261003.log 2>&1
echo [%DATE% %TIME%] spiking s1 done rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
"%PY%" _curriculum_spiking.py 2 12 >> spk_s2_20261003.log 2>&1
echo [%DATE% %TIME%] SLOT C2 spiking chain end rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
