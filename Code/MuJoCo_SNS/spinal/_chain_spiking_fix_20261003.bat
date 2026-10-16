@echo off
rem SLOT C: spiking-mirror curriculum stage 1 (remainder) then a short
rem stage 2. WMI-detached 2026-10-03.
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe
echo [%DATE% %TIME%] SLOT C spiking chain start >> CURRICULUM_RUN_20261003.log
"%PY%" _fail_stale_20261003.py spiking >> CURRICULUM_RUN_20261003.log 2>&1
"%PY%" _curriculum_spiking.py 1 17 >> spk_s1_20261003.log 2>&1
echo [%DATE% %TIME%] spiking s1 done rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
"%PY%" _curriculum_spiking.py 2 12 >> spk_s2_20261003.log 2>&1
echo [%DATE% %TIME%] SLOT C spiking chain end rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
