@echo off
rem SLOT B v2: pruned air (remainder) + walk + BPA curricula. No
rem fail_stale inside chains (2026-10-03 lesson).
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe
echo [%DATE% %TIME%] SLOT B2 pruned+BPA chain start >> CURRICULUM_RUN_20261003.log
"%PY%" _fresh_20261003.py pruned air 23 >> pruned_air_fix_20261003.log 2>&1
echo [%DATE% %TIME%] pruned air done rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
"%PY%" _fresh_20261003.py pruned walk 1 40 >> pruned_walk1_fix_20261003.log 2>&1
echo [%DATE% %TIME%] pruned walk1 done rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
cd w2l_mujoco
"%PY%" _curriculum_bpa.py 25 >> bpa_ns_20261003.log 2>&1
echo [%DATE% %TIME%] BPA ns done rc=%ERRORLEVEL% >> ..\CURRICULUM_RUN_20261003.log
"%PY%" _curriculum_bpa_spiking.py 25 >> bpa_spk_20261003.log 2>&1
echo [%DATE% %TIME%] SLOT B2 pruned+BPA chain end rc=%ERRORLEVEL% >> ..\CURRICULUM_RUN_20261003.log
