@echo off
rem SLOT B3: deepen the BPA studies (cheap trials), then pruned walk
rem seed 2, then extend the s3k mutfix study. No fail_stale inside.
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe
echo [%DATE% %TIME%] SLOT B3 start >> CURRICULUM_RUN_20261003.log
cd w2l_mujoco
"%PY%" _curriculum_bpa.py 40 >> bpa_ns_20261003.log 2>&1
echo [%DATE% %TIME%] BPA ns deepened rc=%ERRORLEVEL% >> ..\CURRICULUM_RUN_20261003.log
"%PY%" _curriculum_bpa_spiking.py 40 >> bpa_spk_20261003.log 2>&1
echo [%DATE% %TIME%] BPA spk deepened rc=%ERRORLEVEL% >> ..\CURRICULUM_RUN_20261003.log
cd ..
"%PY%" _fresh_20261003.py pruned walk 2 40 >> pruned_walk2_fix_20261003.log 2>&1
echo [%DATE% %TIME%] pruned walk2 done rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
"%PY%" _resume_20261003.py s3k 15 >> resume_s3k_round2.log 2>&1
echo [%DATE% %TIME%] SLOT B3 end rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
