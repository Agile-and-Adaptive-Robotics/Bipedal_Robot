@echo off
rem FIGS round 2: the two replays that hit the stdout-GC bug (now fixed
rem in _replay_best_20261003.py).
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe
echo [%DATE% %TIME%] FIGS2 start >> CURRICULUM_RUN_20261003.log
"%PY%" _replay_best_20261003.py pruned optuna_pruned_pruned_walk_1_fix_20261003.db pruned_walk_1_fix_20261003 spinal_run_pruned_best_20261003.npz >> figs_replay2.log 2>&1
"%PY%" _figs_20261003.py spinal_run_pruned_best_20261003.npz s3kpruned >> figs_replay2.log 2>&1
"%PY%" _replay_best_20261003.py spiking optuna_spiking.db curr_spk_s1_air_deaff_20261003 spinal_run_spk_s1_best_20261003.npz 1 >> figs_replay2.log 2>&1
"%PY%" _figs_20261003.py spinal_run_spk_s1_best_20261003.npz spiking_mirror_s1 >> figs_replay2.log 2>&1
echo [%DATE% %TIME%] FIGS2 end rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
