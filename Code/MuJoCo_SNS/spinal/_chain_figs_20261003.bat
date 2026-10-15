@echo off
rem FINAL FIGURES chain: replay each finished model's winner (npz) and
rem build the slideshow figures + gifs. Read-only on the dbs.
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe
echo [%DATE% %TIME%] FIGS chain start >> CURRICULUM_RUN_20261003.log
"%PY%" _replay_best_20261003.py pruned optuna_pruned_pruned_walk_1_fix_20261003.db pruned_walk_1_fix_20261003 spinal_run_pruned_best_20261003.npz >> figs_replay.log 2>&1
"%PY%" _figs_20261003.py spinal_run_pruned_best_20261003.npz s3kpruned >> figs_replay.log 2>&1
cd w2l_mujoco
"%PY%" _bpa_runbest_20261003.py ns ..\curriculum_bpa_ns.json ..\bpa_ns_final_20261003.npz >> ..\figs_replay.log 2>&1
"%PY%" _bpa_runbest_20261003.py spk ..\curriculum_bpa_spk.json ..\bpa_spk_final_20261003.npz >> ..\figs_replay.log 2>&1
cd ..
"%PY%" _figs_20261003.py bpa_ns_final_20261003.npz bpa_nonspiking --format bpa >> figs_replay.log 2>&1
"%PY%" _figs_20261003.py bpa_spk_final_20261003.npz bpa_spiking --format bpa >> figs_replay.log 2>&1
"%PY%" _replay_best_20261003.py spiking optuna_spiking.db curr_spk_s1_air_deaff_20261003 spinal_run_spk_s1_best_20261003.npz 1 >> figs_replay.log 2>&1
"%PY%" _figs_20261003.py spinal_run_spk_s1_best_20261003.npz spiking_mirror_s1 >> figs_replay.log 2>&1
"%PY%" _figs_20261003.py spinal_run_spkbase_spiking.npz spiking_mirror_ground >> figs_replay.log 2>&1
echo [%DATE% %TIME%] FIGS chain end rc=%ERRORLEVEL% >> CURRICULUM_RUN_20261003.log
