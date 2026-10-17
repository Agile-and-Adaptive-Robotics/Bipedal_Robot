@echo off
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
echo FULLRUN START %date% %time% >> reports_20261016\drive_sweep_human\fullrun.log
D:\Anaconda\envs\myo\python.exe -u drive_sweep_full_20261016.py reports_20261016\drive_sweep_human\winner_g1.5_a1_d0.json reports_20261016\drive_sweep_human\winner_g1.5_a1_d0_full.npz >> reports_20261016\drive_sweep_human\fullrun.log 2>&1
D:\Anaconda\envs\myo\python.exe -u drive_sweep_full_20261016.py reports_20261016\drive_sweep_human\runnerup_g2.5_a1_d0.json reports_20261016\drive_sweep_human\runnerup_g2.5_a1_d0_full.npz >> reports_20261016\drive_sweep_human\fullrun.log 2>&1
echo FULLRUN EXITCODE %ERRORLEVEL% %date% %time% >> reports_20261016\drive_sweep_human\fullrun.log
