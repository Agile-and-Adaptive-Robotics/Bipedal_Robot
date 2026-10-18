@echo off
echo START %date% %time% >> D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\drive_sweep_human\sweep.log
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
if not exist reports_20261016\drive_sweep_human mkdir reports_20261016\drive_sweep_human
D:\Anaconda\envs\myo\python.exe -u drive_sweep_human_20261016.py >> reports_20261016\drive_sweep_human\sweep.log 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\drive_sweep_human\sweep.log
