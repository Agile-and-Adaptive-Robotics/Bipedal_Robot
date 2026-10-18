@echo off
echo START %date% %time% >> D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\wean_human\r1_task.log
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
D:\Anaconda\envs\myo\python.exe -u wean_round_launch_20261016.py 1 25 4 >> D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\wean_human\r1_task.log 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\wean_human\r1_task.log
