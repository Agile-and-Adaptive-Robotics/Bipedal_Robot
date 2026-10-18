@echo off
echo START %date% %time% >> D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\wean_human\probe2_task.log
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
D:\Anaconda\envs\myo\python.exe -u wean_probe_20261016.py reports_20261016\wean_human\round8_best3.json 1.0 0.7 0.55 0.45 0.38 0.3 0.2 0.1 0.0 >> D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\wean_human\probe2_task.log 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\wean_human\probe2_task.log
