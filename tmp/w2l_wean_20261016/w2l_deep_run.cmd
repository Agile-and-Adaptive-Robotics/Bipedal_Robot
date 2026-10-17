@echo off
rem W2L wean DEEP extension runner (WMI-detached, see w2l_wean_run.cmd)
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_wean
D:\Anaconda\envs\myo\python.exe -u _w2l_wean_deep.py > deep.log 2>&1
