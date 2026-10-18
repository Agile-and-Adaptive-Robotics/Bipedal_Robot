@echo off
rem W2L wean campaign runner (2026-10-16): invoked via WMI Win32_Process.Create
rem so it escapes the sshd job object (nohup/start/Start-Process all die with
rem the ssh session on this box -- measured today).
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_wean
D:\Anaconda\envs\myo\python.exe -u _w2l_wean_campaign.py > campaign.log 2>&1
