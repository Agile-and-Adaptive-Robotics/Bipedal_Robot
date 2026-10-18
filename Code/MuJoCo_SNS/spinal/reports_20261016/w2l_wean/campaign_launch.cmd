@echo off
rem W2L wean campaign launcher (2026-10-16): detach via Windows `start`
rem (git-bash nohup dies with the ssh console on this box -- measured).
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\w2l_wean
start "w2l_wean" /min cmd /c "D:\Anaconda\envs\myo\python.exe -u _w2l_wean_campaign.py > campaign.log 2>&1"
echo LAUNCHED_VIA_START
