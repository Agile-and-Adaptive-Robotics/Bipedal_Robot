@echo off
rem Double-click to watch the gait2392 spinal walk in a 3D viewer.
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" runner.py --view %*
echo.
echo ===== done: results in spinal_run.png / spinal_run.npz =====
pause
