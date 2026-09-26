@echo off
rem Detached launcher: w2lvar curriculum stage 3 (18 trials), goal4 wiring.
rem Env selected per goal4_wiring.md: script pins AARL_NET/AARL_NPZ itself;
rem AARL_NET is also set here inside the launched command line (harmless,
rem identical value). Concurrent with the syn6 chain (own db + npz).
cd /d "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
set "AARL_NET=w2lvar"
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum_w2lvar.py 3 18 > "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs\curr_w2lvar_s3.log" 2>&1
set "EC=%ERRORLEVEL%"
>"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs\curr_w2lvar_s3.done" (
 echo %EC%
 powershell -NoProfile -Command "Get-Content -LiteralPath 'D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs\curr_w2lvar_s3.log' -Tail 3"
)
