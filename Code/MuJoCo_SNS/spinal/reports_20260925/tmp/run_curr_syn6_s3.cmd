@echo off
rem Detached driver: curriculum stage 3, syn6 variant (reports_20260925 campaign).
rem Runs 18 TPE trials of curr_syn6_s3_balance (fresh study; seeds
rem vest_prop=0.0 / rig_scale=1.0 - stage-2 winner has neither key, so the
rem chain-merge is a no-op here by design), then writes the done-marker via
rem write_done_curr_syn6_s3.ps1.
rem Variant selection per goal4_wiring.md: the script pins AARL_NET/AARL_NPZ
rem itself at main() (_curriculum_syn6.py:252-253); AARL_NET is ALSO set here
rem inside the launched command line as belt-and-suspenders (same value).
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set AARL_NET=syn6
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum_syn6.py 3 18 > "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs\curr_syn6_s3.log" 2>&1
set EC=%ERRORLEVEL%
powershell -NoProfile -ExecutionPolicy Bypass -File "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\tmp\write_done_curr_syn6_s3.ps1" -ExitCode %EC%
