@echo off
rem Detached driver: curriculum stage 2, syn6 variant (reports_20260925 campaign).
rem Runs 18 TPE trials of curr_syn6_s2_air_aff (fresh study; seeds from the
rem stage-1 winner json on the empty-study branch), then writes the
rem done-marker via write_done_curr_syn6_s2.ps1.
rem Variant selection per goal4_wiring.md: the script pins AARL_NET/AARL_NPZ
rem itself at main() (_curriculum_syn6.py:252-253); AARL_NET is ALSO set here
rem inside the launched command line as belt-and-suspenders (same value).
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set AARL_NET=syn6
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum_syn6.py 2 18 > "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs\curr_syn6_s2.log" 2>&1
set EC=%ERRORLEVEL%
powershell -NoProfile -ExecutionPolicy Bypass -File "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\tmp\write_done_curr_syn6_s2.ps1" -ExitCode %EC%
