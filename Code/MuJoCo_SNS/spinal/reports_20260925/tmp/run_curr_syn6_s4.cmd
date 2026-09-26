@echo off
rem Detached driver: curriculum stage 4, syn6 variant (reports_20260925 campaign).
rem Runs 22 TPE trials of curr_syn6_s4_walk_nocontact, resuming the study that
rem holds the 1-trial route-proof smoke from the wiring session (trial 0,
rem score -254.936); study non-empty -> no new seed enqueued. Writes the
rem done-marker via write_done_curr_syn6_s4.ps1.
rem Variant selection per goal4_wiring.md: the script pins AARL_NET/AARL_NPZ
rem itself at main() (_curriculum_syn6.py:252-253); AARL_NET is ALSO set here
rem inside the launched command line as belt-and-suspenders (same value).
cd /d D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set AARL_NET=syn6
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" _curriculum_syn6.py 4 22 > "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs\curr_syn6_s4.log" 2>&1
set EC=%ERRORLEVEL%
powershell -NoProfile -ExecutionPolicy Bypass -File "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\tmp\write_done_curr_syn6_s4.ps1" -ExitCode %EC%
