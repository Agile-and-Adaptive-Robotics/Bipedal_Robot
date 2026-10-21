@echo off
rem 2026-10-21 stage-3Q FIXED-OBJECTIVE study - syn6 (the family whose
rem 2026-10-20 s3p winner gamed the score by sagging onto the slack
rem rig: feet 29.8% BW, COM 0.575 m). REAL budget: 60 trials x 7 evals.
rem New study name curr_syn6_s3q_perturb - the recorded
rem curr_syn6_s3p_perturb study is untouched. Seeds: trial 0 = quiet
rem incumbent defaults, trial 1 = the s3p "wean" winner params
rem (curriculum_s3p_syn6.json) - the pose the rig-share sentinel must
rem now expose.
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal
set PYTHONIOENCODING=utf-8
echo LAUNCH-START %date% %time% >> s3q_syn6_20261021.log
D:\Anaconda\envs\myo\python.exe -u _curriculum_s3p.py syn6 60 >> s3q_syn6_20261021.log 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> s3q_syn6_20261021.log
