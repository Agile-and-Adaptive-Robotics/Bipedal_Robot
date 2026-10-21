@echo off
rem 2026-10-20 wave-2 balance+perturbation batteries: W2L Hill walker
rem (latch feet-clear AND step feet-loaded, rig 0.25) + BPA walker
rem (march drive, support 0.70/0.85). 7 cases each.
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco
set PY=D:\Anaconda\envs\myo\python.exe
set LOG=balance_w2l_20261020.log
echo START %date% %time% >> %LOG%
echo === HILL latch rig0.25 === >> %LOG%
%PY% -u test_w2l_ground.py --phase=stand --perturb=1 --perturb-rig=0.25 --dur_stand=12 >> %LOG% 2>&1
echo === HILL step rig0.25 === >> %LOG%
%PY% -u test_w2l_ground.py --phase=stand --perturb=1 --perturb-rig=0.25 --dur_stand=12 --sregime=step >> %LOG% 2>&1
echo === HILL step rig0.0 (unrigged attempt) === >> %LOG%
%PY% -u test_w2l_ground.py --phase=stand --perturb=1 --perturb-rig=0.0 --dur_stand=12 --sregime=step >> %LOG% 2>&1
echo === BPA march support0.70 === >> %LOG%
%PY% -u test_bpa_stepping.py --balance=1 --dur=12 --hold=3 --lower=1 --stand-te=3 --stand-tf=4 --support=0.70 --pushN=40 --pushAxis=+x >> %LOG% 2>&1
%PY% -u test_bpa_stepping.py --balance=1 --dur=12 --hold=3 --lower=1 --stand-te=3 --stand-tf=4 --support=0.70 --pushN=40 --pushAxis=-x >> %LOG% 2>&1
%PY% -u test_bpa_stepping.py --balance=1 --dur=12 --hold=3 --lower=1 --stand-te=3 --stand-tf=4 --support=0.70 --pushN=40 --pushAxis=+y >> %LOG% 2>&1
%PY% -u test_bpa_stepping.py --balance=1 --dur=12 --hold=3 --lower=1 --stand-te=3 --stand-tf=4 --support=0.70 --pushN=40 --pushAxis=-y >> %LOG% 2>&1
%PY% -u test_bpa_stepping.py --balance=1 --dur=12 --hold=3 --lower=1 --stand-te=3 --stand-tf=4 --support=0.70 --platDx=0.05 >> %LOG% 2>&1
%PY% -u test_bpa_stepping.py --balance=1 --dur=12 --hold=3 --lower=1 --stand-te=3 --stand-tf=4 --support=0.70 --platTilt=5 >> %LOG% 2>&1
echo === BPA march support0.85 (push_f only) === >> %LOG%
%PY% -u test_bpa_stepping.py --balance=1 --dur=12 --hold=3 --lower=1 --stand-te=3 --stand-tf=4 --support=0.85 --pushN=40 --pushAxis=+x >> %LOG% 2>&1
echo EXITCODE %ERRORLEVEL% %date% %time% >> %LOG%
