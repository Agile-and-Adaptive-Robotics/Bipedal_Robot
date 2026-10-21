@echo off
rem 2026-10-21 resume session, W2L Hill robot: Ib ablation (ibnA=0) +
rem fresh intact-stand confirmation, both sregimes at rig 0.25.
cd /d D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco
set PY=D:\Anaconda\envs\myo\python.exe
set LOG=balance_w2l_ibablate_20261021.log
echo START %date% %time% >> %LOG%
echo === HILL latch rig0.25 INTACT (fresh confirm) === >> %LOG%
%PY% -u test_w2l_ground.py --phase=stand --dur-stand=12 >> %LOG% 2>&1
echo RC %ERRORLEVEL% >> %LOG%
echo === HILL step rig0.25 INTACT (fresh confirm) === >> %LOG%
%PY% -u test_w2l_ground.py --phase=stand --dur-stand=12 --sregime=step >> %LOG% 2>&1
echo RC %ERRORLEVEL% >> %LOG%
echo === HILL latch rig0.25 IB-ABLATED (ibnA=0) === >> %LOG%
%PY% -u test_w2l_ground.py --phase=stand --dur-stand=12 --ibnA=0 >> %LOG% 2>&1
echo RC %ERRORLEVEL% >> %LOG%
echo === HILL step rig0.25 IB-ABLATED (ibnA=0) === >> %LOG%
%PY% -u test_w2l_ground.py --phase=stand --dur-stand=12 --sregime=step --ibnA=0 >> %LOG% 2>&1
echo RC %ERRORLEVEL% >> %LOG%
echo === HILL latch rig0.25 IB-ABLATED push battery (perturb) === >> %LOG%
%PY% -u test_w2l_ground.py --phase=stand --perturb=1 --perturb-rig=0.25 --dur-stand=12 --ibnA=0 >> %LOG% 2>&1
echo RC %ERRORLEVEL% >> %LOG%
echo EXITCODE 0 %date% %time% >> %LOG%
