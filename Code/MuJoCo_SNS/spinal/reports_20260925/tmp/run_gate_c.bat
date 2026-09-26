@echo off 
set AARL_NET=w2lvar 
set AARL_NPZ=w2lvar_smoke.npz 
"C:\Users\Ben Bolen\.conda\envs\myo\python.exe" runner.py --no-ground --drive 2.5 %*
