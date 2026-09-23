@echo off
rem Walker env relay (mujoco-sns-walker plugin).
rem Picks the SNS conda python across the three AARL machines, sets
rem CONDA_PREFIX (MuJoCo reads it), and forwards to a plugin script.
rem Usage: walker.cmd <script.py> [args...]
setlocal
set "PY=%AARL_PYTHON%"
if "%PY%"=="" set "PY=C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
if not exist "%PY%" set "PY=D:\Anaconda\envs\myo\python.exe"
if not exist "%PY%" set "PY=C:\Users\Ben\.anaconda3\envs\myoconv\python.exe"
if not exist "%PY%" (
  echo [walker] no SNS python found - set AARL_PYTHON to the env's python.exe 1>&2
  exit /b 1
)
for %%I in ("%PY%") do set "CONDA_PREFIX=%%~dpI"
if "%CONDA_PREFIX:~-1%"=="\" set "CONDA_PREFIX=%CONDA_PREFIX:~0,-1%"
"%PY%" "%~dp0%~nx1" %2 %3 %4 %5 %6 %7 %8 %9
endlocal
