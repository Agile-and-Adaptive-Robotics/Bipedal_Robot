@echo off
rem Routing27 campaign launcher (easteregg2)
rem Usage: run_routing27.cmd [plan^|validate^|full^|smoke] [surrogateBudget] [patternBudget]
rem Blocks while MATLAB runs; launch detached with:
rem   start "" /b cmd /c "D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27\run_routing27.cmd full" ^> console.log 2^>^&1
setlocal
set MODE=%1
if "%MODE%"=="" set MODE=full
if not "%2"=="" set ROUTING27_SURROGATE=%2
if not "%3"=="" set ROUTING27_PATTERN=%3
set ROUTING27_MODE=%MODE%
set BASE=D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27
if not exist "%BASE%\Results" mkdir "%BASE%\Results"
echo [%date% %time%] launching routing27 %MODE% >> "%BASE%\Results\routing27_console_%MODE%.log"
"D:\Program Files\MATLAB\R2025a\bin\matlab.exe" -batch "run('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Routing27/Opt_run_Routing27.m')" >> "%BASE%\Results\routing27_console_%MODE%.log" 2>&1
echo [%date% %time%] routing27 %MODE% finished, exit %ERRORLEVEL% >> "%BASE%\Results\routing27_console_%MODE%.log"
endlocal
