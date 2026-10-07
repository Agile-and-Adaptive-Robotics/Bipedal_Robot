@echo off
rem Routing27_v2 campaign launcher (easteregg2) - penalty restructure
rem Usage: run_routing27_v2.cmd [plan^|validate^|full^|smoke] [surrogateBudget] [patternBudget]
rem Blocks while MATLAB runs; launch detached with:
rem   start "" /b cmd /c "D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\run_routing27_v2.cmd full" ^> console_v2.log 2^>^&1
setlocal
set MODE=%1
if "%MODE%"=="" set MODE=full
if not "%2"=="" set ROUTING27_SURROGATE=%2
if not "%3"=="" set ROUTING27_PATTERN=%3
set ROUTING27_MODE=%MODE%
set BASE=D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2
if not exist "%BASE%\Results" mkdir "%BASE%\Results"
echo [%date% %time%] launching routing27_v2 %MODE% >> "%BASE%\Results\routing27_console_%MODE%.log"
"D:\Program Files\MATLAB\R2025a\bin\matlab.exe" -batch "run('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Routing27_v2/Opt_run_Routing27.m')" >> "%BASE%\Results\routing27_console_%MODE%.log" 2>&1
echo [%date% %time%] routing27_v2 %MODE% finished, exit %ERRORLEVEL% >> "%BASE%\Results\routing27_console_%MODE%.log"
endlocal
