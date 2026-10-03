@echo off
rem Sequential headless AnimatSimulator runs. Usage: run_all.cmd <outdir-tag> <file1.asim> <file2.asim> ...
setlocal
set "BINDIR=D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin"
set "WORK=D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002"
set "TAG=%~1"
shift
:loop
if "%~1"=="" goto end
set "F=%~1"
set "BASE=%~n1"
echo === RUN %TAG% %BASE% ===
"%BINDIR%\AnimatSimulator.exe" "%F%" > "%WORK%\logs\%TAG%_%BASE%.log" 2>&1
echo exitcode=%ERRORLEVEL%
shift
goto loop
:end
endlocal
