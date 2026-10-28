@echo off
powershell.exe -NoProfile -File "%~dp0Launch_AH001_Replay.ps1" %*
if errorlevel 1 pause
