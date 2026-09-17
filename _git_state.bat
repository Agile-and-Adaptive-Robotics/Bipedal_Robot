@echo off
set GIT="C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.5\resources\app\git\cmd\git.exe"
cd /d D:\Github\Bipedal_Robot
echo === STATUS COUNT ===
%GIT% status --short 2>&1 | find /c /v ""
echo === HEAD ===
%GIT% log --oneline -3
echo === A8746F9 ===
%GIT% log --oneline a8746f9 -1 2>&1
echo === REMOTES ===
%GIT% remote -v
echo === BRANCH ===
%GIT% branch --show-current
echo === GIT SIZE ===
dir /s /b .git 2>nul | findstr /c:"." > nul
for /f %%s in ('dir /s /b .git 2^>nul ^| find /c /v ""') do echo .git entries: %%s
