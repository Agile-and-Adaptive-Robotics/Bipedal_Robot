@echo off
set GIT="C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.5\resources\app\git\cmd\git.exe"
cd /d D:\Github\Bipedal_Robot
echo === NPZ IN TIP COMMIT ===
%GIT% show --stat a8746f9 2>&1 | find /i /c ".npz"
echo === NPZ EVER (all commits touching npz) ===
%GIT% log --oneline -- "*.npz" 2>&1
echo === BRANCH VS REMOTE ===
%GIT% status -sb 2>&1
echo === TOP BLOBS IN TIP ===
%GIT% ls-tree -r -l a8746f9 2>&1 | sort /r /+55 | more +0 2>nul
echo === ALL BRANCHES ===
%GIT% branch -a 2>&1
