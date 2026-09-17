@echo off
set GIT="C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.5\resources\app\git\cmd\git.exe"
set BK=D:\Github\rewrite_backup
cd /d D:\Github\Bipedal_Robot
if not exist %BK% mkdir %BK%
echo === SAVE DIRTY LIST ===
%GIT% status --short > %BK%\pre_filter_status.txt
type %BK%\pre_filter_status.txt | find /c /v ""
echo === COPY DIRTY FILES ===
powershell -NoProfile -Command "Get-Content '%BK%\pre_filter_status.txt' | ForEach-Object { $s=$_.Substring(0,2).Trim(); $p=$_.Substring(3).Trim(); if ($p -like '*(*)') { $p = $p -replace ' -> .*$','' }; $src = Join-Path 'D:\Github\Bipedal_Robot' $p; if (Test-Path $src -PathType Leaf) { $dst = Join-Path '%BK%\worktree' $p; New-Item -ItemType Directory -Force -Path (Split-Path $dst) | Out-Null; Copy-Item $src $dst -Force } }"
powershell -NoProfile -Command "(Get-ChildItem -Recurse -File '%BK%\worktree' | Measure-Object).Count"
echo === BUNDLE ===
%GIT% bundle create %BK%\pre_filter.bundle --all 2>&1
dir %BK%\pre_filter.bundle | findstr bundle
