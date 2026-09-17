@echo off
set GIT="C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.5\resources\app\git\cmd\git.exe"
cd /d D:\Github\Bipedal_Robot
echo === TRACKED FILES UNDER SADb_audit ===
%GIT% ls-files SADb_audit 2>&1 | find /c /v ""
echo === DIRTY SADb_audit ENTRIES ===
%GIT% status --short SADb_audit 2>&1 | find /c /v ""
%GIT% status --short SADb_audit 2>&1
echo === BIG SUBDIRS ===
powershell -NoProfile -Command "Get-ChildItem SADb_audit -Directory | ForEach-Object { $s=(Get-ChildItem $_.FullName -Recurse -File -ErrorAction SilentlyContinue | Measure-Object Length -Sum).Sum/1MB; if ($s -gt 20) { '{0,10:N0} MB  {1}' -f $s, $_.Name } }"
echo === BIG FILES AT ROOT ===
powershell -NoProfile -Command "Get-ChildItem SADb_audit -File | Sort-Object Length -Descending | Select-Object -First 6 | ForEach-Object { '{0,10:N1} MB  {1}' -f ($_.Length/1MB), $_.Name }"
