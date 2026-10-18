# Move ZCode_Skills from D:\GitHub to C:\Users\Ben Bolen\Documents\GitHub (2026-10-17, easteregg2)
$src = 'D:\GitHub\ZCode_Skills'
$dst = 'C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills'

if (Test-Path -LiteralPath $dst) { Write-Output 'ABORT: destination already exists'; exit 1 }
if (-not (Test-Path -LiteralPath "$src\.git")) { Write-Output 'ABORT: source is not a git repo'; exit 1 }

robocopy $src $dst /E /MOVE /R:2 /W:2 /NFL /NDL /NP
$rc = $LASTEXITCODE
Write-Output "ROBOCOPY_RC=$rc"
if ($rc -ge 8) { Write-Output 'ABORT: robocopy failed (>=8)'; exit 1 }

if (Test-Path -LiteralPath $src) {
  Write-Output 'WARN: source remnants remain:'
  Get-ChildItem -LiteralPath $src -Recurse -Force -ErrorAction SilentlyContinue | Select-Object -First 10 | ForEach-Object { Write-Output ('  ' + $_.FullName) }
} else {
  Write-Output 'SOURCE-REMOVED'
}

Write-Output 'GIT_HEAD:'
& git -C $dst rev-parse HEAD
Write-Output 'GIT_STATUS:'
& git -C $dst status --porcelain
Write-Output 'MOVE-OK'
