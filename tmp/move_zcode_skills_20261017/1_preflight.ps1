# Preflight inventory for the ZCode_Skills D:->C: move (2026-10-17, easteregg2)
$old = 'D:\GitHub\ZCode_Skills'

function Inv($p) {
  if (Test-Path -LiteralPath $p) {
    Write-Output "== $p"
    Get-ChildItem -LiteralPath $p -Force | ForEach-Object {
      Write-Output ("  " + $_.Name + " | LinkType=" + $_.LinkType + " | Target=" + ($_.Target -join ','))
    }
  } else { Write-Output "== $p MISSING" }
}

Inv 'C:\Users\Ben Bolen\.agents\skills'
Inv 'C:\Users\Ben Bolen\.zcode\skills'
Inv 'C:\Users\Ben Bolen\.zcode\cli\plugins\cache\dev-zcode-skills'

# Any reparse points INSIDE the repo itself? (robocopy traversal decision)
Write-Output "== reparse points inside $old (recursive)"
$reps = Get-ChildItem -LiteralPath $old -Recurse -Force -ErrorAction SilentlyContinue | Where-Object { $_.LinkType }
if ($reps) { $reps | ForEach-Object { Write-Output ("  " + $_.FullName + " -> " + ($_.Target -join ',')) } } else { Write-Output "  none" }
