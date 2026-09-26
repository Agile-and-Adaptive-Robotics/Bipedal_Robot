# Create scratch copies of Tutorials3 for the optimizations, then warm-start-edit the 4a gait scenario.
$ErrorActionPreference = 'Stop'
$TUT   = 'C:\Users\Ben Bolen\Documents\SCONE\Tutorials3'
$GDIR  = Join-Path $env:TEMP 'scone_opt_gait'
$BDIR  = Join-Path $env:TEMP 'scone_opt_balance'

robocopy $TUT $GDIR /E /NFL /NDL /NJH /NJS | Out-Null
robocopy $TUT $BDIR /E /NFL /NDL /NJH /NJS | Out-Null
Write-Host ("scratch gait   : " + $GDIR)
Write-Host ("scratch balance: " + $BDIR)

# --- 4a warm-start edit: insert init block after the min_progress line ---
$gait = Join-Path $GDIR 'Tutorial 4a - Gait - Hyfydy.scone'
$raw = [System.IO.File]::ReadAllText($gait)
$anchor = "min_progress = 1e-4"
if ($raw -notmatch [regex]::Escape($anchor)) { throw 'anchor not found in gait scenario' }
$init = "$anchor`r`n`t`r`n`tinit { file = par/H0918GaitRS2Hfd4.par std_factor = 2 use_best_as_mean = 1 }"
$raw = $raw.Replace($anchor, $init)   # first (only) occurrence in CmaOptimizer header
[System.IO.File]::WriteAllText($gait, $raw)
Write-Host '--- scratch 4a scenario head ---'
Get-Content $gait -TotalCount 9
