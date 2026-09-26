# bipedal-pm SessionStart banner. Read-only; silent on any failure; always exits 0.
# Prints a 3-6 line reminder: dissertation priority, latest dissertation edit,
# CHATGPT file ages, and any known root strays. Full audit lives in /hygiene.
$ErrorActionPreference = 'SilentlyContinue'

function Resolve-Repo {
    $cands = @()
    if ($env:BIPEDAL_REPO) { $cands += $env:BIPEDAL_REPO }
    $cands += @((Get-Location).Path, 'D:\Github\Bipedal_Robot', 'D:\GitHub\Bipedal_Robot', 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot')
    foreach ($c in $cands) {
        if ($c -and (Test-Path -LiteralPath (Join-Path $c 'AGENTS.md')) -and (Test-Path -LiteralPath (Join-Path $c 'Code'))) { return $c }
    }
    return $null
}

$repo = Resolve-Repo
if ($repo) {
    $out = New-Object System.Collections.Generic.List[string]
    $diss = Join-Path $repo 'Documentation\Reports and Papers\Dissertation'

    $line = '[bipedal-pm] PRIORITY 1: dissertation (deadline Sept 2026'
    $dd = Join-Path $diss 'Notes\DEFENSE_DATE.txt'
    if (Test-Path -LiteralPath $dd) {
        $txt = (Get-Content -LiteralPath $dd -TotalCount 1 -ErrorAction SilentlyContinue)
        if ($txt) {
            $d = [datetime]::MinValue
            if ([datetime]::TryParse(([string]$txt).Trim(), [ref]$d)) {
                $days = ($d - (Get-Date)).Days
                $line += ('; defense in {0} d' -f $days)
            }
        }
    }
    $line += ') - route every task through it (/diss /digest /handoff /hygiene)'
    $out.Add($line)

    $tex = Get-ChildItem -LiteralPath (Join-Path $diss 'ProofFinal') -Recurse -Filter *.tex -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending | Select-Object -First 1
    if ($tex) {
        $age = ((Get-Date) - $tex.LastWriteTime).Days
        $ageStr = $(if ($age -le 0) { 'today' } else { '{0} d ago' -f $age })
        $out.Add(('[bipedal-pm] Latest dissertation edit: {0} ({1})' -f $tex.Name, $ageStr))
    }

    foreach ($f in @('CHATGPT_HANDOFF.md', 'CHATGPT_REPORT.md')) {
        $p = Join-Path $repo $f
        if (Test-Path -LiteralPath $p) {
            $age = ((Get-Date) - (Get-Item -LiteralPath $p).LastWriteTime).Days
            $ageStr = $(if ($age -le 0) { 'today' } else { '{0} d ago' -f $age })
            $out.Add(('[bipedal-pm] {0}: updated {1}' -f $f, $ageStr))
        }
    }

    $flags = @('0', '-p', '.DS_Store', 'MUJOCO_LOG.TXT', 'opensim.log', 'slprj', 'temp', '.zcode_tmp')
    $found = @($flags | Where-Object { Test-Path -LiteralPath (Join-Path $repo $_) })
    if ($found.Count -gt 0) {
        $out.Add(('[bipedal-pm] Root strays present: {0} - run /hygiene for the audit + disposition' -f ($found -join ', ')))
    }
    $out | ForEach-Object { Write-Output $_ }
}
exit 0
