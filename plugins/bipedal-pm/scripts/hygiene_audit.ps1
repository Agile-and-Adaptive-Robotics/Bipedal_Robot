<#
  bipedal-pm hygiene_audit.ps1 - READ-ONLY. Never deletes, moves, or renames anything.

  Modes:
    audit                  full audit (default): unexpected repo-root entries vs the
                           AGENTS.md canonical map, suspicious names repo-wide,
                           zero-byte files, large files, advisory notes.
    refs -Name <name>      reference-impact grep for a proposed move/rename target
                           (text extensions only; skips .git, reparse points, >20 MB files).

  Repo-root resolution: BIPEDAL_REPO env -> cwd -> known machine candidates
  (EB475WS4 D:\Github\Bipedal_Robot, easteregg2 D:\GitHub\Bipedal_Robot,
   laptop C:\Users\Ben\Documents\GitHub\Bipedal_Robot).
  Exit codes: 0 ok, 1 usage error, 2 repo not found.
#>
param(
    [Parameter(Position = 0)][string]$Mode = 'audit',
    [string]$Name = '',
    [long]$LargeBytes = 50MB,
    [int]$MaxDepth = 5,
    [int]$MaxList = 40
)
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

function Get-RepoFiles {
    # Stack walk so directory reparse points (junctions) are never followed and
    # skip-listed directory names are pruned before descending.
    param([string]$Root, [string[]]$SkipNames, [int]$MaxDepth)
    $stack = New-Object System.Collections.Stack
    $rootItem = Get-Item -LiteralPath $Root -Force
    $stack.Push(@($rootItem, 0))
    while ($stack.Count -gt 0) {
        $entry = $stack.Pop(); $dir = $entry[0]; $depth = $entry[1]
        $children = Get-ChildItem -LiteralPath $dir.FullName -Force -ErrorAction SilentlyContinue
        if (-not $children) { continue }
        foreach ($c in $children) {
            if ($c.PSIsContainer) {
                if ($c.Attributes -band [IO.FileAttributes]::ReparsePoint) { continue }
                if ($SkipNames -contains $c.Name) { continue }
                if ($depth -lt $MaxDepth) { $stack.Push(@($c, $depth + 1)) }
            } else {
                $c
            }
        }
    }
}

function Fmt-Size([long]$bytes) {
    if ($bytes -ge 1GB) { return ('{0:N1} GB' -f ($bytes / 1GB)) }
    if ($bytes -ge 1MB) { return ('{0:N1} MB' -f ($bytes / 1MB)) }
    if ($bytes -ge 1KB) { return ('{0:N1} KB' -f ($bytes / 1KB)) }
    return ('{0} B' -f $bytes)
}

$repo = Resolve-Repo
if (-not $repo) {
    Write-Output '[bipedal-pm] ERROR: repo root not found (set BIPEDAL_REPO).'
    exit 2
}
$stamp = Get-Date -Format 'yyyy-MM-dd HH:mm'
Write-Output ('== bipedal-pm hygiene @' + $env:COMPUTERNAME + '  ' + $stamp + '  repo=' + $repo + ' (READ-ONLY) ==')

if ($Mode -eq 'refs') {
    if (-not $Name) { Write-Output 'refs mode needs -Name <file-or-folder-name>'; exit 1 }
    if ($Name -match '^[\d\W]{1,3}$') {
        Write-Output ('[refs] WARNING: name "' + $Name + '" is too generic for content grep; treat as filename-only artifact and decide disposition directly.')
    }
    $skip = @('.git', '__pycache__', 'slprj', '.zcode_tmp', 'node_modules')
    $exts = @('.m', '.py', '.md', '.json', '.bat', '.ps1', '.cmd', '.tex', '.xml', '.html',
              '.aproj', '.asim', '.cfg', '.yml', '.yaml', '.bib', '.cls', '.sty', '.gitignore', '.zcodeignore')
    $all = @(Get-RepoFiles -Root $repo -SkipNames $skip -MaxDepth $MaxDepth)
    $cands = $all | Where-Object {
        ($exts -contains $_.Extension.ToLower()) -and ($_.Length -lt 20MB)
    }
    $hits = @($cands | Select-String -SimpleMatch -Pattern $Name -ErrorAction SilentlyContinue)
    if ($hits.Count -eq 0) {
        Write-Output ('[refs] "' + $Name + '": 0 text references found in ' + $cands.Count + ' candidate files.')
    } else {
        $byFile = $hits | Group-Object Path | Sort-Object Count -Descending
        Write-Output ('[refs] "' + $Name + '": ' + $hits.Count + ' hit line(s) in ' + $byFile.Count + ' file(s) of ' + $cands.Count + ' scanned:')
        $i = 0
        foreach ($g in $byFile) {
            if ($i -ge $MaxList) { Write-Output ('  ... ' + ($byFile.Count - $MaxList) + ' more file(s)'); break }
            $tag = ''
            if ($g.Name -match 'AGENTS\.md$') { $tag = '   <- AGENTS.md directory map: MUST update on any move' }
            if ($g.Name -match '\.vscode') { $tag = '   <- .vscode: per-machine config' }
            Write-Output ('  {0,4}x  {1}{2}' -f $g.Count, $g.Name, $tag)
            $i++
        }
    }
    Write-Output '[refs] Protocol: fix EVERY hit class in the same pass, update the AGENTS.md map, then smoke-verify one affected script. Ben approves the move before it happens.'
    exit 0
}

# ---------------- audit mode ----------------
$skip = @('.git', '__pycache__', 'slprj', '.zcode_tmp', 'node_modules')
$all = @(Get-RepoFiles -Root $repo -SkipNames $skip -MaxDepth $MaxDepth)

# --- 1. repo-root inventory vs canonical map (AGENTS.md directory map) ---
$expectedRoot = @(
    'AGENTS.md', 'README.md', 'PROJECT_INSTRUCTIONS.md', 'CHATGPT_HANDOFF.md', 'CHATGPT_REPORT.md',
    'Code', 'Documentation', 'Neuromechanical_Models', 'Pictures', 'Solid_Models', 'Testing_Data',
    'SADb_audit', 'plugins',
    '.git', '.gitignore', '.vscode', '.zcode', '.zcodeignore', '.zcode_tmp', '.zcode-plugins'
)
$hint = @{
    '0'                                   = 'connectome NODES/EDGES layout dump (2026-09-24, 26n/46e, edge TARGETS missing) - keep-but-relocate candidate (spinal\), Ben rules'
    '-p'                                  = 'EMPTY dir - shell "mkdir -p" typo artifact - safe delete candidate'
    '.DS_Store'                           = 'macOS Finder artifact - delete candidate'
    'MUJOCO_LOG.TXT'                      = 'MuJoCo run log at root - delete candidate or gitignore'
    'opensim.log'                         = 'OpenSim run log at root - already gitignored (/opensim.log) - delete candidate'
    'slprj'                               = 'Simulink cache at ROOT (gitignore only covers SNS_Simscape paths) - gitignore + delete candidate'
    'temp'                                = 'root temp dir - review contents, then delete/archive'
    'spring_series.m'                     = 'stray MATLAB file at root - find its owner tree or archive'
    "Jeffrey's practice pcb"              = "personal folder - Ben's call"
    'CHATGPT_REPORT.md.bak_handoff_20260909' = 'dated backup of CHATGPT_REPORT.md - archive candidate'
}
$rootEntries = @(Get-ChildItem -LiteralPath $repo -Force -ErrorAction SilentlyContinue)
$unexpected = @($rootEntries | Where-Object { $expectedRoot -notcontains $_.Name })
if ($unexpected.Count -eq 0) {
    Write-Output '-- Repo root: clean (all entries match the AGENTS.md canonical map)'
} else {
    Write-Output ('-- Repo root: ' + $unexpected.Count + ' entr(ies) NOT in the AGENTS.md canonical map:')
    foreach ($e in $unexpected) {
        $kind = 'DIR '; if (-not $e.PSIsContainer) { $kind = 'FILE' }
        $sz = ''; if (-not $e.PSIsContainer) { $sz = Fmt-Size $e.Length + '  ' }
        $h = $hint[$e.Name]; if (-not $h) { $h = '' }
        if ($e.Name -like '_git_*.bat') { $h = 'one-off git helper .bat - archive candidate' }
        if (-not $h -and $e.Name -match '\.bak') { $h = 'backup file - archive candidate' }
        Write-Output ('   {0}  {1}  {2}  {3}  {4}' -f $kind, $sz.PadRight(9), $e.LastWriteTime.ToString('yyyy-MM-dd'), $e.Name, $h)
    }
}

# --- 2. suspicious names repo-wide ---
$patterns = @('^0$', '^-p$', '^~\$', '\(2\)\.?', '- Copy', '\.DS_Store$', '^Thumbs\.db$', '^desktop\.ini$',
              '\.bak$', '\.orig$', '\.tmp$', '^nul$', '^ Untitled', '^New Folder')
# patterns are already regex strings; join them into one alternation
$rx = [string]::Join('|', $patterns)
$susp = @($all | Where-Object { $_.Name -match $rx })
Write-Output ('-- Suspicious names repo-wide: ' + $susp.Count + ' file(s)/dir(s) matched')
$i = 0
foreach ($f in ($susp | Select-Object -First $MaxList)) {
    Write-Output ('   ' + $f.FullName.Substring($repo.Length + 1) + '  (' + (Fmt-Size $f.Length) + ', ' + $f.LastWriteTime.ToString('yyyy-MM-dd') + ')')
    $i++
}
if ($susp.Count -gt $MaxList) { Write-Output ('   ... ' + ($susp.Count - $MaxList) + ' more') }

# --- 3. zero-byte files ---
$zeros = @($all | Where-Object { $_.Length -eq 0 })
Write-Output ('-- Zero-byte files: ' + $zeros.Count)
$i = 0
foreach ($f in ($zeros | Select-Object -First 15)) {
    Write-Output ('   ' + $f.FullName.Substring($repo.Length + 1) + '  (' + $f.LastWriteTime.ToString('yyyy-MM-dd') + ')')
    $i++
}
if ($zeros.Count -gt 15) { Write-Output ('   ... ' + ($zeros.Count - 15) + ' more') }

# --- 4. large files ---
$big = @($all | Where-Object { $_.Length -gt $LargeBytes } | Sort-Object Length -Descending)
Write-Output ('-- Large files (> ' + (Fmt-Size $LargeBytes) + '): ' + $big.Count)
$i = 0
foreach ($f in ($big | Select-Object -First 30)) {
    $disp = 'data/result - if needed long-term, propose moving OUT of repo (D:\temp\<name> or archive); ask Ben'
    if ($f.FullName -match 'myosuite_gait2392|LocalBuild|\.zcode') { $disp = 'regenerable (gitignored) - may delete locally' }
    Write-Output ('   {0,9}  {1}  {2}' -f (Fmt-Size $f.Length), $f.FullName.Substring($repo.Length + 1), $disp)
    $i++
}
if ($big.Count -gt 30) { Write-Output ('   ... ' + ($big.Count - 30) + ' more') }

# --- 5. advisories ---
$codeSpaces = @($all | Where-Object { $_.FullName -match '\\Code\\' -and $_.Name -match ' ' })
Write-Output ('-- Advisory: ' + $codeSpaces.Count + ' file(s) under Code\ contain spaces in the name (convention drift; mass renames = POST-DEFENSE)')
Write-Output '-- Advisory: .zcode_tmp and temp contents were skipped from the walk - review them separately if flagged at root.'
Write-Output '-- This audit never modifies anything. Dispositions require Ben''s explicit per-item approval; any move runs the refs-impact protocol first.'
exit 0
