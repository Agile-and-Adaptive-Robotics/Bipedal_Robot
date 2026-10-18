$ErrorActionPreference = 'Stop'
function Get-Words($text) {
  # strip LaTeX comments, then commands/braces, then tokenize
  $t = ($text -split "`n") | Where-Object { $_ -notmatch '^\s*%' } | Out-String
  $t = $t -replace '\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})*', ' '
  $t = $t -replace '[{}~$&\\]', ' '
  return [regex]::Matches($t.ToLower(), '[a-z0-9]+') | ForEach-Object { $_.Value }
}
$ms = Get-Content -Raw 'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Isometric_Force_Test\Final\Final_proofread\manuscript.v5\actuators-4207821-edited.tex'
$mw = @(Get-Words $ms)
$msGrams = New-Object 'System.Collections.Generic.HashSet[string]'
for ($i = 0; $i -le $mw.Count - 10; $i++) {
  [void]$msGrams.Add(($mw[$i..($i+9)] -join ' '))
}
Write-Output ("manuscript words: " + $mw.Count + "; 10-grams: " + $msGrams.Count)

$intro = Get-Content 'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\10-introduction.tex'
# changed lines (new side of diff): 8,10,12,14,18,20,24,27,28,29,30,31,35,37,39,41,43,45,47,49
$changed = @(8,10,12,14,18,20,24,27,28,29,30,31,35,37,39,41,43,45,47,49)
$hits = 0
foreach ($ln in $changed) {
  $w = @(Get-Words $intro[$ln-1])
  for ($i = 0; $i -le $w.Count - 10; $i++) {
    $g = ($w[$i..($i+9)] -join ' ')
    if ($msGrams.Contains($g)) {
      $hits++
      Write-Output ("LINE " + $ln + " 10-gram match: " + $g)
    }
  }
}
Write-Output ("total 10-gram matches in changed lines: " + $hits)
# also 7-gram sweep over the same changed lines (external-verbatim proxy)
$ms7 = New-Object 'System.Collections.Generic.HashSet[string]'
for ($i = 0; $i -le $mw.Count - 7; $i++) { [void]$ms7.Add(($mw[$i..($i+6)] -join ' ')) }
$h7 = 0
foreach ($ln in $changed) {
  $w = @(Get-Words $intro[$ln-1])
  for ($i = 0; $i -le $w.Count - 7; $i++) {
    $g = ($w[$i..($i+6)] -join ' ')
    if ($ms7.Contains($g)) { $h7++; Write-Output ("LINE " + $ln + " 7-gram match: " + $g) }
  }
}
Write-Output ("total 7-gram matches in changed lines: " + $h7)
