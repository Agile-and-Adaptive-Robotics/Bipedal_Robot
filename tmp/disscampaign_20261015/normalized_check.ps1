$ol_root = 'D:\Github\Bipedal_Robot\tmp\disscampaign_20261015\overleaf_pull'
$pf_root = 'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal'

# Normalize: strip trailing whitespace per line, drop trailing blank lines, drop graphicspath lines (known env diff)
function Normalize($path) {
  $t = [System.IO.File]::ReadAllText($path)
  $t = $t -replace "`r`n", "`n"
  $lines = $t -split "`n"
  $out = @()
  foreach ($l in $lines) {
    $trim = $l.TrimEnd()
    if ($trim -match '^\\graphicspath') { continue }
    $out += $trim
  }
  while ($out.Count -gt 0 -and $out[-1] -eq '') { $out = $out[0..($out.Count-2)] }
  return ($out -join "`n")
}

$pairs = @(
  @{n='main.tex'; ol="$ol_root\main.tex"; pf="$pf_root\main.tex"},
  @{n='10-introduction.tex'; ol="$ol_root\chapters\10-introduction.tex"; pf="$pf_root\chapters\10-introduction.tex"},
  @{n='15-background.tex'; ol="$ol_root\chapters\15-background.tex"; pf="$pf_root\chapters\15-background.tex"},
  @{n='20-methods.tex'; ol="$ol_root\chapters\20-methods.tex"; pf="$pf_root\chapters\20-methods.tex"},
  @{n='30-results.tex'; ol="$ol_root\chapters\30-results.tex"; pf="$pf_root\chapters\30-results.tex"},
  @{n='40-discussion.tex'; ol="$ol_root\chapters\40-discussion.tex"; pf="$pf_root\chapters\40-discussion.tex"},
  @{n='50-futurework.tex'; ol="$ol_root\chapters\50-futurework.tex"; pf="$pf_root\chapters\50-futurework.tex"},
  @{n='60-conclusion.tex'; ol="$ol_root\chapters\60-conclusion.tex"; pf="$pf_root\chapters\60-conclusion.tex"},
  @{n='90-bib.tex'; ol="$ol_root\chapters\90-bib.tex"; pf="$pf_root\chapters\90-bib.tex"},
  @{n='92-AppendixA.tex'; ol="$ol_root\chapters\92-AppendixA.tex"; pf="$pf_root\chapters\92-AppendixA.tex"},
  @{n='93-AppendixB.tex'; ol="$ol_root\chapters\93-AppendixB.tex"; pf="$pf_root\chapters\93-AppendixB.tex"},
  @{n='94-AppendixC.tex'; ol="$ol_root\chapters\94-AppendixC.tex"; pf="$pf_root\chapters\94-AppendixC.tex"},
  @{n='thesis.bib'; ol="$ol_root\thesis.bib"; pf="$pf_root\thesis.bib"},
  @{n='BolenFrontiers22.bib'; ol="$ol_root\BolenFrontiers22.bib"; pf="$pf_root\BolenFrontiers22.bib"},
  @{n='psutex.cls'; ol="$ol_root\psutex.cls"; pf="$pf_root\psutex.cls"}
)

$allMatch = $true
foreach ($p in $pairs) {
  $a = Normalize $p.pf
  $b = Normalize $p.ol
  if ($a -ceq $b) {
    Write-Output ('MATCH(normalized)  ' + $p.n)
  } else {
    $allMatch = $false
    Write-Output ('*** MISMATCH  ' + $p.n + ' ***')
    # show first differing line
    $al = $a -split "`n"; $bl = $b -split "`n"
    $n = [Math]::Min($al.Count, $bl.Count)
    for ($i = 0; $i -lt $n; $i++) {
      if ($al[$i] -cne $bl[$i]) {
        Write-Output ('  first diff at normalized line ' + ($i+1))
        Write-Output ('  PF: ' + $al[$i].Substring(0, [Math]::Min(160, $al[$i].Length)))
        Write-Output ('  OL: ' + $bl[$i].Substring(0, [Math]::Min(160, $bl[$i].Length)))
        break
      }
    }
    if ($al.Count -ne $bl.Count) { Write-Output ('  line counts differ: PF ' + $al.Count + ' vs OL ' + $bl.Count) }
  }
}
Write-Output ('ALL_MATCH=' + $allMatch)
