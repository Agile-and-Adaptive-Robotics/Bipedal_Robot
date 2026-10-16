$ol_root = 'D:\Github\Bipedal_Robot\tmp\disscampaign_20261015\overleaf_pull'
$pf_root = 'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal'

$files = @(
  @{n='main.tex'; ol="$ol_root\main.tex"; pf="$pf_root\main.tex"},
  @{n='10-introduction.tex'; ol="$ol_root\chapters\10-introduction.tex"; pf="$pf_root\chapters\10-introduction.tex"},
  @{n='15-background.tex'; ol="$ol_root\chapters\15-background.tex"; pf="$pf_root\chapters\15-background.tex"},
  @{n='20-methods.tex'; ol="$ol_root\chapters\20-methods.tex"; pf="$pf_root\chapters\20-methods.tex"},
  @{n='30-results.tex'; ol="$ol_root\chapters\30-results.tex"; pf="$pf_root\chapters\30-results.tex"},
  @{n='40-discussion.tex'; ol="$ol_root\chapters\40-discussion.tex"; pf="$pf_root\chapters\40-discussion.tex"},
  @{n='50-futurework.tex'; ol="$ol_root\chapters\50-futurework.tex"; pf="$pf_root\chapters\50-futurework.tex"},
  @{n='60-conclusion.tex'; ol="$ol_root\chapters\60-conclusion.tex"; pf="$pf_root\chapters\60-conclusion.tex"},
  @{n='92-AppendixA.tex'; ol="$ol_root\chapters\92-AppendixA.tex"; pf="$pf_root\chapters\92-AppendixA.tex"},
  @{n='93-AppendixB.tex'; ol="$ol_root\chapters\93-AppendixB.tex"; pf="$pf_root\chapters\93-AppendixB.tex"},
  @{n='94-AppendixC.tex'; ol="$ol_root\chapters\94-AppendixC.tex"; pf="$pf_root\chapters\94-AppendixC.tex"}
)

foreach ($f in $files) {
  $ol = [System.IO.File]::ReadAllText($f.ol)
  $pf = [System.IO.File]::ReadAllText($f.pf)
  $olLF = ([regex]::Matches($ol, "`n")).Count
  $pfLF = ([regex]::Matches($pf, "`n")).Count
  $olCR = ([regex]::Matches($ol, "`r")).Count
  $pfCR = ([regex]::Matches($pf, "`r")).Count
  $same = if ($ol -ceq $pf) { 'IDENTICAL' } else { 'DIFFER' }
  Write-Output ("{0}: {1} | OL {2} LF / {3} CR | PF {4} LF / {5} CR" -f $f.n, $same, $olLF, $olCR, $pfLF, $pfCR)
}
