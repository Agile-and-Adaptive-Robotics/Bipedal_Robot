# Style/mechanical scan of ProofFinal chapters — read-only, outputs file:line hits.
$ch = 'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters'
$main = 'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\main.tex'
$files = Get-ChildItem "$ch\*.tex" | Sort-Object Name
$files += Get-Item $main

$patterns = [ordered]@{
  'campaign'      = 'campaign'
  'of-record'     = 'of record'
  'year-in-prose' = '(20(1[5-9]|2[0-6])[-_][0-1][0-9])'   # 2015-2026 with -MM (dates like 2026-09)
  'phi-macro'     = '\\phi\b'
  'diameter-mac'  = '\\diameter'
  'plant-word'    = '\bplant\b'
  'em-dash'       = '---'
  'pool-word'     = '\bpool\b'
  'SNN'           = '\bSNN\b'
  'DoggyDeux'     = 'DoggyDeux'
  'MuscleMutt'    = 'Muscle Mutt'
  'point1ms'      = '0\.1\s*\\?(qty|si)?\s*\\?ms'
  'neuron406'     = '406'
  'realtime-claim'= 'faster than real time'
  'bare-cite'     = '\\cite\{'
  'fit-set'       = 'fit set'
  'roadmap-fig'   = 'fig:roadmap'
  'pipeline-fig'  = 'fig:pipeline'
  'harness'       = 'harness'
  '42cm-flexor'   = '42\.0'
  '41.5-flexor'   = '41\.5'
  'back-calc'     = 'back-calculat'
}

foreach ($pat in $patterns.Keys) {
  Write-Output ("=== " + $pat + " ===")
  $hits = 0
  foreach ($f in $files) {
    $ln = 0
    foreach ($line in [IO.File]::ReadLines($f.FullName)) {
      $ln++
      if ($line -match $patterns[$pat]) {
        $hits++
        $short = if ($line.Length -gt 160) { $line.Substring(0,160) + '...' } else { $line }
        Write-Output ("  {0}:{1}: {2}" -f $f.Name, $ln, $short.Trim())
        if ($hits -ge 25) { Write-Output '  ...capped at 25'; break }
      }
    }
    if ($hits -ge 25) { break }
  }
  if ($hits -eq 0) { Write-Output '  (no hits)' }
}
