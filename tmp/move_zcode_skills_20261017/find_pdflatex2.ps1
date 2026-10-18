$roots = @('C:\Program Files', 'C:\Users\Ben Bolen\AppData\Local\Programs', 'D:\Program Files', 'D:\')
$found = @()
foreach ($r in $roots) {
  if (-not (Test-Path $r)) { continue }
  $hits = Get-ChildItem -LiteralPath $r -Filter 'pdflatex.exe' -Recurse -Depth 4 -ErrorAction SilentlyContinue | Select-Object -First 3
  foreach ($h in $hits) { $found += $h.FullName; Write-Output ("HIT " + $h.FullName) }
}
if ($found.Count -eq 0) { Write-Output 'NO-PDFLATEX-FOUND' }
