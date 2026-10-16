$diffDir = 'D:\Github\Bipedal_Robot\tmp\disscampaign_20261015\diffs'
foreach ($df in Get-ChildItem $diffDir -Filter '*.diff' | Sort-Object Name) {
  $lines = Get-Content $df.FullName | Where-Object { $_ -match '^[+-]' -and $_ -notmatch '^(\+\+\+|---)' }
  if (-not $lines) { continue }
  $residual = @()
  $prevMinus = $null
  foreach ($l in $lines) {
    $content = $l.Substring(1)
    if ($content -match '^\\graphicspath') { $residual += '[graphicspath] ' + $content; continue }
    if ($content.Trim() -eq '') { $residual += '[blank/whitespace-only line]'; continue }
    $residual += '[SUBSTANTIVE] ' + $content.Substring(0, [Math]::Min(90, $content.Length))
  }
  Write-Output ('=== ' + $df.Name + ' ===')
  $residual | ForEach-Object { Write-Output ('  ' + $_) }
}
Write-Output 'AUDIT_DONE'
