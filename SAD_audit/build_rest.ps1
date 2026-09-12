$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$gap = Get-Content "$out\personalSAD_not_in_airtable_corrected.csv" -Encoding UTF8 | ConvertFrom-Csv
$used = Get-Content "$out\airtable_created50_ids.csv" -Encoding UTF8 | ConvertFrom-Csv
$usedKeys = @($used | ForEach-Object { $_.zotero_key })
$skip = @('C6252AEY', 'QDN8VDSK')
$full = Get-Content "$out\personal_SAD_items_full.json" -Raw -Encoding UTF8 | ConvertFrom-Json
$byKey = @{}; foreach ($i in @($full)) { if ($i -and $i.key) { $byKey[$i.key] = $i } }
function Get-Year([string]$d) { if (-not $d) { return $null }; $m = [regex]::Match($d, '(19|20)\d{2}'); if ($m.Success) { return [int]$m.Value }; return $null }
function Get-Author($item) {
  $auths = @($item.data.creators | Where-Object { $_.creatorType -eq 'author' })
  $seen = New-Object System.Collections.Generic.HashSet[string]
  $names = New-Object System.Collections.ArrayList
  foreach ($a in $auths) {
    $ln = $a.lastName
    if (-not $ln -and $a.name) { $parts = $a.name.Trim() -split '\s+'; if ($parts.Count -gt 0) { $ln = $parts[-1] } }
    if (-not $ln) { continue }
    $k = $ln.Trim().ToLower()
    if (-not $seen.Contains($k)) { [void]$seen.Add($k); [void]$names.Add($ln.Trim()) }
  }
  if ($names.Count -eq 0) { return '' }
  if ($names.Count -eq 1) { return $names[0] }
  if ($names.Count -eq 2) { return ($names[0] + ' and ' + $names[1]) }
  return ($names[0] + ' et al.')
}
$rows = New-Object System.Collections.ArrayList
foreach ($g in $gap) {
  if ($usedKeys -contains $g.key -or $skip -contains $g.key) { continue }
  $it = $byKey[$g.key]; if (-not $it) { continue }
  $title = $it.data.title
  if ([string]::IsNullOrWhiteSpace($title)) { continue }
  $doi = $it.data.DOI; if (-not $doi) { $doi = '' }
  [void]$rows.Add([pscustomobject]@{ zotero_key = $g.key; title = $title; author = (Get-Author $it); year = (Get-Year $it.data.date); DOI = $doi.Trim() })
}
Write-Output ("remaining to import: " + $rows.Count)
$rows | Export-Csv "$out\airtable_rest_import.csv" -NoTypeInformation -Encoding utf8
# dump for transcription, 50 per chunk marker
$n = 0
foreach ($r in $rows) {
  if ($n % 50 -eq 0) { Write-Output ("===== CHUNK " + ([int]($n / 50) + 1) + " =====") }
  $y = 'null'; if ($r.year) { $y = $r.year }
  Write-Output ($r.zotero_key + ' ||| ' + $r.title + ' ||| ' + $r.author + ' ||| ' + $y + ' ||| ' + $r.DOI)
  $n++
}
Write-Output "DONE"
