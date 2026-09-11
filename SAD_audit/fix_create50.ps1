$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$rows = Get-Content "$out\airtable_create50.csv" -Encoding UTF8 | ConvertFrom-Csv
$rows = @($rows | Where-Object { $_.zotero_key -ne 'C6252AEY' -and $_.zotero_key -ne 'QDN8VDSK' })
foreach ($r in $rows) { if ($r.zotero_key -eq 'MCU88JT6') { $r.author = 'Shik et al.' } }

$gap = Get-Content "$out\personalSAD_not_in_airtable_corrected.csv" -Encoding UTF8 | ConvertFrom-Csv
$used = @($rows | ForEach-Object { $_.zotero_key })
$full = Get-Content "$out\personal_SAD_items_full.json" -Raw -Encoding UTF8 | ConvertFrom-Json
$byKey = @{}
foreach ($i in @($full)) { if ($i -and $i.key) { $byKey[$i.key] = $i } }
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
$need = 50 - $rows.Count
$added = 0
$skip = @('C6252AEY', 'QDN8VDSK')  # mangled Zotero records - curate manually, never import
foreach ($g in $gap) {
  if ($added -ge $need) { break }
  if ($skip -contains $g.key -or $used -contains $g.key) { continue }
  $it = $byKey[$g.key]
  if (-not $it) { continue }
  $title = $it.data.title
  if ([string]::IsNullOrWhiteSpace($title)) { continue }
  # skip obviously mangled titles (PDF-filename style, fused words) for manual curation later
  if ($title -match '^[A-Z][a-z]+ ?[A-Z][a-z]+ ?- ') { continue }
  $doi = $it.data.DOI; if (-not $doi) { $doi = '' }
  $rows += [pscustomobject]@{ zotero_key = $g.key; title = $title; author = (Get-Author $it); year = (Get-Year $it.data.date); DOI = $doi.Trim() }
  $added++
}
$rows = @($rows | Select-Object -First 50)
$rows | Export-Csv "$out\airtable_create50.csv" -NoTypeInformation -Encoding utf8
Write-Output ("final batch: " + $rows.Count + " records")
foreach ($r in $rows | Select-Object -Last 4) { Write-Output ($r.zotero_key + " | " + $r.author + " | " + $r.year + " | " + $r.DOI + " | " + $r.title.Substring(0, [Math]::Min(80, $r.title.Length))) }
Write-Output "DONE"
