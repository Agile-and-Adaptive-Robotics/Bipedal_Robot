$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"

$gap = Get-Content "$out\personalSAD_not_in_airtable_corrected.csv" -Encoding UTF8 | ConvertFrom-Csv
Write-Output ("gap list rows: " + @($gap).Count)

# full personal items for creators
$full = Get-Content "$out\personal_SAD_items_full.json" -Raw -Encoding UTF8 | ConvertFrom-Json
$full = @($full)
$byKey = @{}
foreach ($i in $full) { if ($i -and $i.key) { $byKey[$i.key] = $i } }
Write-Output ("personal full items: " + $full.Count)

function Get-Year([string]$d) {
  if (-not $d) { return $null }
  $m = [regex]::Match($d, '(19|20)\d{2}')
  if ($m.Success) { return [int]$m.Value }
  return $null
}
function Get-Author($item) {
  $auths = @($item.data.creators | Where-Object { $_.creatorType -eq 'author' })
  # fetch double-recorded creators (single "name" variant + firstName/lastName variant) -> dedupe by surname
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
$n = 0
foreach ($g in $gap) {
  if ($n -ge 50) { break }
  $it = $byKey[$g.key]
  if (-not $it) { continue }
  $title = $it.data.title
  if ([string]::IsNullOrWhiteSpace($title)) { continue }
  $doi = $it.data.DOI; if (-not $doi) { $doi = '' }
  [void]$rows.Add([pscustomobject]@{
    zotero_key = $g.key; title = $title; author = (Get-Author $it)
    year = (Get-Year $it.data.date); DOI = $doi.Trim()
  })
  $n++
}
$rows | Export-Csv "$out\airtable_create50.csv" -NoTypeInformation -Encoding utf8
Write-Output ("prepared create batch: " + $rows.Count + " records")
foreach ($r in $rows) {
  Write-Output ($r.zotero_key + " | " + $r.author + " | " + $r.year + " | " + $r.DOI + " | " + $r.title.Substring(0, [Math]::Min(80, $r.title.Length)))
}
Write-Output "DONE"
