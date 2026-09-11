$ErrorActionPreference = "Continue"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
function Norm-Loose([string]$s) {
  if ([string]::IsNullOrWhiteSpace($s)) { return "" }
  $t = ([regex]::Replace($s, '\s+', ' ')).Trim().ToLower()
  $t = [regex]::Replace($t, '[^a-z0-9]+', ' ')
  return ([regex]::Replace($t, '\s+', ' ')).Trim()
}

# Airtable records WITH file attachments
$atWithFile = New-Object System.Collections.Generic.HashSet[string]
foreach ($id in (Get-Content "$out\airtable_with_attachment_ids.csv" -Encoding UTF8)) { $t = $id.Trim(); if ($t) { [void]$atWithFile.Add($t) } }
Write-Output ("airtable records with files: " + $atWithFile.Count)

# AT id -> DOI (patched slim)
$slim = Get-Content "$out\airtable_papers_slim.csv" -Encoding UTF8 | ConvertFrom-Csv
$atDoiById = @{}; $atIdByLoose = @{}
foreach ($s in $slim) {
  if ($s.DOI) { $atDoiById[$s.id] = $s.DOI.Trim().ToLower() }
  $lt = Norm-Loose $s.title; if ($lt) { $atIdByLoose[$lt] = $s.id }
}

# AARL group: parents with file attachments (guard null parents = standalone attachments)
$gatt = New-Object System.Collections.Generic.HashSet[string]
$start = 0; $tmpJson = "$env:TEMP\z_gatt2.json"
while ($true) {
  curl.exe -s -o $tmpJson "http://localhost:23119/api/groups/73551placeholder" # placeholder never used
  break
}
# real loop
$start = 0
while ($true) {
  curl.exe -s -o $tmpJson "http://localhost:23119/api/groups/735051/items?format=json&itemType=attachment&limit=100&start=$start"
  $r = Get-Content $tmpJson -Raw -Encoding UTF8
  $items = $r | ConvertFrom-Json
  $c = @($items).Count
  if ($c -eq 0) { break }
  foreach ($it in @($items)) {
    $par = $it.data.parentItem
    if ($par -and ($it.data.linkMode -eq 'imported_file' -or $it.data.linkMode -eq 'imported_url')) { [void]$gatt.Add($par) }
  }
  if ($c -lt 100) { break }
  $start += 100
}
Write-Output ("AARL group items with file attachments: " + $gatt.Count)

# group Biology key by loose title
$bio = Get-Content "$out\group_Biology_items_slim.csv" -Encoding UTF8 | ConvertFrom-Csv
$bioKeyByLoose = @{}
foreach ($b in $bio) { $lt = Norm-Loose $b.title; if ($lt -and -not $bioKeyByLoose.ContainsKey($lt)) { $bioKeyByLoose[$lt] = $b.key } }

$elig = Get-Content "$out\personal_SAD_attachment_purge_eligible.csv" -Encoding UTF8 | ConvertFrom-Csv
$per = Get-Content "$out\personal_SAD_items_slim.csv" -Encoding UTF8 | ConvertFrom-Csv
$perBy = @{}; foreach ($p in $per) { $perBy[$p.key] = $p }

$rows = New-Object System.Collections.ArrayList
$okBoth = 0; $over = New-Object System.Collections.ArrayList
foreach ($pk in (@($elig | Select-Object -ExpandProperty parentKey -Unique))) {
  $p = $perBy[$pk]; if (-not $p) { continue }
  $lt = Norm-Loose $p.title
  $pd = ''; if ($p.DOI) { $pd = $p.DOI.Trim().ToLower() }
  $atId = ''
  foreach ($s in $slim) { if ($pd -ne '' -and $s.DOI -and $s.DOI.Trim().ToLower() -eq $pd) { $atId = $s.id; break } }
  if ($atId -eq '' -and $lt -and $atIdByLoose.ContainsKey($lt)) { $atId = $atIdByLoose[$lt] }
  $atFile = ($atId -ne '' -and $atWithFile.Contains($atId))
  $gk = ''; if ($lt -and $bioKeyByLoose.ContainsKey($lt)) { $gk = $bioKeyByLoose[$lt] }
  $aarlFile = ($gk -ne '' -and $gatt.Contains($gk))
  $verdict = 'OK (files in both)'
  if (-not ($atFile -and $aarlFile)) { $verdict = 'OVER-PURGED'; [void]$over.Add($p.title) } else { $okBoth++ }
  [void]$rows.Add([pscustomobject]@{ parentKey = $pk; title = $p.title; aarlFile = $aarlFile; airtableFile = $atFile; verdict = $verdict })
}
$rows | Export-Csv "$out\purge_audit_corrected.csv" -NoTypeInformation -Encoding utf8
Write-Output ("verdict: OK-both = " + $okBoth + " / 80 ; OVER-PURGED = " + $over.Count)
$over | ForEach-Object { Write-Output ("  OVER: " + $_) }
Write-Output "DONE"
