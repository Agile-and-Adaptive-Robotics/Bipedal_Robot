$ErrorActionPreference = "Continue"
# 1. Has the desktop synced the web-API deletions?
curl.exe -s -o "$env:TEMP\z_local1.json" -w "local 375JYEPB: %{http_code}`n" "http://localhost:23119/api/users/0/items/375JYEPB"
curl.exe -s -o "$env:TEMP\z_local2.json" -w "local QJTQVEJR: %{http_code}`n" "http://localhost:23119/api/users/0/items/QJTQVEJR"

# 2. Notes still present in personal library? (count via local API)
$n = 0; $start = 0
while ($true) {
  curl.exe -s -o "$env:TEMP\z_notes.json" "http://localhost:23119/api/users/0/items?format=json&itemType=note&limit=100&start=$start"
  $r = Get-Content "$env:TEMP\z_notes.json" -Raw -Encoding UTF8
  $items = $r | ConvertFrom-Json
  $c = @($items).Count
  $n += $c
  if ($c -lt 100) { break }
  $start += 100
}
Write-Output ("personal notes (local): " + $n)

# 3. What exactly was deleted: content-type/filename summary
$elig = Get-Content "D:\Github\Bipedal_Robot\SAD_audit\personal_SAD_attachment_purge_eligible.csv" -Encoding UTF8 | ConvertFrom-Csv
Write-Output ("purged attachments: " + @($elig).Count)
$elig | Group-Object linkMode | ForEach-Object { Write-Output ("  " + $_.Name + ": " + $_.Count) }
$noPdf = @($elig | Where-Object { -not ($_.filename -like '*.pdf') })
Write-Output ("  non-PDF filenames: " + $noPdf.Count)
$noPdf | Select-Object -First 10 | ForEach-Object { Write-Output ("    " + $_.attKey + " " + $_.filename + " (parent " + $_.parentKey + ")") }

# 4. Do the AARL group copies have FILE attachments? (local API, group 735051)
$gatt = @{}; $start = 0
while ($true) {
  curl.exe -s -o "$env:TEMP\z_gatt.json" "http://localhost:23119/api/groups/735051/items?format=json&itemType=attachment&limit=100&start=$start"
  $r = Get-Content "$env:TEMP\z_gatt.json" -Raw -Encoding UTF8
  $items = $r | ConvertFrom-Json
  $c = @($items).Count
  if ($c -eq 0) { break }
  foreach ($it in @($items)) {
    if ($it.data.linkMode -eq 'imported_file' -or $it.data.linkMode -eq 'imported_url') {
      if (-not $gatt.ContainsKey($it.data.parentItem)) { $gatt[$it.data.parentItem] = 0 }
      $gatt[$it.data.parentItem]++
    }
  }
  if ($c -lt 100) { break }
  $start += 100
}
Write-Output ("AARL group items WITH file attachments: " + $gatt.Count)

# map purged personal parents -> group items by loose title
function Norm-Loose([string]$s) {
  if ([string]::IsNullOrWhiteSpace($s)) { return "" }
  $t = ([regex]::Replace($s, '\s+', ' ')).Trim().ToLower()
  $t = [regex]::Replace($t, '[^a-z0-9]+', ' ')
  return ([regex]::Replace($t, '\s+', ' ')).Trim()
}
$per = Get-Content "D:\Github\Bipedal_Robot\SAD_audit\personal_SAD_items_slim.csv" -Encoding UTF8 | ConvertFrom-Csv
$perBy = @{}; foreach ($p in $per) { $perBy[$p.key] = $p }
$bio = Get-Content "D:\Github\Bipedal_Robot\SAD_audit\group_Biology_items_slim.csv" -Encoding UTF8 | ConvertFrom-Csv
$bioLoose = @{}; foreach ($b in $bio) { $lt = Norm-Loose $b.title; if ($lt -and -not $bioLoose.ContainsKey($lt)) { $bioLoose[$lt] = $b.key } }
$parents = @($elig | Select-Object -ExpandProperty parentKey -Unique)
$covered = 0; $uncovered = New-Object System.Collections.ArrayList
foreach ($pk in $parents) {
  $p = $perBy[$pk]; if (-not $p) { continue }
  $lt = Norm-Loose $p.title
  $gk = ''; if ($lt -and $bioLoose.ContainsKey($lt)) { $gk = $bioLoose[$lt] }
  if ($gk -ne '' -and $gatt.ContainsKey($gk)) { $covered++ } else { [void]$uncovered.Add($p.title) }
}
Write-Output ("purged-parent items whose AARL copy HAS a file attachment: " + $covered + " / " + $parents.Count)
Write-Output ("titles with NO file anywhere in AARL (first 15):")
$uncovered | Select-Object -First 15 | ForEach-Object { Write-Output ("  - " + $_.Substring(0, [Math]::Min(85, $_.Length))) }
