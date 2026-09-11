$ErrorActionPreference = "Stop"
$base = "http://localhost:23119/api"
$out  = "D:\Github\Bipedal_Robot\SAD_audit"
New-Item -ItemType Directory -Force -Path $out | Out-Null

function Get-AllPages([string]$url) {
  $all = @(); $start = 0
  while ($true) {
    $sep = if ($url.Contains("?")) { "&" } else { "?" }
    $page = $null
    for ($try = 1; $try -le 3 -and -not $page; $try++) {
      $json = & curl.exe -s -m 60 ("$url$sep" + "limit=100&start=$start")
      $page = @(($json -join "`n" | ConvertFrom-Json) | Where-Object { $_ })
    }
    if (-not $page -or $page.Count -eq 0) { break }
    $all += $page
    if ($page.Count -lt 100) { break }
    $start += 100
    if ($start -gt 20000) { throw "pagination guard hit for $url" }
  }
  return $all
}

# --- Personal library: all collections, locate SAD + descendants
$colls = @(Get-AllPages "$base/users/0/collections")
$colls | ConvertTo-Json -Depth 10 | Out-File "$out\personal_all_collections.json" -Encoding utf8
$sad = $colls | Where-Object { $_.data.name -eq "Sensory Afferent Database" }
if (-not $sad) { throw "Sensory Afferent Database collection not found in personal library" }
Write-Output ("SAD collection key: " + $sad.key)

$keys = @($sad.key)
do {
  $new = @($colls | Where-Object { $keys -contains $_.data.parentCollection -and $keys -notcontains $_.key } | ForEach-Object { $_.key })
  $keys += $new
} while ($new.Count -gt 0)
Write-Output ("SAD collection tree size: " + $keys.Count)

$allItems = @{}
foreach ($k in $keys) {
  $it = @(Get-AllPages "$base/users/0/collections/$k/items/top")
  foreach ($i in $it) { if ($i -and $i.key) { $allItems[$i.key] = $i } }
}
Write-Output ("Personal SAD unique top-level items: " + $allItems.Count)
$allItems.Values | ConvertTo-Json -Depth 10 | Out-File "$out\personal_SAD_items_full.json" -Encoding utf8
$slim = $allItems.Values | ForEach-Object {
  [pscustomobject]@{ key = $_.key; type = $_.data.itemType; title = $_.data.title; DOI = $_.data.DOI; date = $_.data.date; pub = $_.data.publicationTitle }
}
$slim | ConvertTo-Json -Depth 5 | Out-File "$out\personal_SAD_items_slim.json" -Encoding utf8
$slim | Export-Csv "$out\personal_SAD_items_slim.csv" -NoTypeInformation -Encoding utf8

# --- AARL group library 735051
try {
  $gc = @(Get-AllPages "$base/groups/735051/collections")
  $gc | ConvertTo-Json -Depth 10 | Out-File "$out\group_collections.json" -Encoding utf8
  Write-Output ("AARL group collections: " + $gc.Count)

  $gItems = @{}
  foreach ($c in $gc) {
    if (-not $c.key) { continue }
    $it = @(Get-AllPages ("$base/groups/735051/collections/" + $c.key + "/items/top"))
    foreach ($i in $it) { if ($i -and $i.key) { $gItems[$i.key] = $i } }
  }
  # plus items not in any collection
  $all = @(Get-AllPages "$base/groups/735051/items/top")
  foreach ($i in $all) { if ($i -and $i.key) { $gItems[$i.key] = $i } }
  Write-Output ("AARL group unique top-level items: " + $gItems.Count)
  $gItems.Values | ConvertTo-Json -Depth 10 | Out-File "$out\group_items_full.json" -Encoding utf8
  $gSlim = $gItems.Values | ForEach-Object {
    [pscustomobject]@{ key = $_.key; type = $_.data.itemType; title = $_.data.title; DOI = $_.data.DOI; date = $_.data.date; pub = $_.data.publicationTitle }
  }
  $gSlim | Export-Csv "$out\group_items_slim.csv" -NoTypeInformation -Encoding utf8
} catch {
  Write-Output ("GROUP FETCH FAILED: " + $_.Exception.Message)
}
Write-Output "DONE"
