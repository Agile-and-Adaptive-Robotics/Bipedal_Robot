$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"

$colls = Get-Content "$out\group_collections.json" -Raw | ConvertFrom-Json
$colls = @($colls)
Write-Output ("collections parsed: " + $colls.Count)

# key -> collection map
$cmap = @{}
foreach ($c in $colls) { if ($c -and $c.key) { $cmap[$c.key] = $c } }

# reference CSV of the whole tree
$rows = foreach ($c in $colls) { if ($c -and $c.key) { [pscustomobject]@{ key = $c.key; name = $c.data.name; parent = $c.data.parentCollection } } }
$rows | Export-Csv "$out\group_collections_tree.csv" -NoTypeInformation -Encoding utf8

# find collection named exactly 'Biology' (case-insensitive)
$bioKeys = @()
foreach ($c in $colls) {
  if ($c -and $c.data -and $c.data.name -and $c.data.name.Trim().ToLower() -eq 'biology') { $bioKeys += $c.key }
}
if ($bioKeys.Count -eq 0) {
  Write-Output "No collection named exactly 'Biology'; candidates containing 'biolog':"
  foreach ($c in $colls) { if ($c -and $c.data -and $c.data.name -and $c.data.name.ToLower().Contains('biolog')) { Write-Output ("  '" + $c.data.name + "' (" + $c.key + ")") } }
  exit 1
}
Write-Output ("Biology root key: " + ($bioKeys -join ","))

# descendants
$keySet = @{}
foreach ($k in $bioKeys) { $keySet[$k] = 1 }
$changed = $true
while ($changed) {
  $changed = $false
  foreach ($c in $colls) {
    if ($c -and $c.key -and -not $keySet.ContainsKey($c.key)) {
      $par = $c.data.parentCollection
      if ($par -is [string] -and $par.Length -gt 0 -and $keySet.ContainsKey($par)) {
        $keySet[$c.key] = 1; $changed = $true
      }
    }
  }
}
Write-Output ("Biology subtree collection count: " + $keySet.Count)

# items
$items = Get-Content "$out\group_items_full.json" -Raw | ConvertFrom-Json
$items = @($items)
Write-Output ("items parsed: " + $items.Count)

# per-collection counts + subtree items
$collCount = @{}
$inBio = New-Object System.Collections.ArrayList
foreach ($i in $items) {
  if (-not $i -or -not $i.key) { continue }
  $cks = @($i.data.collections)
  $isBio = $false
  foreach ($ck in $cks) {
    if ($collCount.ContainsKey($ck)) { $collCount[$ck]++ } else { $collCount[$ck] = 1 }
    if ($keySet.ContainsKey($ck)) { $isBio = $true }
  }
  if ($isBio) { [void]$inBio.Add($i) }
}
Write-Output ("Unique items in Biology subtree: " + $inBio.Count)

# tree listing with counts, depth-indented, only for subtree
$tree = New-Object System.Collections.ArrayList
foreach ($c in $colls) {
  if ($c -and $c.key -and $keySet.ContainsKey($c.key)) {
    $depth = 0; $cur = $c
    while ($cur.data.parentCollection -is [string] -and $cur.data.parentCollection.Length -gt 0) {
      $depth++
      if (-not $cmap.ContainsKey($cur.data.parentCollection)) { break }
      $cur = $cmap[$cur.data.parentCollection]
      if ($depth -gt 20) { break }
    }
    $n = 0; if ($collCount.ContainsKey($c.key)) { $n = $collCount[$c.key] }
    [void]$tree.Add([pscustomobject]@{ depth = $depth; name = $c.data.name; key = $c.key; items = $n })
  }
}
$tree = $tree | Sort-Object depth, name
$tree | Export-Csv "$out\group_Biology_tree.csv" -NoTypeInformation -Encoding utf8

# slim CSV of subtree items
$slim = New-Object System.Collections.ArrayList
foreach ($i in $inBio) {
  [void]$slim.Add([pscustomobject]@{
    key = $i.key; type = $i.data.itemType; title = $i.data.title; DOI = $i.data.DOI
    date = $i.data.date; pub = $i.data.publicationTitle; collections = (@($i.data.collections) -join ";")
  })
}
$slim | Export-Csv "$out\group_Biology_items_slim.csv" -NoTypeInformation -Encoding utf8

# personal SAD vs Biology subtree
$p = Import-Csv "$out\personal_SAD_items_slim.csv"
$bd = @{}; $bt = @{}
foreach ($s in $slim) {
  if ($s.DOI) { $bd[$s.DOI.Trim().ToLower()] = 1 }
  if ($s.title) { $bt[$s.title.Trim().ToLower()] = 1 }
}
$miss = New-Object System.Collections.ArrayList
$hit = 0
foreach ($x in $p) {
  $found = $false
  if ($x.DOI -and $x.DOI.Trim() -and $bd.ContainsKey($x.DOI.Trim().ToLower())) { $found = $true }
  if (-not $found -and $x.title -and $x.title.Trim() -and $bt.ContainsKey($x.title.Trim().ToLower())) { $found = $true }
  if ($found) { $hit++ } else { [void]$miss.Add($x) }
}
Write-Output ("Personal SAD in Biology subtree: " + $hit + " / " + $p.Count)
Write-Output ("Personal SAD NOT in Biology subtree: " + $miss.Count)
$miss | Select-Object key, type, title, DOI, date | Export-Csv "$out\preview_personal_not_in_Biology.csv" -NoTypeInformation -Encoding utf8
Write-Output "DONE"
