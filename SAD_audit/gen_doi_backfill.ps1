$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"

function Norm-Strict([string]$s) {
  if ([string]::IsNullOrWhiteSpace($s)) { return "" }
  return ([regex]::Replace($s, '\s+', ' ')).Trim().ToLower()
}
function Norm-Loose([string]$s) {
  if ([string]::IsNullOrWhiteSpace($s)) { return "" }
  $t = ([regex]::Replace($s, '\s+', ' ')).Trim().ToLower()
  $t = [regex]::Replace($t, '[^a-z0-9]+', ' ')
  return ([regex]::Replace($t, '\s+', ' ')).Trim()
}
function Read-CsvUtf8([string]$path) {
  return (Get-Content -LiteralPath $path -Encoding UTF8 | ConvertFrom-Csv)
}

$at  = Read-CsvUtf8 "$out\airtable_papers_slim.csv"
$bio = Read-CsvUtf8 "$out\group_Biology_items_slim.csv"
$grp = Read-CsvUtf8 "$out\group_items_slim.csv"

# Biology lookups: normalized title -> LIST of rows (to detect ambiguous collisions)
$strictAll = @{}; $looseAll = @{}
foreach ($b in $bio) {
  $ks = Norm-Strict $b.title; if ($ks) { if (-not $strictAll.ContainsKey($ks)) { $strictAll[$ks] = New-Object System.Collections.ArrayList }; [void]$strictAll[$ks].Add($b) }
  $kl = Norm-Loose  $b.title; if ($kl) { if (-not $looseAll.ContainsKey($kl))  { $looseAll[$kl]  = New-Object System.Collections.ArrayList }; [void]$looseAll[$kl].Add($b) }
}
$grpStrict = @{}; $grpLoose = @{}
foreach ($g in $grp) {
  $ks = Norm-Strict $g.title; if ($ks -and -not $grpStrict.ContainsKey($ks)) { $grpStrict[$ks] = $g }
  $kl = Norm-Loose  $g.title; if ($kl -and -not $grpLoose.ContainsKey($kl))  { $grpLoose[$kl]  = $g }
}

# Manual adjudications from the audit (airtable id -> zotero key), incl. Geyer via group
$adj = @{
  'rec8WdUONYgfvw1zw' = 'HBAEC9DM'; 'recHUhzNl3CEaF8pE' = 'HJUK3CTI';
  'recPoyqO9O9txmwm5' = 'PUVMGVWH'; 'recUQzJIupLCDV5g1' = 'HQZAC8XN';
  'recfiz9tzgn955QdA' = 'FIN3YLAZ'; 'reclwbsReGhWM71WF' = 'BGU8UZJE';
  'recxGoi3bAS1glsgD' = 'JA7MEXI6'; 'recDT62aen4FMrVhd' = 'INYJVQYV'
}
$byKey = @{}
foreach ($b in $bio) { $byKey[$b.key] = $b }
foreach ($g in $grp) { if (-not $byKey.ContainsKey($g.key)) { $byKey[$g.key] = $g } }

function Pick-Doi($rows) {
  # prefer a row that actually has a DOI
  foreach ($r in $rows) { if ($r.DOI -and $r.DOI.Trim()) { return $r } }
  return $rows[0]
}

$res = New-Object System.Collections.ArrayList
$ambiguous = 0
foreach ($a in $at) {
  $row = [pscustomobject]@{ id = $a.id; title = $a.title; match = ''; zotero_key = ''; DOI = ''; flag = '' }
  if ($adj.ContainsKey($a.id)) {
    $z = $byKey[$adj[$a.id]]
    $row.match = 'adjudicated'; $row.zotero_key = $z.key; $row.DOI = $z.DOI
  } else {
    $ks = Norm-Strict $a.title
    if ($ks -and $strictAll.ContainsKey($ks)) {
      $rows = @($strictAll[$ks])
      $dois = @($rows | Where-Object { $_.DOI -and $_.DOI.Trim() } | ForEach-Object { $_.DOI.Trim().ToLower() } | Select-Object -Unique)
      $z = Pick-Doi $rows
      $row.match = 'strict'; $row.zotero_key = $z.key; $row.DOI = $z.DOI
      if ($dois.Count -gt 1) { $row.flag = 'AMBIGUOUS: multiple DOIs for same normalized title'; $ambiguous++ }
    } else {
      $kl = Norm-Loose $a.title
      if ($kl -and $looseAll.ContainsKey($kl)) {
        $rows = @($looseAll[$kl])
        $dois = @($rows | Where-Object { $_.DOI -and $_.DOI.Trim() } | ForEach-Object { $_.DOI.Trim().ToLower() } | Select-Object -Unique)
        $z = Pick-Doi $rows
        $row.match = 'loose'; $row.zotero_key = $z.key; $row.DOI = $z.DOI
        if ($dois.Count -gt 1) { $row.flag = 'AMBIGUOUS: multiple DOIs for same normalized title'; $ambiguous++ }
      } else {
        $row.match = 'none'
      }
    }
  }
  [void]$res.Add($row)
}
$res | Export-Csv "$out\airtable_doi_backfill.csv" -NoTypeInformation -Encoding utf8

$have = @($res | Where-Object { $_.DOI -and $_.DOI.Trim() })
Write-Output ("records with DOI to write: " + $have.Count + " / " + @($at).Count + "   ambiguous flagged: " + $ambiguous)
Write-Output "--- non-strict matches (review these pairs) ---"
foreach ($r in $res) {
  if ($r.match -ne 'strict' -and $r.match -ne 'none') {
    $ztitle = ''
    if ($row) {}
    if ($r.zotero_key -and $byKey.ContainsKey($r.zotero_key)) { $ztitle = $byKey[$r.zotero_key].title }
    Write-Output ("  [" + $r.match + "] " + $r.id + " " + $r.flag)
    Write-Output ("     AT : " + $r.title)
    Write-Output ("     ZOT: " + $ztitle + "  DOI=" + $r.DOI)
  }
}
Write-Output "--- no match (leave blank) ---"
foreach ($r in $res) { if ($r.match -eq 'none') { Write-Output ("  " + $r.id + " :: " + $r.title) } }
Write-Output "DONE"
