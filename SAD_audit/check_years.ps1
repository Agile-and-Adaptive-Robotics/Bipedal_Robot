$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
function Read-CsvUtf8([string]$path) { return (Get-Content -LiteralPath $path -Encoding UTF8 | ConvertFrom-Csv) }

$bf  = Read-CsvUtf8 "$out\airtable_doi_backfill.csv"
$bio = Read-CsvUtf8 "$out\group_Biology_items_slim.csv"
$grp = Read-CsvUtf8 "$out\group_items_slim.csv"
$byKey = @{}
foreach ($b in $bio) { $byKey[$b.key] = $b }
foreach ($g in $grp) { if (-not $byKey.ContainsKey($g.key)) { $byKey[$g.key] = $g } }

function Get-Year([string]$d) {
  if (-not $d) { return $null }
  $m = [regex]::Match($d, '(19|20)\d{2}')
  if ($m.Success) { return [int]$m.Value }
  return $null
}

$mismatch = 0
foreach ($r in $bf) {
  if (-not $r.zotero_key -or -not $byKey.ContainsKey($r.zotero_key)) { continue }
  $zy = Get-Year $byKey[$r.zotero_key].date
  $ay = $null
  if ($r.title) {}
  # airtable year lives in the slim csv, join by id
}
$at = Read-CsvUtf8 "$out\airtable_papers_slim.csv"
$atYear = @{}
foreach ($a in $at) { if ($a.year) { $atYear[$a.id] = [int]$a.year } }
foreach ($r in $bf) {
  if (-not $r.zotero_key -or -not $byKey.ContainsKey($r.zotero_key)) { continue }
  if (-not $atYear.ContainsKey($r.id)) { continue }
  $zy = Get-Year $byKey[$r.zotero_key].date
  if ($null -ne $zy -and $zy -ne $atYear[$r.id]) {
    $mismatch++
    Write-Output ("YEAR  " + $r.id + "  AT=" + $atYear[$r.id] + "  Zotero=" + $zy + "  zdate='" + $byKey[$r.zotero_key].date + "'  :: " + $r.title.Substring(0, [Math]::Min(70, $r.title.Length)))
  }
}
Write-Output ("mismatches: " + $mismatch)
Write-Output "DONE"
