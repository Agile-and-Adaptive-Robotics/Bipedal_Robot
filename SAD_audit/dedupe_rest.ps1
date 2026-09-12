$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
function Norm-Loose([string]$s) {
  if ([string]::IsNullOrWhiteSpace($s)) { return "" }
  $t = ([regex]::Replace($s, '\s+', ' ')).Trim().ToLower()
  $t = [regex]::Replace($t, '[^a-z0-9]+', ' ')
  return ([regex]::Replace($t, '\s+', ' ')).Trim()
}
$slim = Get-Content "$out\airtable_papers_slim.csv" -Encoding UTF8 | ConvertFrom-Csv
$tset = @{}; $dset = @{}
foreach ($s in $slim) {
  $lt = Norm-Loose $s.title; if ($lt) { $tset[$lt] = 1 }
  if ($s.DOI) { $dset[$s.DOI.Trim().ToLower()] = 1 }
}
$rest = Get-Content "$out\airtable_rest_import.csv" -Encoding UTF8 | ConvertFrom-Csv
$keep = New-Object System.Collections.ArrayList
$dups = 0
foreach ($r in $rest) {
  $lt = Norm-Loose $r.title
  $d = ''; if ($r.DOI) { $d = $r.DOI.Trim().ToLower() }
  if (($d -and $dset.ContainsKey($d)) -or ($lt -and $tset.ContainsKey($lt))) {
    $dups++
    Write-Output ("DUP-DROP: " + $r.zotero_key + " :: " + $r.title.Substring(0, [Math]::Min(70, $r.title.Length)))
  } else {
    [void]$keep.Add($r)
    if ($lt) { $tset[$lt] = 1 }
    if ($d) { $dset[$d] = 1 }
  }
}
Write-Output ("kept: " + $keep.Count + "  dropped as existing dups: " + $dups)
$keep | Export-Csv "$out\airtable_rest_import_clean.csv" -NoTypeInformation -Encoding utf8
