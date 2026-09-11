$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"

# Normalization, consistent with analyze_biology.ps1 (trim+lower) plus:
#   strict = collapse all whitespace (Airtable titles are line-wrapped pastes)
#   loose  = strict + strip everything non-alphanumeric (catches punctuation drift)
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

$at   = Read-CsvUtf8 "$out\airtable_papers_slim.csv"
$bio  = Read-CsvUtf8 "$out\group_Biology_items_slim.csv"
$per  = Read-CsvUtf8 "$out\personal_SAD_items_slim.csv"
$tree = Read-CsvUtf8 "$out\group_Biology_tree.csv"

Write-Output ("airtable papers:  " + @($at).Count)
Write-Output ("biology items:    " + @($bio).Count)
Write-Output ("personal SAD:     " + @($per).Count)

# ---- build lookups -----------------------------------------------------------
function Build-Lookups($rows) {
  $strict = @{}; $loose = @{}
  foreach ($r in $rows) {
    $ks = Norm-Strict $r.title
    if ($ks -and -not $strict.ContainsKey($ks)) { $strict[$ks] = $r }
    $kl = Norm-Loose $r.title
    if ($kl -and -not $loose.ContainsKey($kl)) { $loose[$kl] = $r }
  }
  return @{ strict = $strict; loose = $loose }
}
$bioL = Build-Lookups $bio
$perL = Build-Lookups $per
$atL  = Build-Lookups $at

function Find-Match($title, $lookups) {
  # returns 'strict' | 'loose' | $null
  $ks = Norm-Strict $title; $kl = Norm-Loose $title
  if ($ks -and $lookups.strict.ContainsKey($ks)) { return 'strict' }
  if ($kl -and $lookups.loose.ContainsKey($kl)) { return 'loose' }
  return $null
}

# ---- airtable -> Biology / personal ------------------------------------------
$atNotBio = New-Object System.Collections.ArrayList
$nStrict = 0; $nLoose = 0
foreach ($a in $at) {
  $m = Find-Match $a.title $bioL
  if ($m -eq 'strict') { $nStrict++ }
  elseif ($m -eq 'loose') { $nLoose++ }
  else { [void]$atNotBio.Add($a) }
}
Write-Output ""
Write-Output ("== Airtable -> AARL Biology ==")
Write-Output ("matched strict: " + $nStrict + "  loose: " + $nLoose + "  NOT FOUND: " + $atNotBio.Count + " / " + @($at).Count)
foreach ($a in $atNotBio) {
  Write-Output ("  MISS  [" + $a.id + "] " + $a.author + " " + $a.year + " :: " + $a.title)
}
$atNotBio | Select-Object id, title, author, year |
  Export-Csv "$out\airtable_not_in_Biology.csv" -NoTypeInformation -Encoding utf8

# ---- airtable -> personal SAD -------------------------------------------------
$nPerS = 0; $nPerL = 0; $atNotPer = New-Object System.Collections.ArrayList
foreach ($a in $at) {
  $m = Find-Match $a.title $perL
  if ($m -eq 'strict') { $nPerS++ }
  elseif ($m -eq 'loose') { $nPerL++ }
  else { [void]$atNotPer.Add($a) }
}
Write-Output ""
Write-Output ("== Airtable -> personal SAD (strict view) ==")
Write-Output ("matched strict: " + $nPerS + "  loose: " + $nPerL + "  not in personal SAD: " + $atNotPer.Count + " / " + @($at).Count)
$atNotPer | Select-Object id, title, author, year |
  Export-Csv "$out\airtable_not_in_personalSAD.csv" -NoTypeInformation -Encoding utf8

# ---- Biology -> airtable ------------------------------------------------------
$bioNotAt = New-Object System.Collections.ArrayList
foreach ($b in $bio) {
  $m = Find-Match $b.title $atL
  if (-not $m) { [void]$bioNotAt.Add($b) }
}
$bioHit = @($bio).Count - $bioNotAt.Count
Write-Output ""
Write-Output ("== AARL Biology -> Airtable ==")
Write-Output ("in Airtable: " + $bioHit + " / " + @($bio).Count + "   MISSING from Airtable: " + $bioNotAt.Count)
$bioNotAt | Select-Object key, type, title, DOI, date, pub, collections |
  Export-Csv "$out\Biology_not_in_airtable.csv" -NoTypeInformation -Encoding utf8

# ---- personal -> airtable ------------------------------------------------------
$perNotAt = New-Object System.Collections.ArrayList
foreach ($p in $per) {
  $m = Find-Match $p.title $atL
  if (-not $m) { [void]$perNotAt.Add($p) }
}
Write-Output ""
Write-Output ("== personal SAD -> Airtable ==")
Write-Output ("in Airtable: " + (@($per).Count - $perNotAt.Count) + " / " + @($per).Count + "   MISSING from Airtable: " + $perNotAt.Count)
$perNotAt | Select-Object key, type, title, DOI, date, pub |
  Export-Csv "$out\personalSAD_not_in_airtable.csv" -NoTypeInformation -Encoding utf8

# ---- per-Biology-collection gap table -----------------------------------------
$nameOf = @{}
foreach ($t in $tree) { if (-not $nameOf.ContainsKey($t.key)) { $nameOf[$t.key] = $t.name } }
$gapByColl = @{}; $totByColl = @{}
foreach ($b in $bio) {
  $cks = @($b.collections -split ';' | Where-Object { $_ })
  $isMiss = $false
  if (-not (Find-Match $b.title $atL)) { $isMiss = $true }
  foreach ($ck in $cks) {
    if ($totByColl.ContainsKey($ck)) { $totByColl[$ck]++ } else { $totByColl[$ck] = 1 }
    if ($isMiss) {
      if ($gapByColl.ContainsKey($ck)) { $gapByColl[$ck]++ } else { $gapByColl[$ck] = 1 }
    }
  }
}
$rows = New-Object System.Collections.ArrayList
foreach ($ck in $totByColl.Keys) {
  $n = $nameOf[$ck]; if (-not $n) { $n = $ck }
  $g = 0; if ($gapByColl.ContainsKey($ck)) { $g = $gapByColl[$ck] }
  [void]$rows.Add([pscustomobject]@{ collection = $n; key = $ck; items = $totByColl[$ck]; missing = $g; onAirtable = $totByColl[$ck] - $g })
}
$rows = $rows | Sort-Object -Property @{ Expression = 'items'; Descending = $true }
Write-Output ""
Write-Output "== Per-Biology-collection Airtable coverage (an item can sit in several collections) =="
$rows | Format-Table -AutoSize | Out-String -Width 120 | Write-Output

Write-Output "DONE"
