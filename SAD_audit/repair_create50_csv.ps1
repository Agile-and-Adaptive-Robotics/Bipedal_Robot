$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$p = "$out\airtable_create50.csv"
$rows = Get-Content $p -Encoding UTF8 | ConvertFrom-Csv

# Repair cp437-misread double-encoding: mojibake -> cp437 bytes -> UTF-8 decode
function Repair-Mojibake([string]$s) {
  if (-not $s) { return $s }
  if ($s -match '[├┐ΓÇô]' -or $s.Contains([char]0x251C)) {
    $bytes = [System.Text.Encoding]::GetEncoding(437).GetBytes($s)
    return [System.Text.Encoding]::UTF8.GetString($bytes)
  }
  return $s
}
$fixed = 0
foreach ($r in $rows) {
  $a = Repair-Mojibake $r.author
  $t = Repair-Mojibake $r.title
  if ($a -ne $r.author) { $r.author = $a; $fixed++ }
  if ($t -ne $r.title) { $r.title = $t; $fixed++ }
  # record the DOI casing actually written to Airtable for this record
  if ($r.zotero_key -eq '3SB4QSIF') { $r.DOI = '10.1523/JNEUROSCI.1082-18.2018' }
}
$rows | Export-Csv $p -NoTypeInformation -Encoding utf8
Write-Output ("cells repaired: " + $fixed)
foreach ($r in $rows) {
  if ($r.zotero_key -in @('7WF6KLGD','G358BIIV','3SB4QSIF')) {
    Write-Output ($r.zotero_key + " | " + $r.author + " | " + $r.DOI + " | " + $r.title.Substring(0, [Math]::Min(60, $r.title.Length)))
  }
}
