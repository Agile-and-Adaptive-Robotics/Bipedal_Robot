$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$p = "$out\airtable_create50.csv"
$rows = Get-Content $p -Encoding UTF8 | ConvertFrom-Csv

$mojiChars = @([char]0x251C, [char]0x2510, [char]0x0393)  # box-drawing + Greek Gamma: cp437-misread artifacts
function Has-Mojibake([string]$s) {
  foreach ($c in $mojiChars) { if ($s.IndexOf($c) -ge 0) { return $true } }
  return $false
}
$fixed = 0
foreach ($r in $rows) {
  foreach ($prop in @('author', 'title')) {
    $s = $r.$prop
    if ($s -and (Has-Mojibake $s)) {
      $bytes = [System.Text.Encoding]::GetEncoding(437).GetBytes($s)
      $r.$prop = [System.Text.Encoding]::UTF8.GetString($bytes)
      $fixed++
    }
  }
}
$rows | Export-Csv $p -NoTypeInformation -Encoding utf8
Write-Output ("cells repaired this pass: " + $fixed)
# verify bytes for the known case
$line = (Get-Content $p -Encoding UTF8) | Where-Object { $_ -like '*3SB4QSIF*' } | Select-Object -First 1
$i = $line.IndexOf('Left')
Write-Output ('3SB4QSIF title bytes after repair: ' + (([System.Text.Encoding]::UTF8.GetBytes($line.Substring($i, 6)) | ForEach-Object { '{0:X2}' -f $_ }) -join ' ') + '   (want 2D somewhere / no CE93)')
