$ErrorActionPreference = "Continue"
$ApiKey = $args[0]
$uid = 631450
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$elig = Get-Content "$out\personal_SAD_attachment_purge_eligible.csv" -Encoding UTF8 | ConvertFrom-Csv
$keys = @($elig | Select-Object -ExpandProperty attKey)
Write-Output ("keys: " + $keys.Count + "  types: " + ($keys | Select-Object -First 3 | ForEach-Object { $_.GetType().Name }) -join ' ')
$tmpJson = "$env:TEMP\zot_dbg_batch.json"
$verify = @{}
for ($i = 0; $i -lt $keys.Count; $i += 50) {
  $hi = [Math]::Min($i + 49, $keys.Count - 1)
  $batch = $keys[$i..$hi]
  $q = ($batch -join ',')
  curl.exe -s -o $tmpJson -w "batch $i status: %{http_code}`n" -H "Authorization: Bearer $ApiKey" "https://api.zotero.org/users/$uid/items?format=json&itemKey=$q&limit=50"
  $fi = Get-Item $tmpJson
  Write-Output ("  file bytes: " + $fi.Length)
  $r = Get-Content $tmpJson -Raw -Encoding UTF8
  Write-Output ("  raw head: " + $r.Substring(0, [Math]::Min(100, $r.Length)).Replace("`n", " | ").Replace("`r", ""))
  $parsed = @($r | ConvertFrom-Json)
  Write-Output ("  parsed: " + $parsed.Count)
  foreach ($it in $parsed) { $verify[$it.key] = $it.data.parentItem }
}
Write-Output ("verify total: " + $verify.Count)
