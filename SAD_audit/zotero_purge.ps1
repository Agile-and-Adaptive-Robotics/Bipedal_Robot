param([string]$ApiKey = "")
$ErrorActionPreference = "Stop"
if (-not $ApiKey) { Write-Output "ERROR: pass -ApiKey"; exit 1 }
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$uid = 631450

$elig = Get-Content "$out\personal_SAD_attachment_purge_eligible.csv" -Encoding UTF8 | ConvertFrom-Csv
Write-Output ("eligible list: " + @($elig).Count + " attachments on " + @($elig | Select-Object -ExpandProperty parentKey -Unique).Count + " items")

# 1. re-verify current parent of every attachment (batches of 50, file-buffered reads)
$verify = @{}   # attKey -> parentKey (current)
$keys = @($elig | Select-Object -ExpandProperty attKey)
$tmpJson = "$env:TEMP\zot_purge_get.json"
for ($i = 0; $i -lt $keys.Count; $i += 50) {
  $batch = $keys[$i..([Math]::Min($i + 49, $keys.Count - 1))]
  $q = ($batch -join ',')
  curl.exe -s -o $tmpJson -H "Authorization: Bearer $ApiKey" "https://api.zotero.org/users/$uid/items?format=json&itemKey=$q&limit=50"
  $r = Get-Content $tmpJson -Raw -Encoding UTF8
  $items = $r | ConvertFrom-Json
  foreach ($it in @($items)) { $verify[$it.key] = $it.data.parentItem }
}
Write-Output ("re-verified current state: " + $verify.Count + " of " + $keys.Count + " attachments found")

$toDelete = New-Object System.Collections.ArrayList
$skipped = 0
foreach ($e in $elig) {
  if ($verify.ContainsKey($e.attKey) -and $verify[$e.attKey] -eq $e.parentKey) {
    [void]$toDelete.Add($e.attKey)
  } else {
    $skipped++
    $cur = ''; if ($verify.ContainsKey($e.attKey)) { $cur = $verify[$e.attKey] }
    Write-Output ("  SKIP drift: " + $e.attKey + " csv-parent=" + $e.parentKey + " current-parent=" + $cur)
  }
}
Write-Output ("to delete: " + $toDelete.Count + "  skipped for drift: " + $skipped)

# 2. DELETE in batches of 50 (moves to trash - recoverable)
$okCount = 0; $failCount = 0
for ($i = 0; $i -lt $toDelete.Count; $i += 50) {
  $batch = $toDelete[$i..([Math]::Min($i + 49, $toDelete.Count - 1))]
  $q = ($batch -join ',')
  $code = curl.exe -s -o "$env:TEMP\zot_del_resp.txt" -w "%{http_code}" -X DELETE -H "Authorization: Bearer $ApiKey" "https://api.zotero.org/users/$uid/items?itemKey=$q"
  if ($code -eq '204') { $okCount += $batch.Count; Write-Output ("  DELETE batch ok (" + $batch.Count + " keys, status " + $code + ")") }
  else { $failCount += $batch.Count; $body = Get-Content "$env:TEMP\zot_del_resp.txt" -Raw; Write-Output ("  DELETE batch FAILED status " + $code + " : " + $body.Substring(0, [Math]::Min(300, $body.Length))) }
}
Write-Output ("deleted (trashed): " + $okCount + "  failed: " + $failCount)

# 3. post-verify: trashed items should now report deletedTime
$sample = @($toDelete | Select-Object -First 10)
if ($sample.Count -gt 0) {
  $q = ($sample -join ',')
  curl.exe -s -o $tmpJson -H "Authorization: Bearer $ApiKey" "https://api.zotero.org/users/$uid/items?format=json&itemKey=$q&limit=50"
  $r = Get-Content $tmpJson -Raw -Encoding UTF8
  $items = $r | ConvertFrom-Json
  $trashed = 0
  foreach ($it in @($items)) { if ($it.data.deleted) { $trashed++ } }
  Write-Output ("post-verify sample: " + $trashed + "/" + $sample.Count + " now carry deletedTime (trash)")
}
Write-Output "DONE"
