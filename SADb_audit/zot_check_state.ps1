$ErrorActionPreference = "Continue"
# check both new attachment items via LOCAL Zotero (synced view) and web API
foreach ($k in @('HMU7XHUC', '5JNA8S9B')) {
  curl.exe -s -o "$env:TEMP\zc_local_$k.json" "http://localhost:23119/api/users/0/items/$k"
  $raw = Get-Content "$env:TEMP\zc_local_$k.json" -Raw -Encoding UTF8
  Write-Output ("=== $k local response head: " + $raw.Substring(0, [Math]::Min(80, $raw.Length)))
  try {
    $it = $raw | ConvertFrom-Json
    Write-Output ("  local md5: '" + $it.data.md5 + "'  mtime: '" + $it.data.mtime + "'  filename: " + $it.data.filename)
  } catch { Write-Output "  parse failed" }
}
