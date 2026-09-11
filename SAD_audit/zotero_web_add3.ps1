param([string]$ApiKey = "")
$ErrorActionPreference = "Stop"
if (-not $ApiKey) { Write-Output "ERROR: pass -ApiKey"; exit 1 }
$out = "D:\Github\Bipedal_Robot\SAD_audit"

# 1. key info
$info = curl.exe -s -H "Authorization: Bearer $ApiKey" "https://api.zotero.org/keys/current"
Write-Output ("KEYINFO: " + $info)
$ki = $info | ConvertFrom-Json
if (-not $ki.access.user.write) { Write-Output "ERROR: key has NO personal-library write access"; exit 1 }
$uid = $ki.userID
Write-Output ("userID: " + $uid)

# 2. add the 3 items (payload file already prepared, UTF-8)
$respFile = "$env:TEMP\zot_add3_resp.json"
$code = curl.exe -s -o $respFile -w "%{http_code}" -X POST "https://api.zotero.org/users/$uid/items" -H "Authorization: Bearer $ApiKey" -H "Content-Type: application/json" --data-binary "@$out\add3_items.json"
Write-Output ("POST status: " + $code)
$resp = Get-Content $respFile -Raw -Encoding UTF8
Write-Output ("RESPONSE: " + $resp.Substring(0, [Math]::Min(1500, $resp.Length)))

# 3. verify: fetch the 3 DOIs back
foreach ($doi in @('10.1113/jphysiol.1971.sp009487','10.1016/j.asd.2015.07.001','10.1016/s1388-2457(03)00120-2')) {
  $r = curl.exe -s -H "Authorization: Bearer $ApiKey" "https://api.zotero.org/users/$uid/items?format=json&q=$doi&doi=true"
  $items = $r | ConvertFrom-Json
  foreach ($i in @($items)) {
    Write-Output ("VERIFY " + $i.key + " :: " + $i.data.title + " :: collections=" + (@($i.data.collections) -join ','))
  }
}
Write-Output "DONE"
