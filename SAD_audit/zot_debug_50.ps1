$ErrorActionPreference = "Continue"
$key = "UDdwYnJ7USmlcAJe6NSG3RbD"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$elig = Get-Content "$out\personal_SAD_attachment_purge_eligible.csv" -Encoding UTF8 | ConvertFrom-Csv
$keys = @($elig | Select-Object -ExpandProperty attKey)
Write-Output ("total keys: " + $keys.Count)
$q = ($keys[0..49] -join ',')
Write-Output ("url key count: " + ($q -split ',').Count + "  url len: " + ("https://api.zotero.org/users/631450/items?format=json&itemKey=" + $q).Length)
curl.exe -s -o "$env:TEMP\zot_50.json" -w "50-status: %{http_code}`n" -H "Authorization: Bearer $key" "https://api.zotero.org/users/631450/items?format=json&itemKey=$q"
$t = Get-Content "$env:TEMP\zot_50.json" -Raw -Encoding UTF8
Write-Output ("resp len: " + $t.Length)
Write-Output ("head: " + $t.Substring(0, [Math]::Min(200, $t.Length)))
try { $items = $t | ConvertFrom-Json; Write-Output ("parsed count: " + @($items).Count) } catch { Write-Output ("parse error: " + $_.Exception.Message) }
