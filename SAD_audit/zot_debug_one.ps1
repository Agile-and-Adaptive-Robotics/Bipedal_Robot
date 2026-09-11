$ErrorActionPreference = "Stop"
$key = $args[0]
# single key
curl.exe -s -o "$env:TEMP\zot_one.json" -w "single status: %{http_code}`n" -H "Authorization: Bearer $key" "https://api.zotero.org/users/631450/items?format=json&itemKey=375JYEPB"
$j = Get-Content "$env:TEMP\zot_one.json" -Raw -Encoding UTF8
Write-Output ("len: " + $j.Length)
Write-Output $j.Substring(0, [Math]::Min(500, $j.Length))
# ten keys
$ten = "375JYEPB,WL5PV4ZW,8ZXJ7A74,KUYRMXC7,MLD55F9G,PEIW67K9,CJY2QX2E,8SQFYAE4,BH794RAT,54QAT6SJ"
curl.exe -s -o "$env:TEMP\zot_ten.json" -w "ten status: %{http_code}`n" -H "Authorization: Bearer $key" "https://api.zotero.org/users/631450/items?format=json&itemKey=$ten"
$t = Get-Content "$env:TEMP\zot_ten.json" -Raw -Encoding UTF8
Write-Output ("ten len: " + $t.Length)
$items = $t | ConvertFrom-Json
Write-Output ("ten parsed count: " + @($items).Count)
