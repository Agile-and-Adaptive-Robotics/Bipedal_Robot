$ErrorActionPreference = "Continue"
$key = "UDdwYnJ7USmlcAJe6NSG3RbD"
$uid = 631450
$pf = "$env:TEMP\zot_dbg_params.txt"
$content = "md5=ae45cd9cd28aa50bbda83a9dbd48fdca&filename=cofer2010_animatlab.pdf&filesize=1656317&mtime=1789158393321"
[System.IO.File]::WriteAllText($pf, $content, [System.Text.Encoding]::ASCII)
Write-Output ("body file: [" + ([System.IO.File]::ReadAllText($pf)) + "]")
$raw = curl.exe -s -X POST -H "Authorization: Bearer $key" -H "If-None-Match: *" -H "Content-Type: application/x-www-form-urlencoded" --data-binary "@$pf" "https://api.zotero.org/users/$uid/items/HMU7XHUC/file"
Write-Output ("RESPONSE: " + $raw)
