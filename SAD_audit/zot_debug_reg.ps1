$ErrorActionPreference = "Continue"
$key = "UDdwYnJ7USmlcAJe6NSG3RbD"
$uid = 631450
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$f = "$out\cofer2010_animatlab.pdf"
$k = 'HMU7XHUC'
$fi = Get-Item $f
$md5 = (Get-FileHash $f -Algorithm MD5).Hash.ToLower()

# A) auth with If-None-Match: md5
$pf = "$env:TEMP\zot_auth_try.txt"
[System.IO.File]::WriteAllText($pf, "md5=$md5&filename=$($fi.Name)&filesize=$($fi.Length)&mtime=$([DateTimeOffset]::Now.ToUnixTimeMilliseconds())", [System.Text.Encoding]::ASCII)
$raw = curl.exe -s -X POST -H "Authorization: Bearer $key" -H "If-None-Match: $md5" -H "Content-Type: application/x-www-form-urlencoded" --data-binary "@$pf" "https://api.zotero.org/users/$uid/items/$k/file"
Write-Output ("A auth(md5-tag): " + $raw)

# B) register with body + md5 tag
$reg = curl.exe -s -o "$env:TEMP\zot_reg_try.txt" -w "%{http_code}" -X POST -H "Authorization: Bearer $key" -H "If-None-Match: $md5" -H "Content-Type: application/x-www-form-urlencoded" --data-binary "@$pf" "https://api.zotero.org/users/$uid/items/$k/file"
Write-Output ("B register(md5-tag + body): " + $reg + "  " + (Get-Content "$env:TEMP\zot_reg_try.txt" -Raw -ErrorAction SilentlyContinue))

# C) does the item actually already have the file? check item md5 via item view
curl.exe -s -o "$env:TEMP\zot_itemchk.json" -H "Authorization: Bearer $key" "https://api.zotero.org/users/$uid/items/$k?format=json"
$it = (Get-Content "$env:TEMP\zot_itemchk.json" -Raw -Encoding UTF8) | ConvertFrom-Json
Write-Output ("C item md5: " + $it.data.md5 + "  mtime: " + $it.data.mtime)
