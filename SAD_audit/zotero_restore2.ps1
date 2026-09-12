$ErrorActionPreference = "Stop"
$key = "UDdwYnJ7USmlcAJe6NSG3RbD"
$uid = 631450
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$tmp = "$env:TEMP\zot_restore"

# 1. download the two PDFs from Airtable CDN
curl.exe -s -o "$tmp_cofer.pdf" -w "cofer dl: %{http_code} %{size_download}`n" "https://v5.airtableusercontent.com/v3/u/57/57/1789171200000/AHreQXB0w_3XgocXqTP_zQ/jyPBxEQPIZP71NXSaTwydgOfEm6kOCumpEaBVtzrjf94MtlOec5Nl0E1NOb5WOmOvklaicEhILkj2jtf140NlNL3B4AEhLMARLMruLDD3kzDyMwu_2qzz89NJZxoWCHTcsufA972bEZThuqK3Q6qB71KCVD_5m84-3jyiKoYLJfwGFdZOzL9cQygziw28RSAu68kyag8zjbksw6B3rLMUg/RC2FC5KJ6q1oiX-e9--ZiCQSXVUxVAB5kpPO7bbgtmo"
Move-Item "$tmp_cofer.pdf" "$out\cofer2010_animatlab.pdf" -Force
curl.exe -s -o "$tmp_bues.pdf" -w "bueschges dl: %{http_code} %{size_download}`n" "https://v5.airtableusercontent.com/v3/u/57/57/1789171200000/ROp1klimqFf1kJKItMrggQ/BOqIHW1pQC-n4_iAEXemjKroFte7cyJ3iTj3rdVSe-FTQWr1bPm88VLP57XsN0KYM6V4SfOV5kgxqqp0rtJtc9GEDRv00P0-BjlwFfDgQOHFIrtJ89OH8-_h6KX70Grhi1kFCAeMxBOHDAO-pxxA_uCkoNLmT40c8el-dzdH5okye61XhMaEldyVKzr1MahBH5ntJEwckajEm3rhn2CfEg/n0hqSuPBN0fbEvafHAJhyoMUn7Ijqfl_4NsiS3WgOUA"
Move-Item "$tmp_bues.pdf" "$out\buschges1995_pilocarpine.pdf" -Force

# 2. create child attachment items on the Zotero parents
$items = @(
  @{ itemType = 'attachment'; parentItem = 'A924UV2P'; linkMode = 'imported_file'; title = 'Cofer et al. - 2010 - AnimatLab A 3D graphics environment for neuromech.pdf'; contentType = 'application/pdf'; filename = 'Cofer et al. - 2010 - AnimatLab A 3D graphics environment for neuromech.pdf'; tags = @(); relations = @{} },
  @{ itemType = 'attachment'; parentItem = '7M7SRMUQ'; linkMode = 'imported_file'; title = 'Büschges et al. - 1995 - Rhythmic Patterns in the Thoracic Nerve Cord of th.pdf'; contentType = 'application/pdf'; filename = 'Büschges et al. - 1995 - Rhythmic Patterns in the Thoracic Nerve Cord of th.pdf'; tags = @(); relations = @{} }
) | ConvertTo-Json -Depth 5
[System.IO.File]::WriteAllText("$out\restore2_items.json", $items, (New-Object System.Text.UTF8Encoding($false)))
$resp = curl.exe -s -X POST "https://api.zotero.org/users/$uid/items" -H "Authorization: Bearer $key" -H "Content-Type: application/json" --data-binary "@$out\restore2_items.json"
[System.IO.File]::WriteAllText("$env:TEMP\zot_restore2_resp.json", $resp, [System.Text.Encoding]::UTF8)
$r = $resp | ConvertFrom-Json
$newKeys = @()
foreach ($p in $r.success.PSObject.Properties) { $newKeys += $p.Value.key; Write-Output ("created attachment item: " + $p.Value.key + " parent=" + $p.Value.data.parentItem) }
if ($newKeys.Count -ne 2) { Write-Output "FAILED to create both items"; exit 1 }

# 3. get upload authorization + POST files (Zotero file upload flow)
$files = @(
  @{ k = $newKeys[0]; f = "$out\cofer2010_animatlab.pdf"; md5 = (Get-FileHash "$out\cofer2010_animatlab.pdf" -Algorithm MD5).Hash.ToLower(); t = (Get-Item "$out\cofer2010_animatlab.pdf").LastWriteTimeUtc.ToString('R') },
  @{ k = $newKeys[1]; f = "$out\buschges1995_pilocarpine.pdf"; md5 = (Get-FileHash "$out\buschges1995_pilocarpine.pdf" -Algorithm MD5).Hash.ToLower(); t = (Get-Item "$out\buschges1995_pilocarpine.pdf").LastWriteTimeUtc.ToString('R') }
)
foreach ($x in $files) {
  $auth = curl.exe -s -X POST -H "Authorization: Bearer $key" -H "Content-Type: application/x-www-form-urlencoded" -d "md5=$($x.md5)&filename=$( [IO.Path]::GetFileName($x.f) )&filesize=$((Get-Item $x.f).Length)&mtime=$([DateTimeOffset]::Now.ToUnixTimeMilliseconds())" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
  $a = $auth | ConvertFrom-Json
  if (-not $a.uploadKey) { Write-Output ("AUTH FAILED for $($x.k): " + $auth.Substring(0, [Math]::Min(300, $auth.Length))); continue }
  $boundary = "----zot$([guid]::NewGuid().ToString('N'))"
  $uploadKey = $a.uploadKey
  # build multipart body exactly as Zotero expects: uploadKey field first, then file
  $bodyFile = "$env:TEMP\zot_upload_body.bin"
  $pre = "--$boundary`r`nContent-Disposition: form-data; name=`"uploadKey`"`r`n`r`n$uploadKey`r`n--$boundary`r`nContent-Disposition: form-data; name=`"file`"; filename=`"$( [IO.Path]::GetFileName($x.f) )`"`r`nContent-Type: application/pdf`r`n`r`n"
  $post = "`r`n--$boundary--`r`n"
  $preB = [System.Text.Encoding]::UTF8.GetBytes($pre)
  $postB = [System.Text.Encoding]::UTF8.GetBytes($post)
  $fileB = [System.IO.File]::ReadAllBytes($x.f)
  $ms = New-Object System.IO.MemoryStream
  $ms.Write($preB, 0, $preB.Length); $ms.Write($fileB, 0, $fileB.Length); $ms.Write($postB, 0, $postB.Length)
  [System.IO.File]::WriteAllBytes($bodyFile, $ms.ToArray())
  $code = curl.exe -s -o "$env:TEMP\zot_up_resp.txt" -w "%{http_code}" -X POST -H "Authorization: Bearer $key" -H "Content-Type: multipart/form-data; boundary=$boundary" --data-binary "@$bodyFile" $a.url
  Write-Output ("upload $($x.k): status $code")
  # 4. register upload
  $reg = curl.exe -s -o "$env:TEMP\zot_reg_resp.txt" -w "%{http_code}" -X POST -H "Authorization: Bearer $key" -H "If-None-Match: $($a.uploadKey)" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
  Write-Output ("register $($x.k): status $reg")
}
Write-Output "DONE"
