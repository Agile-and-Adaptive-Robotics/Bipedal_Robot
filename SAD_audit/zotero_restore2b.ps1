$ErrorActionPreference = "Stop"
$key = "UDdwYnJ7USmlcAJe6NSG3RbD"
$uid = 631450
$out = "D:\Github\Bipedal_Robot\SAD_audit"

$r = Get-Content "$env:TEMP\zot_restore2_resp.json" -Raw -Encoding UTF8 | ConvertFrom-Json
$pairs = @()
foreach ($p in $r.successful.PSObject.Properties) {
  $pairs += @{ k = $p.Value.key; parent = $p.Value.data.parentItem }
}
Write-Output ("created items: " + (($pairs | ForEach-Object { $_.k + '->' + $_.parent }) -join ' ; '))
if ($pairs.Count -ne 2) { Write-Output "unexpected count"; exit 1 }
$cofer = $pairs | Where-Object { $_.parent -eq 'A924UV2P' }
$bues  = $pairs | Where-Object { $_.parent -eq '7M7SRMUQ' }

$jobs = @(
  @{ k = $cofer.k; f = "$out\cofer2010_animatlab.pdf" },
  @{ k = $bues.k;  f = "$out\buschges1995_pilocarpine.pdf" }
)
foreach ($x in $jobs) {
  $fi = Get-Item $x.f
  $md5 = (Get-FileHash $x.f -Algorithm MD5).Hash.ToLower()
  $mtime = [DateTimeOffset]::Now.ToUnixTimeMilliseconds()
  $dfile = "$env:TEMP\zot_d_$($x.k).txt"
  # -d with data read from file avoids shell &-quoting hazards
  "md5=$md5`nfilename=$($fi.Name)`nfilesize=$($fi.Length)`nmtime=$mtime" | Out-File -Encoding ascii "$env:TEMP\zot_d_params.txt"
  $params = Get-Content "$env:TEMP\zot_d_params.txt" -Raw
  $params = $params.Replace("`r`n", "&").TrimEnd("&")
  [System.IO.File]::WriteAllText("$env:TEMP\zot_d_params_amp.txt", $params, [System.Text.Encoding]::ASCII)
  $authRaw = curl.exe -s -X POST -H "Authorization: Bearer $key" -H "If-None-Match: *" -H "Content-Type: application/x-www-form-urlencoded" --data-binary "@$env:TEMP\zot_d_params_amp.txt" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
  [System.IO.File]::WriteAllText("$dfile", $authRaw, [System.Text.Encoding]::UTF8)
  $a = $authRaw | ConvertFrom-Json
  if (-not $a.uploadKey) { Write-Output ("AUTH FAILED $($x.k): " + $authRaw.Substring(0, [Math]::Min(300, $authRaw.Length))); continue }
  Write-Output ("auth ok $($x.k): exists=" + $a.exists + " url=" + ($a.url -ne $null))
  if ($a.exists) { Write-Output ("  file already exists server-side; registering only") }
  if (-not $a.exists) {
    $boundary = "----zot" + [guid]::NewGuid().ToString('N')
    $bodyFile = "$env:TEMP\zot_body_$($x.k).bin"
    $pre = "--$boundary`r`nContent-Disposition: form-data; name=`"uploadKey`"`r`n`r`n$($a.uploadKey)`r`n--$boundary`r`nContent-Disposition: form-data; name=`"file`"; filename=`"$($fi.Name)`"`r`nContent-Type: application/pdf`r`n`r`n"
    $post = "`r`n--$boundary--`r`n"
    $preB = [System.Text.Encoding]::ASCII.GetBytes($pre)
    $postB = [System.Text.Encoding]::ASCII.GetBytes($post)
    $fileB = [System.IO.File]::ReadAllBytes($x.f)
    $ms = New-Object System.IO.MemoryStream
    $ms.Write($preB, 0, $preB.Length); $ms.Write($fileB, 0, $fileB.Length); $ms.Write($postB, 0, $postB.Length)
    [System.IO.File]::WriteAllBytes($bodyFile, $ms.ToArray())
    $code = curl.exe -s -o "$env:TEMP\zot_up_resp.txt" -w "%{http_code}" -X POST -H "Authorization: Bearer $key" -H "Content-Type: multipart/form-data; boundary=$boundary" --data-binary "@$bodyFile" $a.url
    Write-Output ("  upload: status $code")
  }
  $reg = curl.exe -s -o "$env:TEMP\zot_reg_resp.txt" -w "%{http_code}" -X POST -H "Authorization: Bearer $key" -H "If-None-Match: $($a.uploadKey)" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
  Write-Output ("  register: status $reg")
}
# verify both parents now have a file child
foreach ($par in @('A924UV2P', '7M7SRMUQ')) {
  curl.exe -s -o "$env:TEMP\zot_kids_$par.json" "https://api.zotero.org/users/$uid/items/$par/children?format=json"
  $kids = (Get-Content "$env:TEMP\zot_kids_$par.json" -Raw -Encoding UTF8) | ConvertFrom-Json
  foreach ($c in @($kids)) { Write-Output ("VERIFY $par child: " + $c.key + " " + $c.data.itemType + " " + $c.data.title.Substring(0, [Math]::Min(50, $c.data.title.Length))) }
}
Write-Output "DONE"
