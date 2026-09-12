$ErrorActionPreference = "Continue"
$key = "UDdwYnJ7USmlcAJe6NSG3RbD"
$uid = 631450
$out = "D:\Github\Bipedal_Robot\SAD_audit"

$jobs = @(
  @{ k = 'HMU7XHUC'; f = "$out\cofer2010_animatlab.pdf" },
  @{ k = '5JNA8S9B'; f = "$out\buschges1995_pilocarpine.pdf" }
)
foreach ($x in $jobs) {
  $fi = Get-Item $x.f
  $md5 = (Get-FileHash $x.f -Algorithm MD5).Hash.ToLower()
  $mtime = [DateTimeOffset]::Now.ToUnixTimeMilliseconds()
  $pf = "$env:TEMP\zot_params_$($x.k).txt"
  [System.IO.File]::WriteAllText($pf, "md5=$md5&filename=$($fi.Name)&filesize=$($fi.Length)&mtime=$mtime", [System.Text.Encoding]::ASCII)
  $raw = curl.exe -s -X POST -H "Authorization: Bearer $key" -H "If-None-Match: *" -H "Content-Type: application/x-www-form-urlencoded" --data-binary "@$pf" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
  Write-Output ("$($x.k) auth: " + $raw)
  $a = $raw | ConvertFrom-Json
  $tag = $md5
  if ($a.uploadKey) { $tag = $a.uploadKey }
  if (-not $a.exists) {
    Write-Output "$($x.k): needs real upload (not exists case) - running upload"
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
    Write-Output ("  upload status: " + $code)
  }
  $reg = curl.exe -s -o "$env:TEMP\zot_reg_$($x.k).txt" -w "%{http_code}" -X POST -H "Authorization: Bearer $key" -H "If-None-Match: $tag" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
  Write-Output ("$($x.k) register: status $reg  body: " + (Get-Content "$env:TEMP\zot_reg_$($x.k).txt" -Raw -ErrorAction SilentlyContinue))
}
# verify children
foreach ($par in @('A924UV2P', '7M7SRMUQ')) {
  curl.exe -s -o "$env:TEMP\zot_kids_$par.json" "https://api.zotero.org/users/$uid/items/$par/children?format=json"
  $kids = (Get-Content "$env:TEMP\zot_kids_$par.json" -Raw -Encoding UTF8) | ConvertFrom-Json
  foreach ($c in @($kids)) { Write-Output ("VERIFY $par -> " + $c.key + " [" + $c.data.itemType + "] " + $c.data.filename) }
}
Write-Output "DONE"
