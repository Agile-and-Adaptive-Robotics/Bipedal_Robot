$ErrorActionPreference = "Continue"
$key = "UDdwYnJ7USmlcAJe6NSG3RbD"
$uid = 631450
$out = "D:\Github\Bipedal_Robot\SAD_audit"

# Both files already exist in Zotero storage ({exists:1}) -> register only.
# Registration: POST .../file with If-None-Match: <md5>, empty body.
foreach ($x in @(
  @{ k = 'HMU7XHUC'; f = "$out\cofer2010_animatlab.pdf" },
  @{ k = '5JNA8S9B'; f = "$out\buschges1995_pilocarpine.pdf" }
)) {
  $md5 = (Get-FileHash $x.f -Algorithm MD5).Hash.ToLower()
  $reg = curl.exe -s -o "$env:TEMP\zot_reg_$($x.k).txt" -w "%{http_code}" -X POST -H "Authorization: Bearer $key" -H "If-None-Match: $md5" -H "Content-Length: 0" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
  Write-Output ("$($x.k) register(md5): status $reg  body: " + (Get-Content "$env:TEMP\zot_reg_$($x.k).txt" -Raw -ErrorAction SilentlyContinue))
  if ($reg -ne '204') {
    # fall back: re-auth (no If-None-Match) -> get uploadKey -> register with it
    $fi = Get-Item $x.f
    $pf = "$env:TEMP\zot_params_$($x.k).txt"
    [System.IO.File]::WriteAllText($pf, "md5=$md5&filename=$($fi.Name)&filesize=$($fi.Length)&mtime=$([DateTimeOffset]::Now.ToUnixTimeMilliseconds())", [System.Text.Encoding]::ASCII)
    $raw = curl.exe -s -X POST -H "Authorization: Bearer $key" -H "Content-Type: application/x-www-form-urlencoded" --data-binary "@$pf" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
    Write-Output ("  re-auth: " + $raw)
    $a = $raw | ConvertFrom-Json
    if ($a.uploadKey) {
      $reg2 = curl.exe -s -o "$env:TEMP\zot_reg2_$($x.k).txt" -w "%{http_code}" -X POST -H "Authorization: Bearer $key" -H "If-None-Match: $($a.uploadKey)" -H "Content-Length: 0" "https://api.zotero.org/users/$uid/items/$($x.k)/file"
      Write-Output ("  register(uploadKey): status $reg2")
    }
  }
}
# verify children (with auth)
foreach ($par in @('A924UV2P', '7M7SRMUQ')) {
  curl.exe -s -o "$env:TEMP\zot_kids_$par.json" -H "Authorization: Bearer $key" "https://api.zotero.org/users/$uid/items/$par/children?format=json"
  $raw = Get-Content "$env:TEMP\zot_kids_$par.json" -Raw -Encoding UTF8
  $kids = $raw | ConvertFrom-Json
  foreach ($c in @($kids)) { Write-Output ("VERIFY $par -> " + $c.key + " [" + $c.data.itemType + "] " + $c.data.filename) }
}
Write-Output "DONE"
