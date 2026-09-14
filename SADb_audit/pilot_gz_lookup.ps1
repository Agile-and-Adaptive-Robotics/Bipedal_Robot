# Look up Grillner & Zangger 1975 in the AARL group library (735051) via local API
$ErrorActionPreference = 'Continue'
$doi = '10.1016/0006-8993(75)90401-1'
$enc = [uri]::EscapeDataString($doi)
$url = "http://localhost:23119/api/groups/735051/items?format=json&q=$enc&qmode=everything&limit=25"
curl.exe -s -m 15 $url -o pilot_gz_group.json
$raw = Get-Content -Raw pilot_gz_group.json
$items = $null
try { $items = $raw | ConvertFrom-Json } catch { Write-Output "PARSE FAIL"; Write-Output $raw.Substring(0,[Math]::Min(200,$raw.Length)); exit }
if (-not $items -or @($items).Count -eq 0) { Write-Output "NO HITS in group"; exit }
Write-Output ("hits: " + @($items).Count)
foreach ($i in @($items)) {
    Write-Output ("  " + $i.key + " " + $i.data.itemType + " parent=[" + $i.data.parentItem + "] doi=[" + $i.data.DOI + "]")
}
# top-level matches
$tops = @($items | Where-Object { -not $_.data.parentItem -and $_.data.DOI -and ($_.data.DOI -replace '^https?://(dx\.)?doi\.org/','') -ieq $doi })
foreach ($t in $tops) {
    Write-Output ("TOP: " + $t.key + " | " + $t.data.title)
    $abs = $t.data.abstractNote
    if ($abs) { Write-Output ("ABSTRACT: " + $abs) } else { Write-Output "no abstract on item" }
    curl.exe -s -m 15 "http://localhost:23119/api/groups/735051/items/$($t.key)/children?format=json" -o pilot_gz_kids.json
    $kids = (Get-Content -Raw pilot_gz_kids.json) | ConvertFrom-Json
    foreach ($c in @($kids)) {
        $fn = ''
        if ($c.data.itemType -eq 'attachment') { $fn = " filename=" + $c.data.filename + " ct=" + $c.data.contentType + " link=" + $c.data.linkMode }
        Write-Output ("  child " + $c.key + " " + $c.data.itemType + $fn)
    }
}
