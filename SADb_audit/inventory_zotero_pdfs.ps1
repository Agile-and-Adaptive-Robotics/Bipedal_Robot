# Inventory every PDF in personal + AARL group Zotero, match to corpus DOIs,
# copy matches into SADb_audit\pdf_staging\<Surname Year>__<key>.pdf
# PS 5.1-safe: read curl output with -Encoding UTF8, ASCII-safe output strings.
$ErrorActionPreference = 'Continue'
$storageRoot = 'C:\Users\Ben Bolen\Zotero\storage'
$out = New-Object System.Collections.Generic.List[object]
$att = New-Object System.Collections.Generic.List[object]

function Get-Surname($creators) {
    foreach ($c in @($creators)) {
        if ($c.lastName) { return [string]$c.lastName }
        if ($c.name) { $p = ([string]$c.name).Trim().Split(' '); return $p[$p.Count-1] }
    }
    return ''
}
function Clean-Name($s) {
    if (-not $s) { return 'Unknown' }
    $s = $s -replace '[\\/:*?"<>|]', ''
    return $s.Trim()
}

foreach ($scope in @('users/0','groups/735051')) {
    $start = 0
    while ($true) {
        curl.exe -s -m 40 "http://localhost:23119/api/$scope/items?format=json&limit=100&start=$start&itemType=attachment" -o zt_att.json
        $items = $null
        try { $items = (Get-Content -Raw -Encoding UTF8 zt_att.json) | ConvertFrom-Json } catch { }
        if (-not $items -or @($items).Count -eq 0) { break }
        foreach ($it in @($items)) {
            if ($it.data.contentType -ne 'application/pdf') { continue }
            $att.Add([pscustomobject]@{
                scope = $scope; key = $it.key; parent = [string]$it.data.parentItem
                filename = [string]$it.data.filename
                path = Join-Path $storageRoot $it.key
            })
        }
        $start += 100
        if ($start -gt 2500) { break }
    }
}
Write-Output ("pdf attachments found: " + $att.Count)

# parent item lookup: fetch parents in batches by key to get DOI/title/creators
$parents = @{}
$keys = @($att | ForEach-Object { $_.parent } | Where-Object { $_ } | Sort-Object -Unique)
Write-Output ("distinct parents: " + $keys.Count)
$i = 0
while ($i -lt $keys.Count) {
    $batch = @($keys[$i..([Math]::Min($i+48, $keys.Count-1))])
    $keyParam = ($batch | ForEach-Object { "key:" + $_ }) -join '||'
    $enc = [uri]::EscapeDataString($keyParam)
    curl.exe -s -m 60 "http://localhost:23119/api/users/0/items?format=json&limit=100&itemKey=$enc" -o zt_par.json
    $items = $null
    try { $items = (Get-Content -Raw -Encoding UTF8 zt_par.json) | ConvertFrom-Json } catch { }
    if ($items) {
        foreach ($it in @($items)) {
            $sn = Get-Surname $it.data.creators
            $parents[$it.key] = [pscustomobject]@{ doi = ([string]$it.data.DOI).Trim(); title = [string]$it.data.title; surname = $sn; year = [string]$it.data.date }
        }
    }
    # group parents need the group scope
    curl.exe -s -m 60 "http://localhost:23119/api/groups/735051/items?format=json&limit=100&itemKey=$enc" -o zt_par2.json
    $items2 = $null
    try { $items2 = (Get-Content -Raw -Encoding UTF8 zt_par2.json) | ConvertFrom-Json } catch { }
    if ($items2) {
        foreach ($it in @($items2)) {
            if (-not $parents.ContainsKey($it.key)) {
                $sn = Get-Surname $it.data.creators
                $parents[$it.key] = [pscustomobject]@{ doi = ([string]$it.data.DOI).Trim(); title = [string]$it.data.title; surname = $sn; year = [string]$it.data.date }
            }
        }
    }
    $i += 50
    Start-Sleep -Milliseconds 200
}
Write-Output ("parents resolved: " + $parents.Count)

$staging = 'D:\Github\Bipedal_Robot\SADb_audit\pdf_staging'
New-Item -ItemType Directory -Force -Path $staging | Out-Null
$copied = 0
foreach ($a in $att) {
    $p = $parents[$a.parent]
    if (-not $p) { continue }
    $doi = ($p.doi -replace '^https?://(dx\.)?doi\.org/','').ToLower()
    $file = Join-Path $a.path $a.filename
    $exists = Test-Path -LiteralPath $file
    $yr = $p.year
    if ($yr -match '(\d{4})') { $yr = $Matches[1] } else { $yr = '' }
    $out.Add([pscustomobject]@{
        library = $a.scope; zotero_parent = $a.parent; att_key = $a.key
        doi = $doi; surname = (Clean-Name $p.surname); year = $yr
        title = ($p.title -replace '\s+', ' ')
        file_exists = $exists; local_path = $file
    })
    if ($exists -and $doi) {
        $dest = Join-Path $staging ((Clean-Name $p.surname) + '_' + $yr + '__' + $a.key + '.pdf')
        if (-not (Test-Path -LiteralPath $dest)) {
            Copy-Item -LiteralPath $file -Destination $dest -Force
            $copied++
        }
    }
}
$out | Export-Csv -Path 'D:\Github\Bipedal_Robot\SADb_audit\pdf_inventory.csv' -NoTypeInformation -Encoding UTF8
Write-Output ("inventory rows: " + $out.Count + "  copied to staging: " + $copied)
