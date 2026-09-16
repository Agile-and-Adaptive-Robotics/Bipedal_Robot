# v2: paginate ALL top-level items in both libraries once, build key->meta map,
# then join the 2448 PDF attachments locally. Copy matches to pdf_staging.
$ErrorActionPreference = 'Continue'
$storageRoot = 'C:\Users\Ben Bolen\Zotero\storage'

function Get-Surname($creators) {
    foreach ($c in @($creators)) {
        if ($c.lastName) { return [string]$c.lastName }
        if ($c.name) { $p = ([string]$c.name).Trim().Split(' '); return $p[$p.Count-1] }
    }
    return ''
}
function Clean-Name($s) {
    if (-not $s) { return 'Unknown' }
    return ($s -replace '[\\/:*?"<>|]', '').Trim()
}

$meta = @{}
foreach ($scope in @('users/0','groups/735051')) {
    $start = 0
    while ($true) {
        curl.exe -s -m 60 "http://localhost:23119/api/$scope/items?format=json&limit=100&start=$start&itemType=-attachment" -o zt_top.json
        $items = $null
        try { $items = (Get-Content -Raw -Encoding UTF8 zt_top.json) | ConvertFrom-Json } catch { }
        if (-not $items -or @($items).Count -eq 0) { break }
        foreach ($it in @($items)) {
            if ($it.data.parentItem) { continue }
            if (-not $meta.ContainsKey($it.key)) {
                $meta[$it.key] = [pscustomobject]@{
                    doi = ([string]$it.data.DOI).Trim()
                    title = ([string]$it.data.title)
                    surname = (Get-Surname $it.data.creators)
                    year = ([string]$it.data.date)
                }
            }
        }
        $start += 100
        if ($start -gt 5000) { break }
    }
    Write-Output ("$scope top-level mapped: " + $meta.Count)
}

# attachments again
$att = New-Object System.Collections.Generic.List[object]
foreach ($scope in @('users/0','groups/735051')) {
    $start = 0
    while ($true) {
        curl.exe -s -m 60 "http://localhost:23119/api/$scope/items?format=json&limit=100&start=$start&itemType=attachment" -o zt_att.json
        $items = $null
        try { $items = (Get-Content -Raw -Encoding UTF8 zt_att.json) | ConvertFrom-Json } catch { }
        if (-not $items -or @($items).Count -eq 0) { break }
        foreach ($it in @($items)) {
            if ($it.data.contentType -ne 'application/pdf') { continue }
            $att.Add([pscustomobject]@{ scope=$scope; key=$it.key; parent=[string]$it.data.parentItem; filename=[string]$it.data.filename })
        }
        $start += 100
        if ($start -gt 5000) { break }
    }
}
Write-Output ("pdf attachments: " + $att.Count)

$staging = 'D:\Github\Bipedal_Robot\SADb_audit\pdf_staging'
New-Item -ItemType Directory -Force -Path $staging | Out-Null
$out = New-Object System.Collections.Generic.List[object]
$copied = 0
foreach ($a in $att) {
    if (-not $a.parent -or -not $meta.ContainsKey($a.parent)) { continue }
    $p = $meta[$a.parent]
    $doi = ($p.doi -replace '^https?://(dx\.)?doi\.org/','').ToLower()
    $file = Join-Path (Join-Path $storageRoot $a.key) $a.filename
    $exists = Test-Path -LiteralPath $file
    $yr = ''; if ($p.year -match '(\d{4})') { $yr = $Matches[1] }
    $out.Add([pscustomobject]@{
        library = $a.scope; zotero_parent = $a.parent; att_key = $a.key
        doi = $doi; surname = (Clean-Name $p.surname); year = $yr
        title = ($p.title -replace '\s+',' '); file_exists = $exists; local_path = $file
    })
    if ($exists -and $doi) {
        $dest = Join-Path $staging ((Clean-Name $p.surname) + '_' + $yr + '__' + $a.key + '.pdf')
        if (-not (Test-Path -LiteralPath $dest)) { Copy-Item -LiteralPath $file -Destination $dest -Force; $copied++ }
    }
}
$out | Export-Csv -Path 'D:\Github\Bipedal_Robot\SADb_audit\pdf_inventory.csv' -NoTypeInformation -Encoding UTF8
Write-Output ("inventory rows: " + $out.Count + "  newly copied: " + $copied)
