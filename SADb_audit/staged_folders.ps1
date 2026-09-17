# Map every staged PDF to its AARL/personal Zotero folder(s).
$ErrorActionPreference = 'Continue'

function Get-Surname($creators) {
    foreach ($c in @($creators)) {
        if ($c.lastName) { return [string]$c.lastName }
        if ($c.name) { $p = ([string]$c.name).Trim().Split(' '); return $p[$p.Count-1] }
    }
    return ''
}

$collName = @{}
$gc = Get-Content -Raw -Encoding UTF8 'D:\Github\Bipedal_Robot\SADb_audit\group_collections.json' | ConvertFrom-Json
foreach ($c in @($gc)) { $collName[$c.key] = [string]$c.data.name }

$map = @{}
foreach ($scope in @('groups/735051','users/0')) {
    $start = 0
    while ($true) {
        curl.exe -s -m 60 "http://localhost:23119/api/$scope/items?format=json&limit=100&start=$start&itemType=-attachment" -o zt_coll.json
        $items = $null
        try { $items = (Get-Content -Raw -Encoding UTF8 zt_coll.json) | ConvertFrom-Json } catch { }
        if (-not $items -or @($items).Count -eq 0) { break }
        foreach ($it in @($items)) {
            if ($it.data.parentItem) { continue }
            if (-not $map.ContainsKey($it.key)) {
                $sn = Get-Surname $it.data.creators
                $yr = ''; if (([string]$it.data.date) -match '(\d{4})') { $yr = $Matches[1] }
                $cols = @($it.data.collections | ForEach-Object { [string]$_ }) -join ';'
                $map[$it.key] = [pscustomobject]@{
                    doi = (([string]$it.data.DOI) -replace '^https?://(dx\.)?doi\.org/','').ToLower()
                    surname = $sn; year = $yr
                    title = (($it.data.title -replace '\s+',' '))
                    collections = $cols
                }
            }
        }
        $start += 100
        if ($start -gt 5000) { break }
    }
    Write-Output ("$scope mapped cumulative: " + $map.Count)
}

$inv = Import-Csv 'D:\Github\Bipedal_Robot\SADb_audit\pdf_inventory.csv'
$rows = New-Object System.Collections.Generic.List[object]
foreach ($a in $inv) {
    if (-not $map.ContainsKey($a.zotero_parent)) { continue }
    $m = $map[$a.zotero_parent]
    $names = @()
    foreach ($ck in ($m.collections -split ';' | Where-Object { $_ })) {
        if ($collName.ContainsKey($ck)) { $names += $collName[$ck] }
    }
    $rows.Add([pscustomobject]@{
        att_key = $a.att_key; doi = $m.doi; surname = $m.surname; year = $m.year
        title = $m.title; folders = ($names -join ';'); library = $a.library
    })
}
$rows | Export-Csv -Path 'D:\Github\Bipedal_Robot\SADb_audit\staged_with_folders.csv' -NoTypeInformation -Encoding UTF8
Write-Output ("rows with folders: " + $rows.Count)

# summary per folder for staged PDFs
$rows | ForEach-Object { $_.folders -split ';' } | Where-Object { $_ } |
    Group-Object | Sort-Object Count -Descending |
    ForEach-Object { '{0,6}  {1}' -f $_.Count, $_.Name }
