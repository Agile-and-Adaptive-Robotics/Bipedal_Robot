# List items added since 2026-09-15 20:00Z in personal + AARL group libraries.
$out = New-Object System.Collections.Generic.List[object]
foreach ($scope in @('users/0','groups/735051')) {
    $start = 0
    while ($true) {
        curl.exe -s -m 60 "http://localhost:23119/api/$scope/items?format=json&limit=100&start=$start&itemType=-attachment&sort=dateAdded&direction=desc" -o zt_new.json
        $items = $null
        try { $items = (Get-Content -Raw -Encoding UTF8 zt_new.json) | ConvertFrom-Json } catch { }
        if (-not $items -or @($items).Count -eq 0) { break }
        foreach ($it in @($items)) {
            if ($it.data.parentItem) { continue }
            if ($it.data.dateAdded -lt '2026-09-15T20:00:00') { break }
            $sn = ''
            foreach ($c in @($it.data.creators)) {
                if ($c.lastName) { $sn = [string]$c.lastName; break }
                if ($c.name) { $p = ([string]$c.name).Trim().Split(' '); $sn = $p[$p.Count-1]; break }
            }
            $out.Add([pscustomobject]@{
                scope = $scope; key = $it.key; added = $it.data.dateAdded
                itemType = $it.data.itemType
                doi = ([string]$it.data.DOI).Trim()
                surname = $sn; year = ([string]$it.data.date)
                title = (($it.data.title -replace '\s+',' '))
                abstractLen = ([string]$it.data.abstractNote).Length
            })
        }
        $start += 100
        if ($start -gt 2000) { break }
    }
}
$out | Export-Csv -Path 'D:\Github\Bipedal_Robot\SADb_audit\new_zotero_refs.csv' -NoTypeInformation -Encoding UTF8
Write-Output ("new items since 2026-09-15T20:00Z: " + $out.Count)
foreach ($o in $out) { Write-Output ("  " + $o.scope + " " + $o.key + " " + $o.itemType + " doi=[" + $o.doi + "] " + $o.surname + " " + $o.year + " :: " + $o.title.Substring(0, [Math]::Min(70, $o.title.Length))) }
