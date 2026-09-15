# Dump fresh personal-library author data (key, DOI, title, creator surnames).
$out = New-Object System.Collections.Generic.List[object]
$start = 0
while ($true) {
    curl.exe -s -m 30 "http://localhost:23119/api/users/0/items?format=json&limit=100&start=$start&itemType=-attachment" -o zt_page.json
    $items = $null
    try { $items = (Get-Content -Raw -Encoding UTF8 zt_page.json) | ConvertFrom-Json } catch { }
    if (-not $items -or @($items).Count -eq 0) { break }
    foreach ($it in @($items)) {
        if ($it.data.parentItem) { continue }
        if ($it.data.itemType -eq 'note') { continue }
        $surnames = New-Object System.Collections.Generic.List[string]
        foreach ($c in @($it.data.creators)) {
            if ($c.lastName) { $surnames.Add([string]$c.lastName) }
            elseif ($c.name) {
                $parts = ([string]$c.name).Trim().Split(' ')
                if ($parts.Count -gt 0) { $surnames.Add($parts[$parts.Count-1]) }
            }
        }
        if ($surnames.Count -gt 0) {
            $out.Add([pscustomobject]@{ key = $it.key; doi = [string]$it.data.DOI; title = [string]$it.data.title; surnames = $surnames })
        }
    }
    $start += 100
    if ($start -gt 1500) { break }
}
$json = $out | ConvertTo-Json -Depth 4 -Compress
[System.IO.File]::WriteAllText("$env:TEMP\personal_fresh_authors.json", $json)
Write-Output ("dumped items: " + $out.Count)
