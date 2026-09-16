# Batch 5 grounding: Zotero personal library per DOI (qmode=everything).
$ErrorActionPreference = 'Stop'
$pairs = @(
    @{ k='QKVU6KK3'; doi='10.1152/jn.1993.70.3.1102' },
    @{ k='YX6IC7SQ'; doi='' },
    @{ k='N8C5DJV6'; doi='10.1126/science.1210617' },
    @{ k='723RUJG9'; doi='10.3389/fnbot.2017.00037' },
    @{ k='GKEDLR7B'; doi='10.1111/nyas.12055' },
    @{ k='QRNCP58H'; doi='' },
    @{ k='F4V9VS9X'; doi='10.1162/neco.1992.4.3.356' },
    @{ k='9K4DAP7A'; doi='10.3389/fncir.2023.1146449' },
    @{ k='LYPEYIAI'; doi='10.3390/app8010006' },
    @{ k='4JE58YY8'; doi='10.1152/jn.00739.2014' }
)
$out = New-Object System.Collections.Generic.List[string]
foreach ($p in $pairs) {
    if ($p.doi -ne '') {
        $enc = [uri]::EscapeDataString($p.doi)
        $url = "http://localhost:23119/api/users/0/items?format=json&q=$enc&qmode=everything&limit=25"
        curl.exe -s -m 15 $url -o batch5_tmp.json
        $raw = Get-Content -Raw batch5_tmp.json
        $items = $null
        try { $items = $raw | ConvertFrom-Json } catch { }
        if ($items) {
            $tops = @($items | Where-Object { -not $_.data.parentItem -and $_.data.DOI -and ($_.data.DOI -replace '^https?://(dx\.)?doi\.org/','') -ieq $p.doi })
            if ($tops.Count -eq 0) { $out.Add("$($p.k) NOMATCH") } else {
                foreach ($t in $tops) {
                    $abs = $t.data.abstractNote
                    $absLen = 0; if ($abs) { $absLen = $abs.Length }
                    $out.Add("$($p.k) ITEM $($t.key) type=$($t.data.itemType) absLen=$absLen")
                    if ($abs) { $abs | Set-Content -Encoding UTF8 ("batch5_abs_" + $p.k + ".txt") }
                }
            }
            $out.Add('----'); continue
        }
    }
    # DOI-less (or failed): search by title keyword
    $kw = switch ($p.k) {
        'YX6IC7SQ' { 'Adaptive Control Responses to Behavioral Perturbation' }
        'QRNCP58H' { 'Interlimb communication during human walking' }
        default { '' }
    }
    if ($kw -eq '') { $out.Add("$($p.k) SKIP"); $out.Add('----'); continue }
    $enc2 = [uri]::EscapeDataString($kw)
    curl.exe -s -m 15 "http://localhost:23119/api/users/0/items?format=json&q=$enc2&limit=25" -o batch5_tmp.json
    $raw2 = Get-Content -Raw batch5_tmp.json
    $items2 = $null
    try { $items2 = $raw2 | ConvertFrom-Json } catch { }
    if (-not $items2) { $out.Add("$($p.k) TITLE-NOMATCH"); $out.Add('----'); continue }
    $tops2 = @($items2 | Where-Object { -not $_.data.parentItem -and $_.data.title -and $_.data.title.ToLower().Contains($kw.ToLower().Substring(0,[Math]::Min(30,$kw.Length))) })
    if ($tops2.Count -eq 0) { $out.Add("$($p.k) TITLE-NOMATCH"); $out.Add('----'); continue }
    foreach ($t in $tops2) {
        $abs = $t.data.abstractNote
        $absLen = 0; if ($abs) { $absLen = $abs.Length }
        $out.Add("$($p.k) ITEM $($t.key) type=$($t.data.itemType) absLen=$absLen doi=[$($t.data.DOI)]")
        if ($abs) { $abs | Set-Content -Encoding UTF8 ("batch5_abs_" + $p.k + ".txt") }
    }
    $out.Add('----')
}
$out | Set-Content -Encoding UTF8 batch5_zotero_report.txt
Get-Content batch5_zotero_report.txt
