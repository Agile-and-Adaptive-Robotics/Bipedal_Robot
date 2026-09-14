# Batch 2 grounding: Zotero personal library per DOI (qmode=everything).
$ErrorActionPreference = 'Stop'
$pairs = @(
    @{ k='6J8X9RQC'; doi='10.1038/s41467-023-36587-w' },
    @{ k='MKE6485E'; doi='10.1152/jn.00486.2017' },
    @{ k='Q7VR4INI'; doi='10.1523/jneurosci.1654-20.2020' },
    @{ k='7F8N78A9'; doi='10.1016/j.arcontrol.2019.04.004' },
    @{ k='CEQX3Y93'; doi='10.1152/jn.1997.77.6.3311' },
    @{ k='WLGYZ4MX'; doi='10.1016/j.brainresrev.2007.07.016' },
    @{ k='FXXKJEFZ'; doi='10.1111/j.1469-7793.2000.00639.x' },
    @{ k='Y3SF4GJC'; doi='10.1088/1748-3190/aa8290' },
    @{ k='KIG9DEKJ'; doi='10.7554/elife.73424' },
    @{ k='5BM93A93'; doi='10.1007/978-1-4757-0964-3' }
)
$out = New-Object System.Collections.Generic.List[string]
foreach ($p in $pairs) {
    $enc = [uri]::EscapeDataString($p.doi)
    $url = "http://localhost:23119/api/users/0/items?format=json&q=$enc&qmode=everything&limit=25"
    curl.exe -s -m 15 $url -o batch2_tmp.json
    $raw = Get-Content -Raw batch2_tmp.json
    $items = $null
    try { $items = $raw | ConvertFrom-Json } catch { }
    if (-not $items) { $out.Add("$($p.k) ERROR/EMPTY len=$($raw.Length)"); $out.Add('----'); continue }
    $tops = @($items | Where-Object { -not $_.data.parentItem -and $_.data.DOI -and ($_.data.DOI -replace '^https?://(dx\.)?doi\.org/','') -ieq $p.doi })
    if ($tops.Count -eq 0) { $out.Add("$($p.k) NOMATCH"); $out.Add('----'); continue }
    foreach ($t in $tops) {
        $abs = $t.data.abstractNote
        $absLen = 0; if ($abs) { $absLen = $abs.Length }
        $out.Add("$($p.k) ITEM $($t.key) type=$($t.data.itemType) absLen=$absLen")
        if ($abs) { $abs | Set-Content -Encoding UTF8 ("batch2_abs_" + $p.k + ".txt") }
        curl.exe -s -m 15 "http://localhost:23119/api/users/0/items/$($t.key)/children?format=json" -o batch2_tmp2.json
        $kids = (Get-Content -Raw batch2_tmp2.json) | ConvertFrom-Json
        if ($kids) {
            foreach ($c in @($kids)) {
                $extra = ''
                if ($c.data.itemType -eq 'attachment') { $extra = ' ct=' + $c.data.contentType + ' link=' + $c.data.linkMode }
                $out.Add("    child $($c.key) $($c.data.itemType)$extra")
            }
        }
    }
    $out.Add('----')
}
$out | Set-Content -Encoding UTF8 batch2_zotero_report.txt
Get-Content batch2_zotero_report.txt
