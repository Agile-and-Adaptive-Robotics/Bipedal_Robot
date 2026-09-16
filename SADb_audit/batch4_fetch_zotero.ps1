# Batch 4 grounding: Zotero personal library per DOI (qmode=everything).
$ErrorActionPreference = 'Stop'
$pairs = @(
    @{ k='STI84977'; doi='10.1007/bf00230939' },
    @{ k='YSBJLGU7'; doi='10.1016/0006-8993(87)91442-9' },
    @{ k='BI2N2VXZ'; doi='10.1016/j.neuron.2014.02.013' },
    @{ k='FT6FLCCE'; doi='10.1007/978-1-4614-7320-6_49-2' },
    @{ k='CV83DAZY'; doi='10.1113/jphysiol.1957.sp005794' },
    @{ k='EFRSBRUQ'; doi='10.1016/s0165-0173(98)00006-x' },
    @{ k='VRUAGILB'; doi='10.1016/0301-0082(86)90021-3' },
    @{ k='28XSSMP9'; doi='10.1016/j.neunet.2008.03.014' },
    @{ k='UVG6PMDD'; doi='10.1152/jn.00175.2005' },
    @{ k='FDYA43W5'; doi='10.1007/978-3-319-63537-8_15' }
)
$out = New-Object System.Collections.Generic.List[string]
foreach ($p in $pairs) {
    $enc = [uri]::EscapeDataString($p.doi)
    $url = "http://localhost:23119/api/users/0/items?format=json&q=$enc&qmode=everything&limit=25"
    curl.exe -s -m 15 $url -o batch4_tmp.json
    $raw = Get-Content -Raw batch4_tmp.json
    $items = $null
    try { $items = $raw | ConvertFrom-Json } catch { }
    if (-not $items) { $out.Add("$($p.k) ERROR/EMPTY len=$($raw.Length)"); $out.Add('----'); continue }
    $tops = @($items | Where-Object { -not $_.data.parentItem -and $_.data.DOI -and ($_.data.DOI -replace '^https?://(dx\.)?doi\.org/','') -ieq $p.doi })
    if ($tops.Count -eq 0) { $out.Add("$($p.k) NOMATCH"); $out.Add('----'); continue }
    foreach ($t in $tops) {
        $abs = $t.data.abstractNote
        $absLen = 0; if ($abs) { $absLen = $abs.Length }
        $out.Add("$($p.k) ITEM $($t.key) type=$($t.data.itemType) absLen=$absLen")
        if ($abs) { $abs | Set-Content -Encoding UTF8 ("batch4_abs_" + $p.k + ".txt") }
        curl.exe -s -m 15 "http://localhost:23119/api/users/0/items/$($t.key)/children?format=json" -o batch4_tmp2.json
        $kids = (Get-Content -Raw batch4_tmp2.json) | ConvertFrom-Json
        if ($kids) {
            foreach ($c in @($kids)) {
                $extra = ''
                if ($c.data.itemType -eq 'attachment') { $extra = ' ct=' + $c.data.contentType }
                $out.Add("    child $($c.key) $($c.data.itemType)$extra")
            }
        }
    }
    $out.Add('----')
}
$out | Set-Content -Encoding UTF8 batch4_zotero_report.txt
Get-Content batch4_zotero_report.txt
