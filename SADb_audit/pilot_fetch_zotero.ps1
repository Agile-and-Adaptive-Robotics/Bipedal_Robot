# Pilot batch grounding v2: Zotero local API with qmode=everything (q= alone
# does NOT search the DOI field; everything does). Personal first, then AARL
# group 735051. Saves abstractNotes to files, reports children.
$ErrorActionPreference = 'Stop'
$pairs = @(
    @{ k='9KAXASMM'; doi='10.1016/0959-4388(95)80107-3' },
    @{ k='ZMLZC8RV'; doi='10.1113/jphysiol.1910.sp001362' },
    @{ k='H7FU2J6M'; doi='10.1098/rspb.1911.0077' },
    @{ k='C825GTPT'; doi='10.1016/0006-8993(75)90401-1' },
    @{ k='23ZX8RAS'; doi='10.1152/jn.1977.40.4.737' },
    @{ k='RYXH6JSP'; doi='10.1152/jn.1980.44.3.475' },
    @{ k='MZQZRLAD'; doi='10.1113/jphysiol.1981.sp013556' },
    @{ k='8VKY3MVB'; doi='10.1007/bf00204048' },
    @{ k='NUQ2JRWH'; doi='10.1016/s0960-9822(01)00581-4' },
    @{ k='89WM8PFQ'; doi='10.1016/s0966-6362(99)00052-1' }
)
$out = New-Object System.Collections.Generic.List[string]
foreach ($p in $pairs) {
    # search on a DOI substring without parens/slashes issues: use full DOI
    $enc = [uri]::EscapeDataString($p.doi)
    foreach ($scope in @('users/0','groups/735051')) {
        $url = "http://localhost:23119/api/$scope/items?format=json&q=$enc&qmode=everything&limit=25"
        curl.exe -s -m 15 $url -o pilot_tmp.json
        $raw = Get-Content -Raw pilot_tmp.json
        $items = $null
        try { $items = $raw | ConvertFrom-Json } catch { }
        if (-not $items) { $out.Add("$($p.k) $scope ERROR/EMPTY len=$($raw.Length)"); continue }
        $tops = @($items | Where-Object { -not $_.data.parentItem -and $_.data.DOI -and ($_.data.DOI -replace '^https?://(dx\.)?doi\.org/','') -ieq $p.doi })
        if ($tops.Count -eq 0) { $out.Add("$($p.k) $scope NOMATCH"); continue }
        foreach ($t in $tops) {
            $abs = $t.data.abstractNote
            $absLen = 0; if ($abs) { $absLen = $abs.Length }
            $out.Add("$($p.k) $scope ITEM $($t.key) type=$($t.data.itemType) absLen=$absLen")
            if ($abs) { $abs | Set-Content -Encoding UTF8 ("pilot_abs_" + $p.k + "_" + ($scope -replace '/','_') + ".txt") }
            curl.exe -s -m 15 "http://localhost:23119/api/$scope/items/$($t.key)/children?format=json" -o pilot_tmp2.json
            $kids = (Get-Content -Raw pilot_tmp2.json) | ConvertFrom-Json
            if ($kids) {
                foreach ($c in @($kids)) {
                    $ct = $c.data.itemType
                    $extra = ''
                    if ($ct -eq 'attachment') { $extra = ' ct=' + $c.data.contentType + ' link=' + $c.data.linkMode }
                    $out.Add("    child $($c.key) $ct$extra")
                }
            }
        }
    }
    $out.Add('----')
}
$out | Set-Content -Encoding UTF8 pilot_zotero_report.txt
Get-Content pilot_zotero_report.txt
