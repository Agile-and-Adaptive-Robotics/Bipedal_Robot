# Check Zotero for Armstrong 1986 + Bassler 1998 (fresh connector saves?)
$ErrorActionPreference = 'Continue'
$targets = @(
    @{ k='VRUAGILB'; doi='10.1016/0301-0082(86)90021-3' },
    @{ k='EFRSBRUQ'; doi='10.1016/s0165-0173(98)00006-x' }
)
foreach ($t in $targets) {
    $enc = [uri]::EscapeDataString($t.doi)
    foreach ($scope in @('users/0','groups/735051')) {
        curl.exe -s -m 15 "http://localhost:23119/api/$scope/items?format=json&q=$enc&qmode=everything&limit=25" -o b4_tmp.json
        $items = $null
        try { $items = (Get-Content -Raw b4_tmp.json) | ConvertFrom-Json } catch { }
        if (-not $items) { Write-Output "$($t.k) $scope EMPTY"; continue }
        $tops = @($items | Where-Object { -not $_.data.parentItem -and $_.data.DOI -and ($_.data.DOI -replace '^https?://(dx\.)?doi\.org/','') -ieq $t.doi })
        foreach ($top in $tops) {
            Write-Output "$($t.k) $scope ITEM $($top.key) absLen=$(([string]$top.data.abstractNote).Length)"
            curl.exe -s -m 15 "http://localhost:23119/api/$scope/items/$($top.key)/children?format=json" -o b4_kids.json
            $kids = (Get-Content -Raw b4_kids.json) | ConvertFrom-Json
            foreach ($c in @($kids)) {
                Write-Output "   child $($c.key) $($c.data.itemType) $($c.data.contentType)"
            }
        }
    }
}
