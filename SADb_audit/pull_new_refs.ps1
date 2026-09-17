# Pull full metadata (title/creators/year/DOI/abstract/journal) for the new items.
$keys = @('WJ93DY7W','JMJIQ9ZZ','6CEVM6YE','6NIJ7VIE','8YCE28D9','4HP4LJDQ','6Z7RBVRC','PXUSZ7NR','ZYA3GGJF','SUEWLVPW','QG6A6BH4','E3U29CXK','R29RM7JI','WGZSJTQH','TQYBQQUG','TBHKUUAK','GM9GET82','FNQCAAU7','4S96Y3FG','YE5VCDMJ','JEJUI2NS','F55LBFZG','7NTJY5GQ','BLL4S52K','X43H35G9','KUYUIBA7','IZ2IA2ML','JEUZYHH8','DSTCLWMI','PV8Q9QMY','6NNY7PSJ','F6HULPJH','7ZJV4GTK')
$out = New-Object System.Collections.Generic.List[object]
foreach ($k in $keys) {
    curl.exe -s -m 30 "http://localhost:23119/api/groups/735051/items/$k?format=json" -o zt_item.json
    $it = $null
    try { $it = (Get-Content -Raw -Encoding UTF8 zt_item.json) | ConvertFrom-Json } catch { continue }
    $d = $it.data
    $sn = ''
    foreach ($c in @($d.creators)) {
        if ($c.lastName) { $sn = [string]$c.lastName; break }
        if ($c.name) { $p = ([string]$c.name).Trim().Split(' '); $sn = $p[$p.Count-1]; break }
    }
    $ab = [string]$d.abstractNote
    $out.Add([pscustomobject]@{
        key = $k; itemType = $d.itemType
        doi = ([string]$d.DOI).Trim()
        title = (($d.title -replace '\s+',' '))
        surname = $sn; year = ([string]$d.date)
        pub = (($d.publicationTitle -replace '\s+',' '))
        url = ([string]$d.url)
        abstract = $ab
        numChildren = $it.meta.numChildren
    })
}
$out | ConvertTo-Json -Depth 3 | Set-Content -Encoding UTF8 'D:\Github\Bipedal_Robot\SADb_audit\new_refs_full.json'
Write-Output ("pulled: " + $out.Count)
foreach ($o in $out) { Write-Output ("  " + $o.key + " " + $o.itemType + " absLen=" + $o.abstract.Length + " kids=" + $o.numChildren + " :: " + $o.title.Substring(0,[Math]::Min(60,$o.title.Length))) }
