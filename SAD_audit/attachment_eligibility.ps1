$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
function Read-CsvUtf8([string]$path) { return (Get-Content -LiteralPath $path -Encoding UTF8 | ConvertFrom-Csv) }
function Norm-Loose([string]$s) {
  if ([string]::IsNullOrWhiteSpace($s)) { return "" }
  $t = ([regex]::Replace($s, '\s+', ' ')).Trim().ToLower()
  $t = [regex]::Replace($t, '[^a-z0-9]+', ' ')
  return ([regex]::Replace($t, '\s+', ' ')).Trim()
}

# ---- 0. Patch airtable_papers_slim.csv to post-fix state + add DOI column ----
$at = Read-CsvUtf8 "$out\airtable_papers_slim.csv"
$bf = Read-CsvUtf8 "$out\airtable_doi_backfill.csv"
$doiOf = @{}
foreach ($b in $bf) { if ($b.DOI -and $b.DOI.Trim()) { $doiOf[$b.id] = $b.DOI.Trim() } }
$doiOf['recLQre1rXrmJOLU3'] = '10.1113/jphysiol.1971.sp009487'
$doiOf['recdSXsgiaZQQGIQ5'] = '10.1016/j.asd.2015.07.001'
$doiOf['rectZLWJhnmJpn1le'] = '10.1016/s1388-2457(03)00120-2'
$titleFix = @{
  'rec8WdUONYgfvw1zw' = 'A Role for Hip Position in Initiating the Swing-to-Stance Transition in Walking Cats'
  'recHUhzNl3CEaF8pE' = 'Leg Coordination Mechanisms in the Stick Insect Applied to Hexapod Robot Locomotion'
  'recPoyqO9O9txmwm5' = 'Speed dependency in α-motoneuron activity and locomotor modules in human locomotion: indirect evidence for phylogenetically conserved spinal circuits'
  'recfiz9tzgn955QdA' = 'Contribution of Hind Limb Flexor Muscle Afferents to the Timing of Phase Transitions in the Cat Step Cycle'
  'reclwbsReGhWM71WF' = 'Bio-inspired controller achieving forward speed modulation with a 3D bipedal walker'
  'recxGoi3bAS1glsgD' = 'Relative Contribution of Proprioceptive and Vestibular Sensory Systems to Locomotion: Opportunities for Discovery in the Age of Molecular Science'
  'rectZLWJhnmJpn1le' = 'Spinal Cord Pattern Generators for Locomotion'
}
$yearFix = @{ 'rec3237UGETrSD0Uw'='2003'; 'recGXAkc9YqiOxZyT'='2022'; 'recIKY3z8nHIQBBVJ'='2018'; 'recV8nQBbNyZDolul'='2017'; 'recjEWZKD03OC2BYS'='2016'; 'recoT9VR7NIatzV00'='2016'; 'recyRnIe8WN35qmqC'='1995' }
$authFix = @{ 'rec2jytMTqUIEJMaK'='Robertson and Stein'; 'recFcuIXzjebeAqZF'='Hurteau'; 'recGXAkc9YqiOxZyT'='Shevtsova'; 'recQkgW6mi6AezV1b'='Latash'; 'recjEWZKD03OC2BYS'='Shevtsova'; 'recDKSwBYMQO5uNiK'='Shevtsova'; 'recfCLJpSnraphHug'='McCrea and Rybak' }
$patched = foreach ($a in $at) {
  if ($titleFix.ContainsKey($a.id)) { $a.title = $titleFix[$a.id] }
  if ($yearFix.ContainsKey($a.id))  { $a.year = $yearFix[$a.id] }
  if ($authFix.ContainsKey($a.id))  { $a.author = $authFix[$a.id] }
  $doi = ''; if ($doiOf.ContainsKey($a.id)) { $doi = $doiOf[$a.id] }
  [pscustomobject]@{ id = $a.id; title = $a.title; author = $a.author; year = $a.year; animals = $a.animals; DOI = $doi }
}
$patched | Export-Csv "$out\airtable_papers_slim.csv" -NoTypeInformation -Encoding utf8
Write-Output ("slim CSV patched: " + @($patched).Count + " rows, " + @(($patched | Where-Object { $_.DOI })).Count + " with DOI")

# ---- 1. Dissertation cited set ----
$disRoot = "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
$texs = @(Get-ChildItem -LiteralPath $disRoot -Recurse -Include *.tex -File -ErrorAction SilentlyContinue)
$bibs = @(Get-ChildItem -LiteralPath $disRoot -Recurse -Include *.bib -File -ErrorAction SilentlyContinue)
Write-Output ("dissertation files: " + $texs.Count + " .tex, " + $bibs.Count + " .bib under " + $disRoot)
$citeKeys = New-Object System.Collections.Generic.HashSet[string]
foreach ($t in $texs) {
  $raw = Get-Content -LiteralPath $t.FullName -Raw -Encoding UTF8
  foreach ($m in [regex]::Matches($raw, '\\(?:cite|citep|citet|citeauthor|citeyear|citealt|citealp|autocite|parencite|textcite|footcite|footcitetext)[a-z]*\*?(?:\[[^\]]*\])*\{([^}]*)\}')) {
    foreach ($k in ($m.Groups[1].Value -split ',')) { $k2 = $k.Trim(); if ($k2) { [void]$citeKeys.Add($k2) } }
  }
}
Write-Output ("unique cited keys: " + $citeKeys.Count)
# bib key -> title/doi
$bibMap = @{}
foreach ($b in $bibs) {
  $raw = Get-Content -LiteralPath $b.FullName -Raw -Encoding UTF8
  foreach ($m in [regex]::Matches($raw, '@\w+\s*\{\s*([^,\s]+)\s*,(.*?)\n\}', [System.Text.RegularExpressions.RegexOptions]::Singleline)) {
    $key = $m.Groups[1].Value; $body = $m.Groups[2].Value
    $tm = [regex]::Match($body, '(?i)^\s*title\s*=\s*"{0,1}(.+?)"{0,1}\s*,\s*$', [System.Text.RegularExpressions.RegexOptions]::Multiline)
    if (-not $tm.Success) { $tm = [regex]::Match($body, '(?i)title\s*=\s*[{"](.+?)[}"]\s*,') }
    $dm = [regex]::Match($body, '(?i)doi\s*=\s*[{"]([^}"]+)[}"]')
    if (-not $bibMap.ContainsKey($key)) {
      $bibMap[$key] = @{ title = $tm.Groups[1].Value; doi = $dm.Groups[1].Value }
    }
  }
}
Write-Output ("bib entries parsed: " + $bibMap.Count)
$citedDois = New-Object System.Collections.Generic.HashSet[string]
$citedTitles = @{}
$unresolvedKeys = New-Object System.Collections.ArrayList
foreach ($k in $citeKeys) {
  if ($bibMap.ContainsKey($k)) {
    $d = $bibMap[$k].doi; if ($d) { [void]$citedDois.Add($d.Trim().ToLower().TrimStart('https://doi.org/')) }
    $t = $bibMap[$k].title; if ($t) { $citedTitles[(Norm-Loose $t)] = $k }
  } else { [void]$unresolvedKeys.Add($k) }
}
Write-Output ("cited DOIs: " + $citedDois.Count + "  cited titles: " + $citedTitles.Count + "  unresolved keys: " + $unresolvedKeys.Count)
if ($unresolvedKeys.Count -gt 0) { Write-Output ("  unresolved: " + ($unresolvedKeys -join ', ')) }

# ---- 2. Personal attachments via local API (curl.exe; Invoke-RestMethod fails here) ----
$attList = New-Object System.Collections.ArrayList
$start = 0
while ($true) {
  $tmp = "$env:TEMP\zot_att_$start.json"
  $url = "http://localhost:23119/api/users/0/items?format=json&itemType=attachment&limit=100&start=$start"
  curl.exe -s -o $tmp $url
  $page = Get-Content -LiteralPath $tmp -Raw -Encoding UTF8 | ConvertFrom-Json
  $page = @($page)
  foreach ($i in $page) {
    [void]$attList.Add([pscustomobject]@{
      key = $i.key; parent = $i.data.parentItem; linkMode = $i.data.linkMode
      filename = $i.data.filename; contentType = $i.data.contentType
    })
  }
  if ($page.Count -lt 100) { break }
  $start += 100
}
Write-Output ("personal attachments fetched: " + $attList.Count)
$attList | Export-Csv "$out\personal_attachments.csv" -NoTypeInformation -Encoding utf8

# ---- 3. Eligibility ----
$per = Read-CsvUtf8 "$out\personal_SAD_items_slim.csv"
$perBy = @{}
foreach ($p in $per) { $perBy[$p.key] = $p }
$atDoiSet = @{}; $atTitleSet = @{}
foreach ($a in $patched) {
  if ($a.DOI) { $atDoiSet[$a.DOI.Trim().ToLower()] = 1 }
  $lt = Norm-Loose $a.title; if ($lt) { $atTitleSet[$lt] = 1 }
}
$inAirtable = @{}
foreach ($p in $per) {
  $pd = ''; if ($p.DOI) { $pd = $p.DOI.Trim().ToLower() }
  $lt = Norm-Loose $p.title
  if (($pd -and $atDoiSet.ContainsKey($pd)) -or ($lt -and $atTitleSet.ContainsKey($lt))) { $inAirtable[$p.key] = $true }
}
Write-Output ("personal SAD items also on Airtable: " + $inAirtable.Count + " / " + @($per).Count)

$fileAtt = @($attList | Where-Object { $_.linkMode -eq 'imported_file' -or $_.linkMode -eq 'imported_url' })
Write-Output ("file attachments (imported_file/imported_url): " + $fileAtt.Count)

$eligible = New-Object System.Collections.ArrayList
$protected = 0
foreach ($att in $fileAtt) {
  if (-not $att.parent -or -not $perBy.ContainsKey($att.parent)) { continue }
  $p = $perBy[$att.parent]
  $pd = ''; if ($p.DOI) { $pd = $p.DOI.Trim().ToLower() }
  $lt = Norm-Loose $p.title
  $cited = ($pd -and $citedDois.Contains($pd)) -or ($lt -and $citedTitles.ContainsKey($lt))
  if ($cited) { $protected++; continue }
  if ($inAirtable.ContainsKey($p.key)) {
    [void]$eligible.Add([pscustomobject]@{
      parentKey = $p.key; title = $p.title; DOI = $p.DOI
      attKey = $att.key; filename = $att.filename; linkMode = $att.linkMode
    })
  }
}
$eligible | Export-Csv "$out\personal_SAD_attachment_purge_eligible.csv" -NoTypeInformation -Encoding utf8
$parents = @($eligible | Select-Object -ExpandProperty parentKey -Unique)
Write-Output ("ELIGIBLE for attachment removal: " + $eligible.Count + " attachments on " + $parents.Count + " items (dissertation-protected with attachments: " + $protected + ")")
Write-Output "DONE"
