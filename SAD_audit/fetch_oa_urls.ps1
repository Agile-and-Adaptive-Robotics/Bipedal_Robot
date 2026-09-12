$ErrorActionPreference = "Continue"
$out = "D:\Github\Bipedal_Robot\SAD_audit"
$create = Import-Csv "$out\airtable_create50.csv" -Encoding UTF8
$ids    = Import-Csv "$out\airtable_created50_ids.csv" -Encoding UTF8
$idByZ  = @{}; foreach ($i in $ids) { $idByZ[$i.zotero_key] = $i.airtable_id }

$rows = New-Object System.Collections.ArrayList
foreach ($r in $create) {
  $atId = $idByZ[$r.zotero_key]
  $pdf = ''
  if ($r.DOI) {
    $doi = [uri]::EscapeDataString($r.DOI)
    curl.exe -s -o "$env:TEMP\oa_$($r.zotero_key).json" --max-time 30 "https://api.openalex.org/works/doi:$doi"
    $raw = Get-Content "$env:TEMP\oa_$($r.zotero_key).json" -Raw -Encoding UTF8
    $ms = [regex]::Matches($raw, '"pdf_url"\s*:\s*"([^"]+)"')
    if ($ms.Count -gt 0) {
      $cands = @($ms | ForEach-Object { $_.Groups[1].Value.Replace('\/', '/') } | Where-Object { $_ -match '\.pdf' })
      if ($cands.Count -eq 0) { $cands = @($ms | ForEach-Object { $_.Groups[1].Value.Replace('\/', '/') }) }
      # preferPMC/core over generic
      $pdf = @($cands | Where-Object { $_ -match 'pmc\.ncbi|pmcid' } | Select-Object -First 1)
      if (-not $pdf) { $pdf = $cands[0] }
    }
  }
  [void]$rows.Add([pscustomobject]@{ zotero_key = $r.zotero_key; airtable_id = $atId; title = $r.title; DOI = $r.DOI; pdf_url = $pdf })
  Write-Output ($r.zotero_key + " -> " + ($(if ($pdf) { $pdf.Substring(0, [Math]::Min(90, $pdf.Length)) } else { '(no OA pdf)' })))
}
$rows | Export-Csv "$out\airtable_create50_oa_urls.csv" -NoTypeInformation -Encoding utf8
$have = @($rows | Where-Object { $_.pdf_url })
Write-Output ("OA pdf urls found: " + $have.Count + " / 50")
Write-Output "DONE"
