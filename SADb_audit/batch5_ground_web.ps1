# Batch 5 round 2: PubMed for Buford1993 + Hurteau2015; Crossref year check for Rubeo.
$ErrorActionPreference = 'Continue'
function Get-PubMedAbstract($tag, $doi) {
    $term = [uri]::EscapeDataString(($doi + '[DOI]'))
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=$term&retmode=json" -o batch5_pm_s.json
    Start-Sleep -Milliseconds 600
    $j = (Get-Content -Raw batch5_pm_s.json) | ConvertFrom-Json
    $ids = @($j.esearchresult.idlist)
    if ($ids.Count -eq 0) { Write-Output "$tag NO PMID"; return }
    $pmid = $ids[0]
    Write-Output "$tag PMID=$pmid"
    Start-Sleep -Milliseconds 600
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$pmid&rettype=abstract&retmode=text" -o ("batch5_ground_" + $tag + ".txt")
    Write-Output ("  saved " + (Get-Item ("batch5_ground_" + $tag + ".txt")).Length + " bytes")
}
Get-PubMedAbstract 'buford1993' '10.1152/jn.1993.70.3.1102'
Get-PubMedAbstract 'hurteau2015' '10.1152/jn.00739.2014'
curl.exe -s -m 20 "https://api.crossref.org/works/10.3390/app8010006" -o batch5_cr_rubeo.json
$j2 = (Get-Content -Raw batch5_cr_rubeo.json) | ConvertFrom-Json
Write-Output ("rubeo crossref year: " + $j2.message.issued.'date-parts'[0][0] + " ; published: " + $j2.message.published.'date-parts'[0][0] + "-" + $j2.message.published.'date-parts'[0][1] + " ; container: " + $j2.message.'container-title'[0])
