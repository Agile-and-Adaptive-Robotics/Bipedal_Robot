# Batch 2 grounding round 2: PubMed for Chu 2018 + Dubuc 2008
$ErrorActionPreference = 'Continue'
function Get-PubMedAbstract($tag, $doi) {
    $term = [uri]::EscapeDataString(($doi + '[DOI]'))
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=$term&retmode=json" -o batch2_pm_s.json
    Start-Sleep -Milliseconds 500
    $j = (Get-Content -Raw batch2_pm_s.json) | ConvertFrom-Json
    $ids = @($j.esearchresult.idlist)
    if ($ids.Count -eq 0) { Write-Output "$tag NO PMID"; return }
    $pmid = $ids[0]
    Write-Output "$tag PMID=$pmid"
    Start-Sleep -Milliseconds 500
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$pmid&rettype=abstract&retmode=text" -o ("batch2_ground_" + $tag + ".txt")
    Write-Output ("  saved " + (Get-Item ("batch2_ground_" + $tag + ".txt")).Length + " bytes")
}
Get-PubMedAbstract 'chu2018' '10.1152/jn.00486.2017'
Get-PubMedAbstract 'dubuc2008' '10.1016/j.brainresrev.2007.07.016'
