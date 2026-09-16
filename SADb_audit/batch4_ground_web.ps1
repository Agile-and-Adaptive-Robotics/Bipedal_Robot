# Batch 4 grounding round 2: PubMed abstracts (esearch by DOI -> efetch)
$ErrorActionPreference = 'Continue'
function Get-PubMedAbstract($tag, $doi) {
    $term = [uri]::EscapeDataString(($doi + '[DOI]'))
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=$term&retmode=json" -o batch4_pm_s.json
    Start-Sleep -Milliseconds 600
    $j = (Get-Content -Raw batch4_pm_s.json) | ConvertFrom-Json
    $ids = @($j.esearchresult.idlist)
    if ($ids.Count -eq 0) { Write-Output "$tag NO PMID"; return }
    $pmid = $ids[0]
    Write-Output "$tag PMID=$pmid"
    Start-Sleep -Milliseconds 600
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$pmid&rettype=abstract&retmode=text" -o ("batch4_ground_" + $tag + ".txt")
    Write-Output ("  saved " + (Get-Item ("batch4_ground_" + $tag + ".txt")).Length + " bytes")
}
Get-PubMedAbstract 'barbeau1987' '10.1016/0006-8993(87)91442-9'
Get-PubMedAbstract 'bassler1998' '10.1016/s0165-0173(98)00006-x'
Get-PubMedAbstract 'armstrong1986' '10.1016/0301-0082(86)90021-3'
Get-PubMedAbstract 'ijspeert2008' '10.1016/j.neunet.2008.03.014'
Get-PubMedAbstract 'quevedo2005' '10.1152/jn.00175.2005'
