# Fetch remaining grounding texts: PubMed (eutils) for modern papers,
# Europe PMC full text for the two open-access classics.
$ErrorActionPreference = 'Continue'

function Get-PubMedAbstract($tag, $doi) {
    $term = [uri]::EscapeDataString(($doi + '[DOI]'))
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=$term&retmode=json" -o pilot_pm_s.json
    $j = (Get-Content -Raw pilot_pm_s.json) | ConvertFrom-Json
    $ids = @($j.esearchresult.idlist)
    if ($ids.Count -eq 0) { Write-Output "$tag NO PMID"; return }
    $pmid = $ids[0]
    Write-Output "$tag PMID=$pmid"
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$pmid&rettype=abstract&retmode=text" -o ("pilot_ground_" + $tag + ".txt")
    Write-Output ("  saved " + (Get-Item ("pilot_ground_" + $tag + ".txt")).Length + " bytes")
}

Get-PubMedAbstract 'pearson1995' '10.1016/0959-4388(95)80107-3'
Get-PubMedAbstract 'grillner1975' '10.1016/0006-8993(75)90401-1'
Get-PubMedAbstract 'dietz2000' '10.1016/s0966-6362(99)00052-1'
# McCrea 1980: PMID known from the J Neurophysiol page captured in Zotero
curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=7441311&rettype=abstract&retmode=text" -o pilot_ground_mccrea1980.txt
Write-Output ("mccrea1980 saved " + (Get-Item pilot_ground_mccrea1980.txt).Length + " bytes")

# Europe PMC: classics full text
function Get-Epmc($tag, $doi) {
    $q = [uri]::EscapeDataString(('DOI:"' + $doi + '"'))
    curl.exe -s -m 20 "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=$q&format=json" -o pilot_epmc_s.json
    $j = (Get-Content -Raw pilot_epmc_s.json) | ConvertFrom-Json
    $hits = @($j.resultList.result)
    if ($hits.Count -eq 0) { Write-Output "$tag EPMC NOHIT"; return }
    $h = $hits[0]
    Write-Output ("$tag EPMC hit pmcid=" + $h.pmcid + " pmid=" + $h.pmid + " isOpenAccess=" + $h.isOpenAccess)
    if ($h.pmcid) {
        curl.exe -s -m 30 ("https://www.ebi.ac.uk/europepmc/webservices/rest/" + $h.pmcid + "/fullTextXML") -o ("pilot_ground_" + $tag + "_full.xml")
        Write-Output ("  saved " + (Get-Item ("pilot_ground_" + $tag + "_full.xml")).Length + " bytes")
    }
}
Get-Epmc 'sherrington1910' '10.1113/jphysiol.1910.sp001362'
Get-Epmc 'brown1911' '10.1098/rspb.1911.0077'
