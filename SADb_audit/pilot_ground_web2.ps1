# Grounding round 2: PMC full text for the classics, Semantic Scholar abstracts
# for Grillner/Zangger 1975 and McCrea 1980. Pauses to respect eutils rate limits.
$ErrorActionPreference = 'Continue'

function Pause-Short { Start-Sleep -Milliseconds 700 }

# --- PMC search for the classics ---
$targets = @(
    @{ tag='sherrington1910'; doi='10.1113/jphysiol.1910.sp001362' },
    @{ tag='brown1911'; doi='10.1098/rspb.1911.0077' }
)
foreach ($t in $targets) {
    $term = [uri]::EscapeDataString(($t.doi + '[DOI]'))
    curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=$term&retmode=json" -o pilot_pmc_s.json
    Pause-Short
    $j = (Get-Content -Raw pilot_pmc_s.json) | ConvertFrom-Json
    $ids = @($j.esearchresult.idlist)
    Write-Output ("$($t.tag) PMC ids: " + ($ids -join ','))
    if ($ids.Count -gt 0) {
        curl.exe -s -m 30 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=$($ids[0])&retmode=xml" -o ("pilot_ground_" + $t.tag + "_pmc.xml")
        Write-Output ("  saved " + (Get-Item ("pilot_ground_" + $t.tag + "_pmc.xml")).Length + " bytes")
        Pause-Short
    }
}

# --- Semantic Scholar abstracts ---
$ss = @(
    @{ tag='grillner1975'; doi='10.1016/0006-8993(75)90401-1' },
    @{ tag='mccrea1980'; doi='10.1152/jn.1980.44.3.475' },
    @{ tag='pearson1995'; doi='10.1016/0959-4388(95)80107-3' }
)
foreach ($s in $ss) {
    $enc = [uri]::EscapeDataString(('DOI:' + $s.doi))
    curl.exe -s -m 20 "https://api.semanticscholar.org/graph/v1/paper/$enc?fields=title,abstract,year,venue" -o ("pilot_ss_" + $s.tag + ".json")
    Start-Sleep -Milliseconds 1300
}
Get-ChildItem pilot_ss_*.json | ForEach-Object {
    Write-Output ("=== " + $_.Name + " ===")
    $raw = Get-Content -Raw $_.FullName
    Write-Output $raw.Substring(0, [Math]::Min(1600, $raw.Length))
}
