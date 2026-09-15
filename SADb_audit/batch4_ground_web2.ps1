# Batch 4 grounding round 3: Europe PMC core for Armstrong/Bassler PMIDs;
# archive.org search for J Physiol vol 137 (1957).
$ErrorActionPreference = 'Continue'
foreach ($pmid in @('3526411','9639677')) {
    curl.exe -s -m 20 "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:$pmid&resultType=core&format=json" -o ("batch4_epmc_$pmid.json")
    Start-Sleep -Milliseconds 500
    $j = (Get-Content -Raw ("batch4_epmc_$pmid.json")) | ConvertFrom-Json
    $r = $j.resultList.result[0]
    Write-Output "=== PMID $pmid ==="
    Write-Output ("title: " + $r.title)
    $ab = $r.abstractText
    if ($ab) { $ab | Set-Content -Encoding UTF8 ("batch4_epmc_abs_$pmid.txt"); Write-Output ("abstract saved: " + $ab.Length + " chars") } else { Write-Output "abstract: NONE" }
}
# archive.org: J Physiol 1957 vol 137
curl.exe -s -m 30 "https://archive.org/advancedsearch.php?q=title%3A%28journal+of+physiology%29+AND+volume%3A137&fl%5B%5D=identifier&fl%5B%5D=title&rows=10&output=json" -o batch4_ia.json
$j2 = (Get-Content -Raw batch4_ia.json) | ConvertFrom-Json
foreach ($d in $j2.response.docs) { Write-Output ("IA: " + $d.identifier + " | " + $d.title) }
