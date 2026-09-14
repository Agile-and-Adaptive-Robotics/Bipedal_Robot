# Europe PMC core records (abstractText) for the two PMIDs.
$ErrorActionPreference = 'Continue'
foreach ($pmid in @('1148835','7441311')) {
    $ok = $false
    for ($i = 1; $i -le 4 -and -not $ok; $i++) {
        curl.exe -s -m 20 "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:$pmid&resultType=core&format=json" -o ("pilot_epmc_$pmid.json")
        $raw = Get-Content -Raw ("pilot_epmc_$pmid.json")
        if ($raw -and $raw.Contains('resultList')) {
            $ok = $true
            $j = $raw | ConvertFrom-Json
            $r = $j.resultList.result[0]
            Write-Output "=== PMID $pmid ==="
            Write-Output ("title: " + $r.title)
            Write-Output ("journal: " + $r.journalInfo.journal.title + " " + $r.pubYear)
            $ab = $r.abstractText
            if ($ab) { Write-Output ("abstract: " + $ab) } else { Write-Output "abstract: NONE" }
        } else {
            Write-Output "PMID $pmid try $i failed: $($raw.Substring(0,[Math]::Min(120,$raw.Length)))"
            Start-Sleep -Seconds 3
        }
    }
}
