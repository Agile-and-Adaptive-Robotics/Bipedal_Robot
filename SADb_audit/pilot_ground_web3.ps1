# Retry PubMed abstracts after rate-limit pause: Grillner 1975 (PMID 1148835)
Start-Sleep -Seconds 3
curl.exe -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=1148835&rettype=abstract&retmode=text" -o pilot_ground_grillner1975.txt
Write-Output ("grillner1975 saved " + (Get-Item pilot_ground_grillner1975.txt).Length + " bytes")
Get-Content pilot_ground_grillner1975.txt
