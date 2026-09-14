# Search archive.org for J Physiol vol 40 (1910)
curl.exe -s -m 30 "https://archive.org/advancedsearch.php?q=title%3A%28journal+of+physiology%29+AND+volume%3A40&fl%5B%5D=identifier&fl%5B%5D=title&rows=20&output=json" -o pilot_ia_search.json
$j = (Get-Content -Raw pilot_ia_search.json) | ConvertFrom-Json
$docs = $j.response.docs
Write-Output ("IA hits: " + $docs.Count)
foreach ($d in $docs) { Write-Output ("  " + $d.identifier + " | " + $d.title) }
