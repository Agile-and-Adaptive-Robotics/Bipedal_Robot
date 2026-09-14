# 1) APS free first page PNG for McCrea 1980 (abstract is in the image)
curl.exe -s -m 40 -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64)" "https://journals.physiology.org/na101/home/literatum/publisher/physio/journals/content/jn/1980/jn.1980.44.issue-3/jn.1980.44.3.475/production/jn.1980.44.3.475.fp.png_v03" -o pilot_fp_mccrea1980.png
Write-Output ("mccrea fp.png bytes: " + (Get-Item pilot_fp_mccrea1980.png).Length)

# 2) PMC scanned PDF for Sherrington 1910
curl.exe -s -L -m 90 -A "Mozilla/5.0 (Windows NT 10.0; Win64; x64)" "https://pmc.ncbi.nlm.nih.gov/articles/PMC1533734/pdf/jphysiol02243-0033.pdf" -o pilot_sherrington1910_pmc.pdf
Write-Output ("sherrington pmc.pdf bytes: " + (Get-Item pilot_sherrington1910_pmc.pdf).Length)
