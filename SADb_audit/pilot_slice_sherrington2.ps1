# Slice Section III-IV of Sherrington 1910 (journal pp. ~55-80) from volume OCR.
$raw = Get-Content -Raw -Encoding UTF8 pilot_sherrington1910_ia.txt
# article started at 78781; ~3750 OCR chars/page; p.55 is ~27 pages in
$start = 78781 + 100000
$slice = $raw.Substring($start, [Math]::Min(90000, $raw.Length - $start))
$slice | Set-Content -Encoding UTF8 pilot_ground_sherrington1910_s4.txt
Write-Output ("saved from " + $start + ", " + $slice.Length + " chars")
Write-Output ("FIRST 400: " + ($slice.Substring(0,400) -replace '\s+', ' '))
foreach ($kw in @('flexion-phase','extension-phase','umkehr','reflex walking','extensor-thrust','nocip')) {
    $idx = $slice.ToLower().IndexOf($kw.ToLower())
    Write-Output ("kw [$kw]: " + $idx)
}
