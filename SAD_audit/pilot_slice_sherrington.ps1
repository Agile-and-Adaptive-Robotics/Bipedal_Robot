# Slice Sherrington 1910 (J Physiol 40:28-121) out of the whole-volume OCR text.
$raw = Get-Content -Raw -Encoding UTF8 pilot_sherrington1910_ia.txt
$i = $raw.IndexOf('FLEXION-REFLEX')
if ($i -lt 0) { $i = $raw.IndexOf('Flexion-reflex') }
Write-Output "start index: $i"
if ($i -ge 0) {
    # take ~60k chars from the article start (covers intro + methods + early results)
    $slice = $raw.Substring($i, [Math]::Min(60000, $raw.Length - $i))
    $slice | Set-Content -Encoding UTF8 pilot_ground_sherrington1910.txt
    Write-Output ("saved " + $slice.Length + " chars")
    # also show the first 1200 chars as a sanity check
    Write-Output $slice.Substring(0, 1200)
}
