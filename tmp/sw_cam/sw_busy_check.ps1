$p = Get-Process SLDWORKS
$t1 = $p.TotalProcessorTime.TotalSeconds
$r1 = $p.Responding
Start-Sleep -Seconds 8
$p.Refresh()
$t2 = $p.TotalProcessorTime.TotalSeconds
Write-Output ("CPU seconds delta over 8s: " + [math]::Round($t2 - $t1, 2))
Write-Output ("Responding: " + $p.Responding + " (was " + $r1 + ")")
