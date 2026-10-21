Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
  Where-Object { $_.CommandLine -match '_curriculum_s3p' } |
  ForEach-Object { Write-Host ("killing python " + $_.ProcessId); Stop-Process -Id $_.ProcessId -Force }
Get-CimInstance Win32_Process -Filter "Name='cmd.exe'" |
  Where-Object { $_.CommandLine -match 's3p' -and $_.CommandLine -match 'launch' } |
  ForEach-Object { Write-Host ("killing cmd " + $_.ProcessId); Stop-Process -Id $_.ProcessId -Force }
