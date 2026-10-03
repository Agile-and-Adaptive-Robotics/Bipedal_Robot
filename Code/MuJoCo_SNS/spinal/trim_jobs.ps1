$procs = Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
  Where-Object { $_.CommandLine -match 'gait_validate|prune_robust' }
foreach ($p in $procs) { Stop-Process -Id $p.ProcessId -Force; Write-Host ("stopped pid " + $p.ProcessId) }
Write-Host "--- remaining python jobs ---"
Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
  ForEach-Object { Write-Host $_.CommandLine }
