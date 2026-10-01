$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
# kill only the s3k shards (match on command line)
$procs = Get-CimInstance Win32_Process -Filter "Name = 'python.exe'" |
  Where-Object { $_.CommandLine -match 'prune_matrix\.py s3k' }
foreach ($p in $procs) { Stop-Process -Id $p.ProcessId -Force }
Start-Sleep -Seconds 2
"killed: $($procs.Count)"
Remove-Item (Join-Path $wd 'prune_results_s3k.jsonl') -ErrorAction SilentlyContinue
Remove-Item (Join-Path $wd 'prune_s3k_*.log') -ErrorAction SilentlyContinue
foreach ($s in 0,1,2) {
  $cmd = "cmd /c cd /d $wd && D:\Anaconda\envs\myo\python.exe prune_matrix.py s3k $s 3 >> prune_s3k_$s.log 2>&1"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
  Start-Sleep -Milliseconds 300
}
"relaunched 3 s3k shards"
