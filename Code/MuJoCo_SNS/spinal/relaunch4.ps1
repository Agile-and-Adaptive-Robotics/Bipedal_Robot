$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
$py = 'D:\Anaconda\envs\myo\python.exe'
foreach ($j in @(@('air','pruned_air'), @('walk 1','pruned_w1'), @('walk 2','pruned_w2'), @('walk 3','pruned_w3'))) {
  $cmd = "cmd /c cd /d $wd && $py _curriculum_pruned.py $($j[0]) >> $($j[1]).log 2>&1"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
  Start-Sleep -Milliseconds 600
}
Start-Sleep -Seconds 90
Get-Process python -ErrorAction SilentlyContinue | ForEach-Object { try { $_.PriorityClass = 'BelowNormal' } catch {} }
"retune relaunched (3rd): " + (Get-Process python -ErrorAction SilentlyContinue | Measure-Object).Count + " python procs"
