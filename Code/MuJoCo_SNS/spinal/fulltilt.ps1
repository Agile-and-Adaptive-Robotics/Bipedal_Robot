$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
$py = 'D:\Anaconda\envs\myo\python.exe'
# validation (crash-safe: skips already-scored variants)
$cmd = "cmd /c cd /d $wd && $py gait_validate_all.py >> gait_validate.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
# robustness WALK rerun (skips done STAND/PUSH cells)
foreach ($v in 's3k','s3kpruned','w2lvar','syn6') {
  $cmd = "cmd /c cd /d $wd && $py prune_robust.py $v 0 1 >> robust_${v}_walkfix.log 2>&1"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
  Start-Sleep -Milliseconds 300
}
Start-Sleep -Seconds 3
# full tilt: nobody else on the box - restore Normal priority
Get-Process python -ErrorAction SilentlyContinue | ForEach-Object { try { $_.PriorityClass = 'Normal' } catch {} }
"full tilt: " + (Get-Process python -ErrorAction SilentlyContinue | Measure-Object).Count + " python procs at Normal"
