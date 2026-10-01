$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco'
$py = 'D:\Anaconda\envs\myo\python.exe'
$jobs = @(
  @{k='support085'; a='--support=0.85'},
  @{k='support05';  a='--support=0.5'},
  @{k='support100'; a='--support=1.0'},
  @{k='clear002';   a='--clear=0.02'},
  @{k='clear006';   a='--clear=0.06'},
  @{k='hold20';     a='--hold=2.0'}
)
foreach ($j in $jobs) {
  $cmd = "cmd /c cd /d $wd && $py test_li_stepping.py $($j.a) >> li_sw_$($j.k).log 2>&1"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
  Start-Sleep -Milliseconds 400
}
Start-Sleep -Seconds 20
Get-Process python -ErrorAction SilentlyContinue | ForEach-Object { try { $_.PriorityClass = 'BelowNormal' } catch {} }
"sweep launched: 6 jobs"
