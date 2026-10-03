$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco'
$py = 'D:\Anaconda\envs\myo\python.exe'
$cmd = "cmd /c cd /d $wd && $py test_li_stepping.py >> li_drop.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
Start-Sleep -Seconds 2
$p = Get-Process python -ErrorAction SilentlyContinue | Sort-Object StartTime -Descending | Select-Object -First 1
if ($p) { $p.PriorityClass = 'BelowNormal' }
"li drop-protocol gate launched"
