$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco'
$py = 'D:\Anaconda\envs\myo\python.exe'
$cmd = "cmd /c cd /d $wd && $py test_w2l_air.py > w2l_air_rerun.log 2>&1 && $py test_w2l_ground.py > w2l_ground_protocol.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
Start-Sleep -Seconds 15
Get-Process python -ErrorAction SilentlyContinue | ForEach-Object { try { $_.PriorityClass = 'BelowNormal' } catch {} }
"W2L gate chain launched (air rerun -> ground protocol)"
