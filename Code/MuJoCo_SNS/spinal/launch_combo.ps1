$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
$cmd = "cmd /c cd /d $wd && D:\Anaconda\envs\myo\python.exe prune_combo.py >> prune_combo.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
"combo launched"
