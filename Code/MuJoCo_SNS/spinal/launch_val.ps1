$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
$cmd = "cmd /c cd /d $wd && D:\Anaconda\envs\myo\python.exe gait_validate_all.py >> gait_validate.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
"validation launched"
