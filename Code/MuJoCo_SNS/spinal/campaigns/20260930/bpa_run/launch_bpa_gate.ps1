$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\campaigns\20260930\bpa_run'
$py  = 'D:\Anaconda\envs\myo\python.exe'
$pp  = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal;D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS'
$cmd = "cmd /c cd /d $wd && set PYTHONPATH=$pp&& $py test_bpa_stepping.py --dur=20 --hold=6 --lower=2 --kmax-frac=0.6 --cap=0.6 --count=2 >> bpa_gate.log 2>&1"
$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{CommandLine=$cmd; CurrentDirectory=$wd}
Write-Output ("pid=" + $r.ProcessId + " ret=" + $r.ReturnValue)
