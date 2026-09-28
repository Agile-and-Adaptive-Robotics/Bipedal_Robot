$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
Remove-Item (Join-Path $wd 'prune_*.log') -ErrorAction SilentlyContinue
Remove-Item (Join-Path $wd 'prune_*.err') -ErrorAction SilentlyContinue
foreach ($v in 's3k','w2lvar','syn6') {
  foreach ($s in 0,1,2) {
    $cmd = "cmd /c cd /d $wd && D:\Anaconda\envs\myo\python.exe prune_matrix.py $v $s 3 >> prune_$v`_$s.log 2>&1"
    $r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd }
    if ($r.ReturnValue -ne 0) { Write-Host ("FAIL {0} {1}: {2}" -f $v, $s, $r.ReturnValue) }
    Start-Sleep -Milliseconds 300
  }
}
"launched via WMI"
