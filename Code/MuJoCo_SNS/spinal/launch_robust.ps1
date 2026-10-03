$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
$py = 'D:\Anaconda\envs\myo\python.exe'
foreach ($v in 's3k','s3kpruned','w2lvar','syn6') {
  foreach ($s in 0,1) {
    $cmd = "cmd /c cd /d $wd && $py prune_robust.py $v $s 2 >> robust_${v}_$s.log 2>&1"
    Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
    Start-Sleep -Milliseconds 400
  }
}
"robust sweep launched (8 shards)"
