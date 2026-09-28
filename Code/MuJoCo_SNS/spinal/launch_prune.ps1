$py = 'D:\Anaconda\envs\myo\python.exe'
$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
foreach ($v in 's3k','w2lvar','syn6') {
  foreach ($s in 0,1,2) {
    $out = Join-Path $wd ("prune_{0}_{1}.log" -f $v,$s)
    $err = Join-Path $wd ("prune_{0}_{1}.err" -f $v,$s)
    Start-Process -FilePath $py -ArgumentList @('prune_matrix.py',$v,"$s",'3') -WorkingDirectory $wd -WindowStyle Hidden -RedirectStandardOutput $out -RedirectStandardError $err
    Start-Sleep -Milliseconds 300
  }
}
Start-Sleep -Seconds 5
$n = (Get-Process python -ErrorAction SilentlyContinue | Measure-Object).Count
"python processes running: $n"
