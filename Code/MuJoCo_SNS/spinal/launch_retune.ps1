$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
$py = 'D:\Anaconda\envs\myo\python.exe'
$jobs = @(
  @{ n = 'pruned_air';  c = '_curriculum_pruned.py air' },
  @{ n = 'pruned_w1';   c = '_curriculum_pruned.py walk 1' },
  @{ n = 'pruned_w2';   c = '_curriculum_pruned.py walk 2' },
  @{ n = 'pruned_w3';   c = '_curriculum_pruned.py walk 3' }
)
foreach ($j in $jobs) {
  $cmd = "cmd /c cd /d $wd && $py $($j.c) >> $($j.n).log 2>&1"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
  Start-Sleep -Milliseconds 400
}
"launched: " + ($jobs.n -join ', ')
